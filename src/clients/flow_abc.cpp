#include "flow_abc.hpp"

#include "core/forloop.hpp"
#include "operator/advection.hpp"
#include "operator/xblas.hpp"
#include "time/rk3_tvd.hpp"
#include "tools/ioh5.hpp"

using std::string;
using std::to_string;

FlowABC::~FlowABC() {
    //-------------------------------------------------------------------------
    delete (vel_);
    delete (scal_);
    delete (grid_);

    if (prof_ != nullptr) {
        prof_->Disp();
        delete (prof_);
    }

    m_log("Navier Stokes is dead");
    //-------------------------------------------------------------------------
}

void FlowABC::InitParam(ParserArguments* param) {
    //-------------------------------------------------------------------------
    // call the general testcase parameters
    this->TestCase::InitParam(param);

    // the the standard stuffs
    if (param->profile) {
        int comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        string name = string("ABC_Flow_") + to_string(comm_size) + string("ranks");
        prof_       = new Prof(name);
    }

    // setup the grid
    bool period[3]    = {true, true, true};
    grid_             = new Grid(param->init_lvl, period, param->length, MPI_COMM_WORLD, prof_);
    const real_t L[3] = {1.0, 1.0, 1.0};

    // set the levels
    grid_->level_limit(param->level_min, param->level_max);

    // get the fields
    vel_  = new Field("velocity", 3);
    scal_ = new Field("scalar", 1);
    grid_->AddField(vel_);
    grid_->AddField(scal_);

    // temp the velocity as we reimpose it's value right after adapt
    vel_->is_temp(true);

    // setup the flow ring
    lambda_setvalue_t lambda_init = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
        real_t pos[3];
        block->pos(i0, i1, i2, pos);
        block->data(fid, 0).Write(i0, i1, i2)[0] = scalar_ring(pos, param->vr_center, param->vr_radius, param->vr_sigma, param->vr_normal);
    };
    const bidx_t ghost_len[2] = {grid_->interp()->nghost_front(), grid_->interp()->nghost_back()};
    SetValue     flow_ring(lambda_init, ghost_len);
    flow_ring.Profile(prof_);
    flow_ring(grid_, scal_);

    grid_->SetTol(param->refine_tol, param->coarsen_tol);
    grid_->SetRecursiveAdapt(true);
    grid_->Adapt(scal_, &flow_ring);
    //-------------------------------------------------------------------------
}

void FlowABC::Run() {
    m_begin;
    //-------------------------------------------------------------------------
    // time
    lid_t  iter    = 0;
    real_t t_start = 0.0;
    real_t t_final = 1000.0;
    real_t t       = 0.0;

    // SetABSVelocity flow_vel(1.0, 0.5, 0.25, grid_->interp());
    // setup the flow ring
    const real_t      a          = 1.0;
    const real_t b = 0.5;
    const real_t c = 0.25;
    lambda_setvalue_t lambda_abs = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
        real_t pos[3];
        block->pos(i0, i1, i2, pos);
        const real_t x                           = pos[0];
        const real_t y                           = pos[1];
        const real_t z                           = pos[2];
        block->data(fid, 0).Write(i0, i1, i2)[0] = a * sin(2.0 * M_PI * z) + c * cos(2.0 * M_PI * y);
        block->data(fid, 1).Write(i0, i1, i2)[0] = b * sin(2.0 * M_PI * x) + a * cos(2.0 * M_PI * z);
        block->data(fid, 2).Write(i0, i1, i2)[0] = c * sin(2.0 * M_PI * y) + b * cos(2.0 * M_PI * x);
    };
    const bidx_t ghost_len[2] = {grid_->interp()->nghost_front(),grid_->interp()->nghost_back()};
    SetValue flow_vel(lambda_abs, ghost_len);

    Advection<M_WENO_Z, 3> adv(vel_);
    RK3_TVD                rk3(grid_, scal_, &adv, prof_);
    adv.Profile(prof_);

    // let's gooo
    m_profStart(prof_, "run");
    while (t < t_final && iter < iter_max()) {
        m_log("--------------------");
        //................................................
        // adapt the mesh
        if (iter % iter_adapt() == 0) {
            m_log("---- adapt mesh");
            m_profStart(prof_, "adapt");
            grid_->Adapt(scal_);
            m_profStop(prof_, "adapt");

            // reset the velocity
            flow_vel(grid_, vel_);
            m_assert(vel_->ghost_status(ghost_len), "the velocity ghosts must have been computed");
        }
        // we run the first diagnostic
        if (iter == 0) {
            m_profStart(prof_, "diagnostics");
            m_log("---- run diag");
            Diagnostics(t, 0, iter);
            m_profStop(prof_, "diagnostics");
        }

        //................................................
        // get the time-step given the field
        real_t dt = rk3.ComputeDt(&adv, vel_);

        // dump some info
        m_log("RK3 - time = %f - step %d/%d - dt = %e", t, iter, iter_max(), dt);

        //................................................
        // advance in time
        m_log("---- do time-step");
        m_profStart(prof_, "do dt");
        rk3.DoDt(dt, &t);
        iter++;
        m_profStop(prof_, "do dt");

        //................................................
        // diagnostics, dumps, whatever
        if (iter % iter_diag() == 0) {
            m_profStart(prof_, "diagnostics");
            m_log("---- run diag");
            Diagnostics(t, dt, iter);
            m_profStop(prof_, "diagnostics");
        }
    }
    m_profStop(prof_, "run");
    // run the last diag
    if (iter % iter_diag() != 0) {
        Diagnostics(t, 0.0, iter);
    }
    //-------------------------------------------------------------------------
    m_end;
}

void FlowABC::Diagnostics(const real_t time, const real_t dt, const lid_t iter) {
    m_begin;
    //-------------------------------------------------------------------------
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // if the folder does not exist, create it
    struct stat st = {0};
    if (rank == 0 && stat(folder_diag_.c_str(), &st) == -1) {
        mkdir(folder_diag_.c_str(), 0770);
    }

    // get fields diags
    real_t  scal_min[2], scal_max[2];
    BMinMax minmax;
    minmax(grid_, scal_, scal_min, scal_max);

    // open the file
    FILE*   file_error;
    FILE*   file_diag;
    level_t min_level = grid_->MinLevel();
    level_t max_level = grid_->MaxLevel();
    if (rank == 0) {
        file_diag = fopen(string(folder_diag_ + "/diag.data").c_str(), "a+");
        // iter, time, dt, total quad, level min, level max
        fprintf(file_diag, "%6.6d %e %e %ld %d %d", iter, time, dt, grid_->global_num_quadrants(), min_level, max_level);
        for (lda_t ida = 0; ida < scal_->lda(); ida++) {
            fprintf(file_diag, " %e %e", scal_min[ida], scal_max[ida]);
        }
        fprintf(file_diag, "\n");
        fclose(file_diag);
    }

    if (iter % iter_dump() == 0) {
        // dump the vorticity field
        IOH5 dump(folder_diag_);
        grid_->GhostPull(scal_,ghost_len_ioh5);
        // grid_->GhostPull(vel_);
        dump(grid_, scal_, iter);
        // dump(grid_(), vel_(), iter);

        // dump the details
        if (dump_detail()) {
            Field details("detail", 2);
            details.bctype(M_BC_EXTRAP);
            grid_->AddField(&details);
            grid_->StoreDetails(scal_, &details);

            grid_->GhostPull(&details,ghost_len_ioh5);
            // IOH5 dump("data");
            dump(grid_, &details, iter);

            grid_->DeleteField(&details);
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}
