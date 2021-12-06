#include "clients/time_dependent.hpp"

#include "operator/diagnostics.hpp"
#include "time/rk3_tvd.hpp"
#include "tools/ioh5.hpp"

using std::string;
using std::to_string;

TimeDependent::~TimeDependent() {
    //-------------------------------------------------------------------------
    // delete the field
    m_profStart(prof_, "cleanup");
    delete grid_;
    m_profStop(prof_, "cleanup");

    if (!(prof_ == nullptr)) {
        prof_->Disp();
        delete prof_;
    }
    //-------------------------------------------------------------------------
}

void TimeDependent::InitParam(ParserArguments* param) {
    //--------------------------------------------------------------------------
    // call the general testcase parameters
    this->TestCase::InitParam(param);

    // the the standard stuffs
    if (param->profile) {
        int comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        string name = this->name() + to_string(comm_size) + string("ranks") + string("_w") + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT);
        prof_       = new Prof(name);
    }
    m_profStart(prof_, "init");

    // setup the grid
    grid_ = new Grid(param->init_lvl, param->period, param->length, M_GRIDBLOCK, MPI_COMM_WORLD, prof_);
    m_log("grid created at level %d - size = %dx%dx%d - periodic? %dx%dx%d", param->init_lvl,
          param->length[0], param->length[1], param->length[2],
          param->period[0], param->period[1], param->period[2]);

    // set the min/max level
    grid_->level_limit(param->level_min, param->level_max);
    no_adapt_ = param->no_adapt;

    // tolerances
    optimal_tol_       = param->optimal_tol;
    real_t coarsen_tol = (param->coarsen_tol < 0.0) ? 0.0 : ((optimal_tol_) ? (param->refine_tol / pow(2.0, M_WAVELET_N)) : (param->coarsen_tol));
    grid_->SetTol(param->refine_tol, coarsen_tol);

    // call the setup function
    Setup(param);

    m_profStop(prof_, "init");
    //--------------------------------------------------------------------------
}

void TimeDependent::Run() {
    m_begin;
    //--------------------------------------------------------------------------
    // time integration
    iter_t iter = 0;
    real_t t    = tstart_;
    real_t dt   = 0.0;

    //--------------------------------------------------------------------------
    m_profStart(prof_, "diagnostics");
    {
        m_log_level_plus;
        real_t time_now = MPI_Wtime();
        Diagnostics(t, dt, iter, time_now);
        m_log_level_minus;
    }

    // let's gooo
    m_profStart(prof_, "run");
    const real_t wtime_start = MPI_Wtime();
    while (t < tfinal_ && iter < iter_max()) {
        m_log("--------------------------------------------------------------------------------");
        //......................................................................
        m_log("---- do time-step");
        m_profStart(prof_, "do dt");
        {
            m_log_level_plus;
            DoTimeStep(&t, &dt);
            iter++;

            // dump some info
            real_t wtime_now = MPI_Wtime();
            m_log("now -> time = %f/%f - step %d/%d - dt used = %e - wtime = %e min", t, tfinal_, iter, iter_max(), dt, (wtime_now - wtime_start) / 60.0);

            m_log_level_minus;
        }
        m_profStop(prof_, "do dt");

        //......................................................................
        // adapt the mesh
        if (iter % iter_adapt() == 0) {
            if (!no_adapt_) {
                m_log("---- adapt mesh");
                m_profStart(prof_, "adapt");
                {
                    m_log_level_plus;
                    Adapt(t, dt);
                    m_log_level_minus;
                }
                m_profStop(prof_, "adapt");
            }
        }

        //......................................................................
        // diagnostics, dumps, whatever
        if (iter % iter_diag() == 0) {
            m_log("---- run diag");
            m_profStart(prof_, "diagnostics");
            {
                m_log_level_plus;
                real_t time_now = MPI_Wtime();
                Diagnostics(t, dt, iter, time_now - wtime_start);
                m_log_level_minus;
            }
            m_profStop(prof_, "diagnostics");
        }
    }
    m_profStop(prof_, "run");
    //--------------------------------------------------------------------------
    // run the last diag
    m_profStart(prof_, "diagnostics");
    {
        m_log_level_plus;
        real_t time_now = MPI_Wtime();
        Diagnostics(t, dt, iter, time_now);
        m_log_level_minus;
    }
    //--------------------------------------------------------------------------
    m_end;
}
