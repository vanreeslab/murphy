#include "prof.hpp"

#include <mpi.h>

using std::map;
using std::string;

static constexpr rank_t    upper_rank = 1000; // approximates the infinity of procs
static map<rank_t, real_t> t_nu       = {{0, 0.0},
                                   {1, 6.314},
                                   {2, 2.920},
                                   {3, 2.353},
                                   {4, 2.132},
                                   {5, 2.015},
                                   {7, 1.895},
                                   {10, 1.812},
                                   {15, 1.753},
                                   {20, 1.725},
                                   {30, 1.697},
                                   {50, 1.676},
                                   {100, 1.660},
                                   {upper_rank, 1.645}};
/**
 * @brief return the t_nu for 90% confidence interval width based on the interpolation of the above table
 * 
 * @param nu the number of proc-1
 * @return real_t the confidence interval param
 */
real_t t_nu_interp(const rank_t nu) {
    m_assert(nu >= 0, "the nu param = %d must be positive", nu);
    //-------------------------------------------------------------------------
    if (nu == 0) {
        // easy, it's 0
        return 0.0;
    } else if (nu <= 5) {
        // we have an exact entry
        const auto it = t_nu.find(nu);
        return it->second;
    } else if (nu >= upper_rank) {
        // we are too big, it's like a normal distribution
        return 1.645;
    } else {
        // find the right point
        auto         it_up  = t_nu.lower_bound(nu);  // first element >= nu
        const rank_t nu_up  = it_up->first;
        const real_t t_up   = it_up->second;
        auto         it_low = std::prev(it_up,1);  // take the previous one
        const rank_t nu_low = it_low->first;
        const real_t t_low  = it_low->second;
        // m_log("nu_up = %d, nu_low = %d, t_up=%f, t_low=%f",nu_up,nu_low,t_up, t_low);
        return t_low + (t_up-t_low)/(nu_up-nu_low)*(nu-nu_low);
    }
    //-------------------------------------------------------------------------
}

TimerBlock::TimerBlock(string name) {
    name_  = name;
    t0_    = -1.0;
    t1_    = -1.0;
    count_ = 0;
}

TimerBlock::~TimerBlock() {
    for (auto it = children_.begin(); it != children_.end(); ++it) {
        delete it->second;
    }
}

/**
 * @brief start the timer
 * 
 */
void TimerBlock::Start() {
    m_assert(t0_ < -0.5, "the block %s has already been started", name_.c_str());
    count_ += 1;
    t0_ = MPI_Wtime();
}

/**
 * @brief stop the timer
 * 
 */
void TimerBlock::Stop() {
    m_assert(t0_ > -0.5, "the block %s is stopped without being started", name_.c_str());
    // get the time
    t1_ = MPI_Wtime();
    // store it
    real_t dt = t1_ - t0_;
    time_acc_ = time_acc_ + dt;
    // reset to negative for the checks
    t0_ = -1.0;
    t1_ = -1.0;
}

TimerBlock* TimerBlock::AddChild(string child_name) noexcept {
    // find cleanly the child
    auto it = children_.find(child_name);

    if (it == children_.end()) {
        TimerBlock* child     = new TimerBlock(child_name);
        children_[child_name] = child;
        child->SetParent(this);
        return child;
    } else {
        return it->second;
    }
}

/**
 * @brief store the parent pointer
 * 
 * @param parent 
 */
void TimerBlock::SetParent(TimerBlock* parent) {
    parent_ = parent;
    // is_root_ = false;
}

/**
 * @brief display the time accumulated
 * 
 * If the timer is a ghost timer (never called but still created), we return the sum on the children
 * 
 * @return real_t 
 */
real_t TimerBlock::time_acc() const {
    if (count_ > 0) {
        return time_acc_;
    } else {
        real_t sum = 0.0;
        for (auto it = children_.cbegin(); it != children_.cend(); ++it) {
            const TimerBlock* child = it->second;
            sum += child->time_acc();
        }
        return sum;
    }
}

/**
 * @brief display the time for the TimerBlock
 * 
 * @param file 
 * @param level 
 * @param total_time 
 */
void TimerBlock::Disp(FILE* file, const level_t level, const real_t total_time, const lda_t icol) const {
    // check if any proc has called the agent
    int total_count = 0;
    MPI_Allreduce(&count_, &total_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // get the size and useful stuffs
    int comm_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

// setup the displayed name
#ifdef COLOR_PROF
    string shifter = "\033[0;35m";
#else
    string shifter;
#endif
    for (int l = 1; l < level - 1; l++) {
        shifter = shifter + "|   ";
    }
    if (level > 1) {
        shifter = shifter + "|-> ";
    }
// string myname = shifter + "\033[0m" + "\033[1m" + name_ + "\033[0m";
#ifdef COLOR_PROF
    string myname = shifter + "\033[0m" + name_;
#else
    string myname = shifter + name_;
#endif

    //................................................
    // compute my numbers
    if (total_count > 0) {
        real_t scale = 1.0 / comm_size;

        // compute the counters (mean, max, min)
        real_t local_count = count_;
        real_t max_count;
        MPI_Allreduce(&local_count, &max_count, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // compute times passed inside + children
        real_t local_time = time_acc_;
        real_t sum_time;
        MPI_Allreduce(&local_time, &sum_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        real_t mean_time           = sum_time / comm_size;
        real_t mean_time_per_count = sum_time / total_count;
        real_t glob_percent        = mean_time / total_time * 100.0;

        // confidence interval 90% using the t distribution
        real_t sum_timesq;
        real_t local_timesq = (local_time - mean_time) * (local_time - mean_time);
        MPI_Allreduce(&local_timesq, &sum_timesq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        real_t std_time   = sqrt(sum_timesq / (comm_size - 1));
        real_t ci_90_time = std_time / sqrt(comm_size) * t_nu_interp(comm_size - 1);

        // printf the important information
        if (rank == 0) {
#ifdef COLOR_PROF
            // printf("%-25.25s|  %9.4f\t%9.4f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%09.1f\t%9.2f\n", myname.c_str(), glob_percent, loc_percent, mean_time, self_time, mean_time_per_count, min_time_per_count, max_time_per_count, mean_count, mean_bandwidth);
            if (icol == 0) {  // go red
                printf("%-60.60s %s\033[0;31m%09.6f\033[0m %% -> \033[0;31m%07.4f\033[0m [s] +- %07.4f [s] \t\t\t(%.4f [s/call], %.0f calls)\n", myname.c_str(), shifter.c_str(), glob_percent, mean_time, ci_90_time, mean_time_per_count, max_count);
            }
            if (icol == 1) {  // go orange
                printf("%-60.60s %s\033[0;33m%09.6f\033[0m %% -> \033[0m%07.4f\033[0m [s] +- %07.4f [s] \t\t\t(%.4f [s/call], %.0f calls)\n", myname.c_str(), shifter.c_str(), glob_percent, mean_time, ci_90_time, mean_time_per_count, max_count);
            }
            if (icol == 2) {  // go normal
                printf("%-60.60s %s\033[0m%09.6f\033[0m %% -> \033[0m%07.4f\033[0m [s] +- %07.4f [s] \t\t\t(%.4f [s/call], %.0f calls)\n", myname.c_str(), shifter.c_str(), glob_percent, mean_time, ci_90_time, mean_time_per_count, max_count);
            }
#else
            printf("%-60.60s %s%09.6f %% -> %07.4f [s] +- %07.4f [s] \t\t\t(%.4f [s/call], %.0f calls)\n", myname.c_str(), shifter.c_str(), glob_percent, mean_time, ci_90_time, mean_time_per_count, max_count);
#endif
            // printf in the file
            if (file != nullptr) {
                fprintf(file, "%s;%d;%.6f;%.6f;%.6f;%.0f\n", name_.c_str(), level, mean_time, glob_percent, mean_time_per_count, max_count);
            }
        }
    } else if (name_ != "root") {
        // printf the important information
        if (rank == 0) {
            //     printf("%-35.35s| %s\033[0m%09.6f\033[0m %% -> \033[0m%07.4f\033[0m [s] \t\t(mean/call %.4f [s], %.0f calls)\n", myname.c_str(), shifter.c_str(), 0.0, 0.0, 0.0, 0.0);
            if (file != nullptr) {
                fprintf(file, "%s\n", name_.c_str());
            }
        }
    }

    //................................................
    // check that everything is ok for the MPI
#ifndef NDEBUG
    iblock_t nchildren     = children_.size();
    iblock_t nchildren_max = 0;
    iblock_t nchildren_min = 0;
    MPI_Allreduce(&nchildren, &nchildren_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&nchildren, &nchildren_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    m_assert((nchildren_max == nchildren) && (nchildren == nchildren_min), "TimerBlock %s: nchildren do not match: local = %d, max = %d, min = %d", name_.c_str(), nchildren, nchildren_max, nchildren_min);
#endif

    //................................................
    // recursive call to the childrens
    // get the colors
    string max_name, min_name;
    real_t max_time = -1e+15;
    real_t min_time = +1e+15;
    for (auto it = children_.cbegin(); it != children_.cend(); ++it) {
        TimerBlock* child = it->second;
        real_t      ctime = child->time_acc();
        if (ctime > max_time) {
            max_time = ctime;
            max_name = child->name();
        }
        if (ctime < min_time) {
            min_time = ctime;
            min_name = child->name();
        }
    }
    for (auto it = children_.cbegin(); it != children_.cend(); ++it) {
        TimerBlock* child = it->second;
        if (child->name() == max_name && icol == 0) {
            // go red
            child->Disp(file, level + 1, total_time, 0);
        } else if (child->name() == max_name && icol > 0) {
            // go orange
            child->Disp(file, level + 1, total_time, 1);
        } else {
            child->Disp(file, level + 1, total_time, 2);
        }
    }
}

//===============================================================================================================================
//===============================================================================================================================
//===============================================================================================================================

/**
 * @brief Construct a new Prof with a default name
 */
Prof::Prof() : name_("default") {
    // create the current Timer Block "root"
    current_ = new TimerBlock("root");
}

/**
 * @brief Construct a new Prof with a given name
 */
Prof::Prof(const string myname) : name_(myname) {
    current_ = new TimerBlock("root");
}

/**
 * @brief Destroy the Prof
 */
Prof::~Prof() {
    delete current_;
}

/**
 * @brief initialize the timer and move to it
 */
void Prof::Init(string name) {
    current_ = current_->AddChild(name);
}

/**
 * @brief start the timer of the TimerBlock
 */
void Prof::Start(string name) {
    current_->Start();
}

/**
 * @brief stop the timer of the TimerBlock
 */
void Prof::Stop(string name) {
    m_assert(name == current_->name(), "we are trying to stop %s which is not the most recent timer started = %s", name.c_str(), current_->name().c_str());
    current_->Stop();
}

/**
 * @brief go back to the parent of the present timer block
 */
void Prof::Leave(string name) {
    current_ = current_->parent();
}

/**
 * @brief display the whole profiler using 
 * 
 */
void Prof::Disp() const {
    m_assert(current_->name() == "root", "the current TimerBlock is not the root, please stop any current timer firsts");
    int comm_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    FILE*  file;
    string folder = "./prof";

    MPI_Barrier(MPI_COMM_WORLD);

    // //-------------------------------------------------------------------------
    // /** - I/O of the parentality */
    // //-------------------------------------------------------------------------

    // if (rank == 0) {
    //     struct stat st = {0};
    //     if (stat(folder.c_str(), &st) == -1) {
    //         mkdir(folder.c_str(), 0770);
    //     }

    //     string filename = folder + "/" + name_ + "_parent.csv";
    //     file            = fopen(filename.c_str(), "w+");

    //     if (file != nullptr) {
    //         // current_->DumpParentality(file, 0);
    //         fclose(file);
    //     } else {
    //         printf("unable to open file %s !", filename.c_str());
    //     }
    // }

    //-------------------------------------------------------------------------
    /** - do the IO of the timing */
    //-------------------------------------------------------------------------

    if (rank == 0) {
        string filename = "./prof/" + name_ + "_time.csv";
        file            = fopen(filename.c_str(), "w+");
    }

    // get the global timing
    real_t total_time = current_->time_acc();

    // display the header
    if (rank == 0) {
        printf("===================================================================================================================================================\n");
#ifdef COLOR_PROF
        printf("        PROFILER %s --> total time = \033[0;33m%.4f\033[m [s] \n\n", name_.c_str(), total_time);
#else
        printf("        PROFILER %s --> total time = %.4f [s] \n\n", name_.c_str(), total_time);
#endif
    }

    // display root with the total time, root is the only block which is common to everybody
    current_->Disp(file, 0, total_time, 0);

    // display footer
    if (rank == 0) {
#ifdef COLOR_PROF
        printf("===================================================================================================================================================\n");
        printf("WARNING:\n");
        printf("  - times are mean-time with their associated 90%% CI\n");
        printf("  - the percentage might not be consistent are they only reflect rank-0 timing\n");
        printf("legend:\n");
        printf("  - \033[0;31mthis indicates the most expensive step of the most expensive operation\033[0m\n");
        printf("  - \033[0;33mthis indicates the most expensive step of the parent operation\033[0m\n");
        printf("===================================================================================================================================================\n");
#endif

        if (file != nullptr) {
            fclose(file);
        } else {
            printf("unable to open file for profiling !");
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
