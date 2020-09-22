#include "prof.hpp"

#include <mpi.h>

TimerBlock::TimerBlock(string name) {
    name_ = name;
}

/**
 * @brief start the timer
 * 
 */
void TimerBlock::Start() {
    count_ += 1;
    t0_ = MPI_Wtime();
}

/**
 * @brief stop the timer
 * 
 */
void TimerBlock::Stop() {
    // get the time
    t1_ = MPI_Wtime();
    // store it
    real_t dt = t1_ - t0_;
    time_acc_ = time_acc_ + dt;
    time_max_ = m_max(time_max_, dt);
    time_min_ = m_min(time_min_, dt);
}

/**
 * @brief reset the timer
 * 
 */
void TimerBlock::Reset() {
    t1_       = 0.0;
    t0_       = 0.0;
    time_acc_ = 0.0;
    time_max_ = 0.0;
    time_min_ = 0.0;
}

/**
 * @brief adds memory to the timer to compute bandwith
 * 
 */
void TimerBlock::AddMem(size_t mem) {
    memsize_ += mem;
}

/**
 * @brief add a child to the timer 
 * 
 * @param child 
 */
void TimerBlock::AddChild(TimerBlock* child) {
    string childName = child->name();

    map<string, TimerBlock*>::iterator it = children_.find(childName);
    // if it does not already exist
    if (it == children_.end()) {
        children_[childName] = child;
        child->SetParent(this);
    }
}
/**
 * @brief store the dady pointer
 * 
 * @param daddy 
 */
void TimerBlock::SetParent(TimerBlock* parent) {
    parent_  = parent;
    is_root_ = false;
}

/**
 * @brief display the time accumulated. If it's a ghost timer (no calls), we sum the time of the children
 * 
 * @return real_t 
 */
real_t TimerBlock::time_acc() const {
    if (count_ > 0) {
        return time_acc_;
    } else {
        real_t sum = 0.0;
        for (map<string, TimerBlock*>::const_iterator it = children_.begin(); it != children_.end(); it++) {
            const TimerBlock* child = it->second;
            sum += child->time_acc();
        }
        return sum;
    }
}

/**
 * @brief display the min time among all calls. If it's a ghost timer (no calls), we sum the time of the children
 * 
 * @return real_t 
 */
real_t TimerBlock::time_min() const {
    if (count_ > 0) {
        return time_min_;
    } else {
        real_t sum = 0.0;
        for (map<string, TimerBlock*>::const_iterator it = children_.begin(); it != children_.end(); it++) {
            const TimerBlock* child = it->second;
            sum += child->time_min();
        }
        return sum;
    }
}

/**
 * @brief display the max time among all calls. If it's a ghost timer (no calls), we sum the time of the children
 * 
 * @return real_t 
 */
real_t TimerBlock::time_max() const {
    if (count_ > 0) {
        return time_max_;
    } else {
        real_t sum = 0.0;
        for (map<string, TimerBlock*>::const_iterator it = children_.begin(); it != children_.end(); it++) {
            const TimerBlock* child = it->second;
            sum += child->time_max();
        }
        return sum;
    }
}

void TimerBlock::DumpParentality(FILE* file, const int level) {
    fprintf(file, "%d;%s", level, name_.c_str());
    for (map<string, TimerBlock*>::const_iterator it = children_.begin(); it != children_.end(); it++) {
        const TimerBlock* child     = it->second;
        string            childName = child->name();
        fprintf(file, ";%s", childName.c_str());
    }
    fprintf(file, "\n");

    for (map<string, TimerBlock*>::const_iterator it = children_.begin(); it != children_.end(); it++) {
        it->second->DumpParentality(file, level + 1);
    }
}

/**
 * @brief display the time for the TimerBlock
 * 
 * @param file 
 * @param level 
 * @param total_time 
 */
void TimerBlock::Disp(FILE* file, const int level, const real_t total_time) {
    // check if any proc has called the agent
    int total_count;
    MPI_Allreduce(&count_, &total_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    // get the size and usefull stuffs
    int comm_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // setup the displayed name
    string myname = name_;
    if (level > 1) {
        myname = " " + myname;
    }
    for (int l = 1; l < level; l++) {
        myname = "--" + myname;
    }
    // if someone has every call the agent, display it
    if (total_count > 0) {
        real_t scale = 1.0 / comm_size;

        // compute the counters (mean, max, min)
        real_t local_count = count_;
        real_t mean_count, max_count, min_count;
        MPI_Allreduce(&local_count, &mean_count, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_count, &max_count, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&local_count, &min_count, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        mean_count *= scale;

        // compute times passed inside + children
        real_t local_time = time_acc_;
        real_t mean_time, max_time, min_time;
        MPI_Allreduce(&local_time, &mean_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        mean_time *= scale;

        real_t mean_time_per_count = mean_time / mean_count;
        real_t min_time_per_count  = min_time / mean_count;
        real_t max_time_per_count  = max_time / mean_count;
        real_t glob_percent        = mean_time / total_time * 100.0;

        // compute the self time  = time passed inside - children
        real_t sum_child = 0.0;
        for (map<string, TimerBlock*>::iterator it = children_.begin(); it != children_.end(); it++) {
            TimerBlock* child = it->second;
            sum_child += child->time_acc();
        }
        real_t loc_self_time = (this->time_acc() - sum_child);
        real_t self_time;
        real_t self_percent;
        MPI_Allreduce(&loc_self_time, &self_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        self_time *= scale;
        self_percent = self_time / total_time * 100.0;

        // comnpute the time passed inside the daddy
        real_t loc_percent;
        if (parent_ != nullptr) {
            real_t dad_local_time = parent_->time_acc();
            real_t dad_mean_time;
            MPI_Allreduce(&dad_local_time, &dad_mean_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            dad_mean_time *= scale;
            loc_percent = mean_time / dad_mean_time * 100.0;
        } else {
            loc_percent = 100.0;
        }

        // compute the bandwith
        real_t local_band_time    = time_acc_;
        real_t local_band_memsize = (real_t)memsize_;
        real_t band_memsize, band_time, mean_bandwidth;
        MPI_Allreduce(&local_band_time, &band_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_band_memsize, &band_memsize, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        mean_bandwidth = (band_memsize / band_time) / std::pow(10.0, 6.0);

        // printf the important information
        if (rank == 0) {
            printf("%-25.25s|  %9.4f\t%9.4f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%9.6f\t%09.1f\t%9.2f\n", myname.c_str(), glob_percent, loc_percent, mean_time, self_time, mean_time_per_count, min_time_per_count, max_time_per_count, mean_count, mean_bandwidth);
            if (file != nullptr) {
                fprintf(file, "%s;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.6f;%09.0f;%09.2f\n", name_.c_str(), glob_percent, loc_percent, mean_time, self_time, mean_time_per_count, min_time_per_count, max_time_per_count, mean_count, mean_bandwidth);
            }
        }
    } else if (name_ != "root") {
        // printf the important information
        if (rank == 0) {
            printf("%-25.25s|\n", myname.c_str());
            if (file != nullptr) {
                fprintf(file, "%s\n", name_.c_str());
            }
        }
    }
    // recursive call to the childrens
    for (map<string, TimerBlock*>::iterator it = children_.begin(); it != children_.end(); it++) {
        TimerBlock* child = it->second;
        child->Disp(file, level + 1, total_time);
    }
}

//===============================================================================================================================
//===============================================================================================================================
//===============================================================================================================================

Prof::Prof() : name_("default") {
    CreateSingle_("root");
}
Prof::Prof(const string myname) : name_(myname) {
    CreateSingle_("root");
}
Prof::~Prof() {
    for (map<string, TimerBlock*>::iterator it = time_map_.begin(); it != time_map_.end(); it++) {
        TimerBlock* current = it->second;
        delete (current);
    }
}

/**
 * @brief create a TimerBlock inside the #time_map_
 * 
 * @param name the key of the entry
 */
void Prof::CreateSingle_(string name) {
    map<string, TimerBlock*>::iterator it = time_map_.find(name);
    // if it does not already exist
    // if (it != time_map_.end()) {
    //     // time_map_[name]->Reset();
    // } else 
    if (it == time_map_.end()) {
        time_map_[name] = new TimerBlock(name);
        time_map_[name]->Reset();
    }
}

/**
 * @brief create a new TimerBlock with "root" as parent
 * 
 * @param name the TimerBlock name
 */
void Prof::Create(string name) {
    Create(name, "root");
}

/**
 * @brief create a new TimerBlock 
 * 
 * @param child the new TimerBlock
 * @param daddy the dad of the new TimerBlock if it does not exists, it is created
 */
void Prof::Create(string child, string parent) {
    // create a new guy
    CreateSingle_(child);
    // find the daddy agent in the root
    map<string, TimerBlock*>::iterator it = time_map_.find(parent);
    if (it == time_map_.end()) {
        Create(parent);
    }
    time_map_[parent]->AddChild(time_map_[child]);
}

/**
 * @brief start the timer of the TimerBlock
 * 
 * @param name the TimerBlock name
 */
void Prof::Start(string name) {
#ifdef NDEBUG
    time_map_[name]->Start();
#else
    map<string, TimerBlock*>::iterator it = time_map_.find(name);
    if (it != time_map_.end()) {
        time_map_[name]->Start();
    } else {
        string msg = "timer " + name + " not found";
    }
#endif
}

/**
 * @brief stop the timer of the TimerBlock
 * 
 * @param name the TimerBlock name
 */
void Prof::Stop(string name) {
#ifdef NDEBUG
    time_map_[name]->Stop();
#else
    map<string, TimerBlock*>::iterator it = time_map_.find(name);
    if (it != time_map_.end()) {
        time_map_[name]->Stop();
    } else {
        string msg = "timer " + name + " not found";
    }
#endif
}

void Prof::AddMem(string name, size_t mem) {
#ifdef NDEBUG
    time_map_[name]->AddMem(mem);
#else
    map<string, TimerBlock*>::iterator it = time_map_.find(name);
    if (it != time_map_.end()) {
        time_map_[name]->AddMem(mem);
    } else {
        string msg = "timer " + name + " not found";
    }
#endif
}

/**
 * @brief get the accumulated time
 * 
 * @param name 
 * @return real_t 
 */
real_t Prof::Time(const std::string ref) {
    int    comm_size;
    real_t local_total_time = time_map_[ref]->time_acc();
    real_t total_time;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    MPI_Allreduce(&local_total_time, &total_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    total_time /= comm_size;

    return total_time;
}

/**
 * @brief display the whole profiler using 
 * 
 */
void Prof::Disp() {
    this->Disp("root");
}
/**
 * @brief display the profiler using the timer spent in ref as a reference for the global percentage computation
 * 
 * @param ref 
 */
void Prof::Disp(const std::string ref) {
    int comm_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //-------------------------------------------------------------------------
    /** - I/O of the parentality */
    //-------------------------------------------------------------------------
    FILE*  file;
    string folder = "./prof";

    if (rank == 0) {
        struct stat st = {0};
        if (stat(folder.c_str(), &st) == -1) {
            mkdir(folder.c_str(), 0770);
        }

        string filename = folder + "/" + name_ + "_parent.csv";
        file            = fopen(filename.c_str(), "w+");

        if (file != nullptr) {
            time_map_["root"]->DumpParentality(file, 0);
            fclose(file);
        } else {
            printf("unable to open file %s !", filename.c_str());
        }
    }

    //-------------------------------------------------------------------------
    /** - do the IO of the timing */
    //-------------------------------------------------------------------------

    if (rank == 0) {
        string filename = "./prof/" + name_ + "_time.csv";
        file            = fopen(filename.c_str(), "w+");
    }
    // display the header
    if (rank == 0) {
        printf("===================================================================================================================================================\n");
        printf("        PROFILER %s \n", name_.c_str());
        // printf("\t-NAME-   \t\t\t-%% global-\t-%% local-\t-Total time-\t-Self time-\t-time/call-\t-Min tot time-\t-Max tot time-\t-Mean cnt-\n");
        printf("%25s|  %-13s\t%-13s\t%-13s\t%-13s\t%-13s\t%-13s\t%-13s\t%-13s\t%-13s\n", "-NAME-    ", "-% global-", "-% local-", "-Total time-", "-Self time-", "-time/call-", "-Min time-", "-Max time-", "-Mean cnt-", "-(MB/s)-");
    }
    // get the global timing
    real_t total_time = this->Time(ref);

    // display root with the total time
    time_map_["root"]->Disp(file, 0, total_time);
    // display footer
    if (rank == 0) {
        printf("===================================================================================================================================================\n");
        printf("%% global - %% of the total time passed inside or in its children (based on the mean time among processors\n");
        // printf("%% self glob - %% of the total time passed inside (children not included, based on the mean time among processors\n");
        printf("%% local - %% of the dad's time passed inside or in its children (from the mean time among processors\n");
        printf("Total time - the total time spend in that timer (averaged among the processors)\n");
        printf("Self time - the self time spend in that timer = children not included (averaged among the processors)\n");
        printf("Time/call - the total time spend in that timer per call of the timer (averaged among the processors)\n");
        printf("Min time - the min time / call spend in that timer among the processors\n");
        printf("Max time - the max time / call spend in that timer among the processors\n");
        printf("Mean cnt - the total number of time the timer has been called (averaged among the processors)\n");
        printf("===================================================================================================================================================\n");

        if (file != nullptr) {
            fclose(file);
        } else {
            printf("unable to open file for profiling !");
        }
    }
}
