#ifndef SRC_PROF_HPP_
#define SRC_PROF_HPP_
// c headers
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
// cpp headers
#include <cmath>
#include <iostream>
#include <list>
#include <map>
#include <string>

#include "murphy.hpp"

using std::map;
using std::string;

class TimerBlock {
   protected:
    bool   is_root_ = true;
    lid_t  count_   = 0;
    size_t memsize_ = 0;

    real_t t0_       = 0.0;
    real_t t1_       = 0.0;
    real_t time_acc_ = 0.0;
    real_t time_max_ = 0.0;
    real_t time_min_ = 0.0;

    string name_ = "noname";

    TimerBlock*              parent_ = nullptr;
    map<string, TimerBlock*> children_;

   public:
    explicit TimerBlock(string name);

    void Start();
    void Stop();
    void Reset();
    void AddMem(size_t mem);
    void Disp(FILE* file, const int level, const real_t totalTime);

    int    count() const { return count_; }
    bool   is_root() const { return is_root_; }
    string name() const { return name_; }
    real_t time_acc() const;
    real_t time_min() const;
    real_t time_max() const;

    void AddChild(TimerBlock* child);
    void SetParent(TimerBlock* parent);
    void DumpParentality(FILE* file, const int level);
};

class Prof {
   protected:
    map<string, TimerBlock*> time_map_;

    const string name_;
    void         CreateSingle_(string name);

   public:
    Prof();
    explicit Prof(const string myname);
    ~Prof();

    void Create(string name);
    void Create(string child, string parent);

    void Start(string name);
    void Stop(string name);
    void AddMem(string name, size_t mem);

    real_t Time(const std::string ref);

    void Disp();
    void Disp(const std::string ref);
};

#define m_profStart(prof, name)                            \
    ({                                                     \
        Prof*  m_profStart_prof_ = (Prof*)(prof);          \
        std::string m_profStart_name_ = (std::string)(name);         \
        if ((m_profStart_prof_) != nullptr) {              \
            (m_profStart_prof_)->Start(m_profStart_name_); \
        }                                                  \
    })

#define m_profStop(prof, name)                          \
    ({                                                  \
        Prof*  m_profStop_prof_ = (Prof*)(prof);        \
        std::string m_profStop_name_ = (std::string)(name);       \
        if ((m_profStop_prof_) != nullptr) {            \
            (m_profStop_prof_)->Stop(m_profStop_name_); \
        }                                               \
    })

#endif  // SRC_PROF_HPP_
