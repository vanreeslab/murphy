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

class TimerBlock {
   protected:
    // bool        is_root_  = true;      //!< indicate if the block is root
    lid_t       count_    = 0;         //!< the number of times this block has been called
    size_t      memsize_  = 0;         //!< the memory size associated with a memory operation
    real_t      t0_       = -1.0;      //!< temp start time of the block
    real_t      t1_       = -1.0;      //!< temp stop time of the block
    real_t      time_acc_ = 0.0;       //!< accumulator to add the time accumulation
    std::string name_     = "noname";  //!< the default name of the block

    TimerBlock* parent_ = nullptr;  //!< the link to the parent blocks

    std::map<std::string, TimerBlock*> children_;  //!< the link to the children blocks

   public:
    explicit TimerBlock(std::string name);
    ~TimerBlock();

    void Start();
    void Stop();
    // void Reset();
    // void AddMem(size_t mem);
    void Disp(FILE* file, const int level, const real_t totalTime);

    int         count() const { return count_; }
    // bool        is_root() const { return is_root_; }
    std::string name() const { return name_; }
    real_t      time_acc() const;
    real_t      time_min() const;
    real_t      time_max() const;
    TimerBlock* parent() const { return parent_; }

    TimerBlock* AddChild(std::string child_name);

    void SetParent(TimerBlock* parent);
    // void DumpParentality(FILE* file, const int level) const;
};

class Prof {
   protected:
    std::map<std::string, TimerBlock*> time_map_;

    TimerBlock* current_;  //!< this is a pointer to the last TimerBlock

    const std::string name_;

    // void              CreateSingle_(string name);

   public:
    explicit Prof();
    explicit Prof(const std::string myname);
    ~Prof();

    // void Create(string name);
    // void Create(string child, string parent);

    void Start(std::string name);
    void Stop(std::string name);
    // void AddMem(string name, size_t mem);

    // real_t Time(const std::string ref);

    void Disp() const;
    // void Disp(const std::string ref) const;
};

#define m_profStart(prof, name)                              \
    ({                                                       \
        Prof*       m_profStart_prof_ = (Prof*)(prof);       \
        std::string m_profStart_name_ = (std::string)(name); \
        if ((m_profStart_prof_) != nullptr) {                \
            (m_profStart_prof_)->Start(m_profStart_name_);   \
        }                                                    \
    })

#define m_profStop(prof, name)                              \
    ({                                                      \
        Prof*       m_profStop_prof_ = (Prof*)(prof);       \
        std::string m_profStop_name_ = (std::string)(name); \
        if ((m_profStop_prof_) != nullptr) {                \
            (m_profStop_prof_)->Stop(m_profStop_name_);     \
        }                                                   \
    })

// #define m_profCreate(prof, name)                              \
//     ({                                                        \
//         Prof*       m_profCreate_prof_ = (Prof*)(prof);       \
//         std::string m_profCreate_name_ = (std::string)(name); \
//         if ((m_profCreate_prof_) != nullptr) {                \
//             (m_profCreate_prof_)->Create(m_profCreate_name_); \
//         }                                                     \
//     })

// #define m_profCreateParent(prof, name_parent, name_child)                                                        \
//     ({                                                                                                           \
//         Prof*       m_profCreateParent_prof_        = (Prof*)(prof);                                             \
//         std::string m_profCreateParent_name_parent_ = (std::string)(name_parent);                                \
//         std::string m_profCreateParent_name_child_  = (std::string)(name_child);                                 \
//         if ((m_profCreateParent_prof_) != nullptr) {                                                             \
//             (m_profCreateParent_prof_)->Create(m_profCreateParent_name_child_,m_profCreateParent_name_parent_); \
//         }                                                                                                        \
//     })

#endif  // SRC_PROF_HPP_
