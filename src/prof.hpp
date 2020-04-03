/**
 * @file Profiler.hpp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @copyright Copyright © UCLouvain 2019
 * 
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 * 
 * Copyright (C) <2019> <Universite catholique de Louvain (UCLouvain), Belgique>
 * 
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE file.
 * 
 * This program (FLUPS) is free software: 
 * you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program (see COPYING file).  If not, 
 * see <http://www.gnu.org/licenses/>.
 * 
 */

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

class TimerAgent {
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

    TimerAgent*              parent_ = nullptr;
    map<string, TimerAgent*> children_;

   public:
    explicit TimerAgent(string name);

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

    void AddChild(TimerAgent* child);
    void SetParent(TimerAgent* parent);
    void DumpParentality(FILE* file, const int level);
};

class Prof {
   protected:
    map<string, TimerAgent*> time_map_;

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

#endif  // SRC_PROF_HPP_
