################################################################################
# @copyright Copyright © UCLouvain 2019
# 
# FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
# 
# Copyright (C) <2019> <Universite catholique de Louvain (UCLouvain), Belgique>
# 
# List of the contributors to the development of FLUPS, Description and complete License: see LICENSE file.
# 
# This program (FLUPS) is free software: 
# you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program (see COPYING file).  If not, 
# see <http://www.gnu.org/licenses/>.
# 
################################################################################

################################################################################
# ARCH DEPENDENT VARIABLES
ARCH_FILE ?= make_arch/make.docker_gcc

include $(ARCH_FILE)

################################################################################
# FROM HERE, DO NOT TOUCH
#-----------------------------------------------------------------------------
PREFIX ?= ./
NAME := murphy
# library naming
TARGET := $(NAME)

#-----------------------------------------------------------------------------
SRC_DIR := src
TEST_DIR := test
OBJ_DIR := build

## add the headers to the vpaths
INC := -I$(SRC_DIR)

#-----------------------------------------------------------------------------
#---- FFTW
FFTW_INC ?= /usr/include
FFTW_LIB ?= /usr/lib
FFTW_LIBNAME ?= -lfftw3_omp -lfftw3
INC += -I$(FFTW_INC)
LIB += -L$(FFTW_LIB) $(FFTW_LIBNAME) -Wl,-rpath,$(FFTW_LIB)

#---- HDF5
HDF5_INC ?= /usr/include
HDF5_LIB ?= /usr/lib
HDF5_LIBNAME ?= -lhdf5
INC += -I$(HDF5_INC)
LIB += -L$(HDF5_LIB) $(HDF5_LIBNAME) -Wl,-rpath,$(HDF5_LIB)

#---- P4EST
P4EST_INC ?= /usr/include
P4EST_LIB ?= /usr/lib
P4EST_LIBNAME ?= -lsc -lp4est
INC += -I$(P4EST_INC)
LIB += -L$(P4EST_LIB) $(P4EST_LIBNAME) -Wl,-rpath,$(P4EST_LIB)

#---- GTEST
GTEST_INC ?= /usr/include
GTEST_LIB ?= /usr/lib
GTEST_LIBNAME ?= -lgtest

#-----------------------------------------------------------------------------
## add the wanted folders - common folders
SRC := $(notdir $(wildcard $(SRC_DIR)/*.cpp))
HEAD := $(wildcard $(SRC_DIR)/*.hpp)
TSRC := $(notdir $(wildcard $(TEST_DIR)/$(SRC_DIR)/*.cpp))
THEAD := $(wildcard $(TEST_DIR)/$(SRC_DIR)/*.hpp)

## generate object list
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)
OBJ := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
IN := $(SRC:%.cpp=$(OBJ_DIR)/%.in)
TOBJ := $(TSRC:%.cpp=$(TEST_DIR)/$(OBJ_DIR)/%.o)
TDEP := $(TSRC:%.cpp=$(TEST_DIR)/$(OBJ_DIR)/%.d)

################################################################################
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -std=c++11 -fPIC -MMD -c $< -o $@

$(TEST_DIR)/$(OBJ_DIR)/%.o : $(TEST_DIR)/$(SRC_DIR)/%.cpp $(THEAD)
	$(CXX) $(CXXFLAGS) -I$(TEST_DIR)/$(SRC_DIR) $(INC) -I$(GTEST_INC) $(DEF) -std=c++11 -fPIC -MMD -c $< -o $@

$(OBJ_DIR)/%.in : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -std=c++11 -fPIC -MMD -E $< -o $@

################################################################################
default: $(TARGET)

all: $(TARGET)

preproc: $(IN)

$(TARGET): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB)

test: $(TOBJ) $(filter-out $(OBJ_DIR)/main.o,$(OBJ))
	$(CXX) $(LDFLAGS) $^ -o $(TARGET)_$@ $(LIB) -L$(GTEST_LIB) $(GTEST_LIBNAME) -Wl,-rpath,$(GTEST_LIB)

.PHONY: clean
clean:
	@rm -rf $(OBJ_DIR)/*.o
	@rm -rf $(TARGET)
	@rm -rf $(TEST_DIR)/$(OBJ_DIR)/*.o
	@rm -rf $(TARGET)_test

.PHONY: destroy
destroy:
	@rm -rf $(OBJ_DIR)/*.o
	@rm -rf $(OBJ_DIR)/*.d
	@rm -rf $(TARGET)
	@rm -rf $(TEST_DIR)/$(OBJ_DIR)/*.o
	@rm -rf $(TEST_DIR)/$(OBJ_DIR)/*.d
	@rm -rf $(TARGET)_test

.PHONY: logo
info: logo
	$(info prefix = $(PREFIX)/lib )
	$(info compiler = $(shell $(CXX) --version))
	$(info compil. flags = $(CXXFLAGS) $(INC) $(DEF) -fPIC -MMD)
	$(info linker flags = -shared $(LDFLAGS))
	$(info using arch file = $(ARCH_FILE) )
	$(info LGF path = $(LGF_PATH) )
	$(info ------------)
	$(info FFTW:)
	$(info - include: -I$(FFTW_INC) )
	$(info - lib: -L$(FFTW_LIB) $(FFTW_LIBNAME) -Wl,-rpath,$(FFTW_LIB))
	$(info ------------)
	$(info HDF5:)
	$(info - include: -I$(HDF5_INC) )
	$(info - lib: -L$(HDF5_LIB) $(HDF5_LIBNAME) -Wl,-rpath,$(HDF5_LIB))
	$(info ------------)
	$(info LIST OF OBJECTS:)
	$(info - SRC = $(SRC))
	$(info - OBJ = $(OBJ))
	$(info - DEP = $(DEP))
	$(info - test SRC = $(TSRC))
	$(info - test OBJ = $(TOBJ))
	$(info - test DEP = $(TDEP))
	$(info ------------)

.NOTPARALLEL: logo

logo:
	$(info  .----------------. .----------------. .----------------. .----------------. .----------------. .----------------. )
	$(info | .--------------. | .--------------. | .--------------. | .--------------. | .--------------. | .--------------. |)
	$(info | | ____    ____ | | | _____  _____ | | |  _______     | | |   ______     | | |  ____  ____  | | |  ____  ____  | |)
	$(info | ||_   \  /   _|| | ||_   _||_   _|| | | |_   __ \    | | |  |_   __ \   | | | |_   ||   _| | | | |_  _||_  _| | |)
	$(info | |  |   \/   |  | | |  | |    | |  | | |   | |__| |   | | |    | |__| |  | | |   | |__| |   | | |   \ \  / /   | |)
	$(info | |  | |\  /| |  | | |  | '    ' |  | | |   |  __ /    | | |    |  ___/   | | |   |  __  |   | | |    \ \/ /    | |)
	$(info | | _| |_\/_| |_ | | |   \ `--' /   | | |  _| |  \ \_  | | |   _| |_      | | |  _| |  | |_  | | |    _|  |_    | |)
	$(info | ||_____||_____|| | |    `.__.'    | | | |____| |___| | | |  |_____|     | | | |____||____| | | |   |______|   | |)
	$(info | |              | | |              | | |              | | |              | | |              | | |              | |)
	$(info | '--------------' | '--------------' | '--------------' | '--------------' | '--------------' | '--------------' |)
	$(info  '----------------' '----------------' '----------------' '----------------' '----------------' '----------------' )

-include $(DEP)