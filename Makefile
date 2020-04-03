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
BUILDDIR := ./build
SRC_DIR := ./src
OBJ_DIR := ./build

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

#-----------------------------------------------------------------------------
## add the wanted folders - common folders
SRC := $(notdir $(wildcard $(SRC_DIR)/*.cpp))
HEAD := $(wildcard $(SRC_DIR)/*.hpp)
API := $(wildcard $(SRC_DIR)/*.h)

## generate object list
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)
OBJ := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
IN := $(SRC:%.cpp=$(OBJ_DIR)/%.in)

################################################################################
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp $(HEAD) $(API)
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -std=c++11 -fPIC -MMD -c $< -o $@

$(OBJ_DIR)/%.in : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INC) $(DEF) -std=c++11 -fPIC -MMD -E $< -o $@

################################################################################
default: $(TARGET)


# compile static and dynamic lib
all: $(TARGET)

preproc: $(IN)

$(TARGET): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB)

test:
	@echo $(SRC)

clean:
	@rm -f $(OBJ_DIR)/*.o
	@rm -f $(TARGET)

destroy:
	@rm -rf $(OBJ_DIR)/*.o
	@rm -rf $(OBJ_DIR)/*.d
	@rm -rf $(TARGET)
	@rm -rf $(OBJ_DIR)/*

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
	$(info - OBJ A2A = $(OBJ_A2A))
	$(info - OBJ NB = $(OBJ_NB))
	$(info - DEP = $(DEP))
	$(info - LGF_DATA = $(LGF_DATA))
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
