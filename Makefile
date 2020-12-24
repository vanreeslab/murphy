# based on 
# file names: https://www.gnu.org/software/make/manual/html_node/File-Name-Functions.html
#

################################################################################
# ARCH DEPENDENT VARIABLES
ARCH_FILE ?= make_arch/make.docker_gcc

include $(ARCH_FILE)

################################################################################
# FROM HERE, DO NOT CHANGE
#-----------------------------------------------------------------------------
PREFIX ?= ./
NAME := murphy
# library naming
TARGET := $(NAME)

#-----------------------------------------------------------------------------
# get a list of all the source directories + the main one
SRC_DIR := src $(shell find src/** -type d)
# SUB_DIR := $(shell find $(SRC_DIR)/** -type d| sed 's/$(SRC_DIR)\///1')
TEST_DIR := test
OBJ_DIR := build



#-----------------------------------------------------------------------------
# the sources/headers are listed without the folder, vpath will find them
SRC := $(foreach dir,$(SRC_DIR),$(notdir $(wildcard $(dir)/*.cpp)))
HEAD := $(foreach dir,$(SRC_DIR),$(notdir $(wildcard $(dir)/*.hpp)))

# find the test sources, mandatory all in the same folder!
TSRC := $(foreach dir,$(SRC_DIR),$(notdir $(wildcard $(TEST_DIR)/$(dir)/*.cpp)))
THEAD := $(foreach dir,$(SRC_DIR),$(notdir $(wildcard $(TEST_DIR)/$(dir)/*.hpp)))

## generate object list
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)
OBJ := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
IN := $(SRC:%.cpp=$(OBJ_DIR)/%.in)
TIDY := $(SRC:%.cpp=$(OBJ_DIR)/%.tidy)
TOBJ := $(TSRC:%.cpp=$(TEST_DIR)/$(OBJ_DIR)/%.o)
TDEP := $(TSRC:%.cpp=$(TEST_DIR)/$(OBJ_DIR)/%.d)


#-----------------------------------------------------------------------------
# add the folders to the includes and to the vpath

## add the source dirs to the includes flags
INC := $(foreach dir,$(SRC_DIR),-I$(dir))
TINC := $(foreach dir,$(TEST_DIR)/$(SRC_DIR),-I$(dir))

## add them to the VPATH as well (https://www.gnu.org/software/make/manual/html_node/Selective-Search.html)
vpath %.hpp $(SRC_DIR) $(foreach dir,$(SRC_DIR),$(TEST_DIR)/$(dir))
vpath %.cpp $(SRC_DIR) $(foreach dir,$(SRC_DIR),$(TEST_DIR)/$(dir))

#-----------------------------------------------------------------------------
# LIBRARIES
#---- FLUPS
FLUPS_INC ?= /flups/include
FLUPS_LIB ?= /flups/lib
FLUPS_LIBNAME ?= -lflups_a2a
INC += -I$(FLUPS_INC)
LIB += -L$(FLUPS_LIB) $(FLUPS_LIBNAME) -Wl,-rpath,$(FLUPS_LIB)

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

#---- MPI - get flags
MPI_INC := $(shell mpic++ --showme:compile)

################################################################################
$(OBJ_DIR)/%.o : %.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) -std=c++17 -fPIC -MMD -c $< -o $@

$(TEST_DIR)/$(OBJ_DIR)/%.o : %.cpp $(HEAD) $(THEAD)
	$(CXX) $(CXXFLAGS) $(OPTS) $(TINC) $(INC) -I$(GTEST_INC) $(DEF) -std=c++17 -fPIC -MMD -c $< -o $@

$(OBJ_DIR)/%.in : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) -std=c++17 -fPIC -MMD -E $< -o $@

$(OBJ_DIR)/%.tidy : %.cpp $(HEAD)
	clang-tidy $< --format-style=.clang-format --checks=mpi-*,openmp-*,google-*,performance-* -- $(MPI_INC) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) $(MPI_INC) -std=c++17 -fPIC -MMD

################################################################################
default: $(TARGET)

all: $(TARGET)

.PHONY: tidy
tidy: $(TIDY)

.PHONY: preproc 
preproc: $(IN)

$(TARGET): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIB)

.PHONY: test 
test: $(TOBJ) $(filter-out $(OBJ_DIR)/main.o,$(OBJ))
	$(CXX) $(LDFLAGS) $^ -o $(TARGET)_$@ $(LIB) -L$(GTEST_LIB) $(GTEST_LIBNAME) -Wl,-rpath,$(GTEST_LIB)


.PHONY: clean 
clean:
	@rm -rf $(OBJ_DIR)/*.o
	@rm -rf $(OBJ_DIR)/*.tidy
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

.PHONY: logo info
info: logo
	$(info prefix = $(PREFIX)/lib )
	$(info compiler = $(shell $(CXX) --version))
	$(info compil. options = $(OPTS))
	$(info compil. flags = $(CXXFLAGS) $(OPTS) $(INC) $(DEF) -fPIC -MMD)
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
	$(info - SRC  = $(SRC))
	$(info - HEAD = $(HEAD))
	$(info - OBJ  = $(OBJ))
	$(info - DEP  = $(DEP))
	$(info - TEST_DIR = $(TEST_DIR))
	$(info - TEST_DIR = $(TEST_DIR)/$(OBJ_DIR))
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

# mkdir the needed dirs
$(shell   mkdir -p $(OBJ_DIR) $(TEST_DIR)/$(OBJ_DIR))
