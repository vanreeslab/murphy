# Makefile MURPHY
#------------------------------------------------------------------------------
# useful links: 
# - automatic vars: https://www.gnu.org/software/make/manual/html_node/Automatic-Variables.html
# - file names: https://www.gnu.org/software/make/manual/html_node/File-Name-Functions.html
#------------------------------------------------------------------------------

################################################################################
# ARCH DEPENDENT VARIABLES
ARCH_FILE ?= make_arch/make.docker_gcc

#only include the ARCH_FILE when we do not clean or destroy
ifneq ($(MAKECMDGOALS),clean)
include $(ARCH_FILE)
endif

################################################################################
# FROM HERE, DO NOT CHANGE
#-------------------------------------------------------------------------------
PREFIX ?= ./
NAME := murphy
# library naming
TARGET := $(NAME)
# git commit
GIT_COMMIT ?= $(shell git describe --always --dirty)

#-------------------------------------------------------------------------------
# get a list of all the source directories + the main one
SRC_DIR := src $(shell find src/** -type d)
TEST_DIR := test
OBJ_DIR := build

#-------------------------------------------------------------------------------
# the sources/headers are listed without the folder, vpath will find them
SRC := $(foreach dir,$(SRC_DIR),$(notdir $(wildcard $(dir)/*.cpp)))
HEAD := $(foreach dir,$(SRC_DIR),$(notdir $(wildcard $(dir)/*.hpp)))

# find the test sources, mandatory all in the same folder!
TSRC := $(foreach dir,$(SRC_DIR),$(notdir $(wildcard $(TEST_DIR)/$(dir)/*.cpp)))
THEAD := $(foreach dir,$(SRC_DIR),$(notdir $(wildcard $(TEST_DIR)/$(dir)/*.hpp)))

## generate object list
OBJ := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEP := $(SRC:%.cpp=$(OBJ_DIR)/%.d)
# ASS := $(SRC:%.cpp=$(OBJ_DIR)/%.s)
CDB := $(SRC:%.cpp=$(OBJ_DIR)/%.o.json)
TIDY := $(SRC:%.cpp=$(OBJ_DIR)/%.tidy)

TOBJ := $(TSRC:%.cpp=$(TEST_DIR)/$(OBJ_DIR)/%.o)
TDEP := $(TSRC:%.cpp=$(TEST_DIR)/$(OBJ_DIR)/%.d)
TCDB := $(TSRC:%.cpp=$(TEST_DIR)/$(OBJ_DIR)/%.o.json)

#-------------------------------------------------------------------------------
# add the folders to the includes and to the vpath

## add the source dirs to the includes flags
INC := $(foreach dir,$(SRC_DIR),-I$(dir))
TINC := $(foreach dir,$(TEST_DIR)/$(SRC_DIR),-I$(dir))

## add them to the VPATH as well (https://www.gnu.org/software/make/manual/html_node/Selective-Search.html)
vpath %.hpp $(SRC_DIR) $(foreach dir,$(SRC_DIR),$(TEST_DIR)/$(dir))
vpath %.cpp $(SRC_DIR) $(foreach dir,$(SRC_DIR),$(TEST_DIR)/$(dir))

#-------------------------------------------------------------------------------
# LIBRARIES
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

################################################################################
# mandatory flags
M_FLAGS := -std=c++17 -fPIC -DGIT_COMMIT=\"$(GIT_COMMIT)\"

#-------------------------------------------------------------------------------
# compile + dependence + json file
$(OBJ_DIR)/%.o : %.cpp $(HEAD)
ifeq ($(shell $(CXX) -v 2>&1 | grep -c "clang"), 1)
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) $(M_FLAGS) -MMD -MF $(OBJ_DIR)/$*.d -MJ $(OBJ_DIR)/$*.o.json -c $< -o $@
else
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) $(M_FLAGS) -MMD -c $< -o $@
endif

# json only
$(OBJ_DIR)/%.o.json :  %.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) $(M_FLAGS) -MJ $@ -E $< -o $(OBJ_DIR)/$*.ii

#-------------------------------------------------------------------------------
# tests
# compile + dependence + json
$(TEST_DIR)/$(OBJ_DIR)/%.o : %.cpp $(HEAD) $(THEAD)
ifeq ($(shell $(CXX) -v 2>&1 | grep -c "clang"), 1)
	$(CXX) $(CXXFLAGS) $(OPTS) $(TINC) $(INC) -I$(GTEST_INC) $(DEF) $(M_FLAGS) -MMD -MF $(OBJ_DIR)/$*.d -MJ $(OBJ_DIR)/$*.o.json -c $< -o $@
else
	$(CXX) $(CXXFLAGS) $(OPTS) $(TINC) $(INC) -I$(GTEST_INC) $(DEF) $(M_FLAGS) -MMD -c $< -o $@
endif

# json only
$(TEST_DIR)/$(OBJ_DIR)/%.o.json : %.cpp $(HEAD) $(THEAD)
	$(CXX) $(CXXFLAGS) $(OPTS) $(TINC) $(INC) -I$(GTEST_INC) $(DEF) $(M_FLAGS) -MJ $@ -E $< -o $(TEST_DIR)/$(OBJ_DIR)/$*.ii

# clang-tidy files, define the MPI_INC which is only for this target
$(OBJ_DIR)/%.tidy : %.cpp $(HEAD)
	clang-tidy $< --format-style=.clang-format --checks=all*

################################################################################
.PHONY: default
default: $(TARGET)

.PHONY: all
all: $(TARGET) compdb

#-------------------------------------------------------------------------------
# the main target
$(TARGET): $(OBJ)
	$(CXX) $(LDFLAGS) $^ $(LIB) -o $@

#-------------------------------------------------------------------------------
# clang stuffs
.PHONY: tidy
tidy: $(TIDY)

# for the sed commande, see https://sarcasm.github.io/notes/dev/compilation-database.html#clang and https://sed.js.org
.PHONY: compdb
compdb: $(CDB)
	sed -e '1s/^/[\n/' -e '$$s/,$$/\n]/' $^ > compile_commands.json

# build the full compilation data-base, need to know the test libs!!
.PHONY: compdb_full
compdb_full: $(CDB) $(TCDB)
	sed -e '1s/^/[\n/' -e '$$s/,$$/\n]/' $^ > compile_commands.json

#-------------------------------------------------------------------------------
.PHONY: test 
test: $(TOBJ) $(filter-out $(OBJ_DIR)/main.o,$(OBJ))
	$(CXX) $(LDFLAGS) $^ -o $(TARGET)_$@ $(LIB) -L$(GTEST_LIB) $(GTEST_LIBNAME) -Wl,-rpath,$(GTEST_LIB)

#-------------------------------------------------------------------------------
#clean
.PHONY: clean 
clean:
	@rm -rf $(OBJ_DIR)/*
	@rm -rf $(OBJ_DIR)/*
	@rm -rf $(TEST_DIR)/$(OBJ_DIR)/*.o
	@rm -rf $(TARGET)
	@rm -rf $(TARGET)_test

#-------------------------------------------------------------------------------
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
