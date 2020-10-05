
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
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) -std=c++17 -fPIC -MMD -c $< -o $@

$(TEST_DIR)/$(OBJ_DIR)/%.o : $(TEST_DIR)/$(SRC_DIR)/%.cpp $(THEAD)
	$(CXX) $(CXXFLAGS) $(OPTS) -I$(TEST_DIR)/$(SRC_DIR) $(INC) -I$(GTEST_INC) $(DEF) -std=c++17 -fPIC -MMD -c $< -o $@

$(OBJ_DIR)/%.in : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(OPTS) $(INC) $(DEF) -std=c++17 -fPIC -MMD -E $< -o $@

################################################################################
default: $(TARGET)

all: $(TARGET)

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

# mkdir the needed dirs
$(shell   mkdir -p $(OBJ_DIR) $(TEST_DIR)/$(OBJ_DIR))