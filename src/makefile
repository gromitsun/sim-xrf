# makefile
mkfile_path := $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))
XRL_PATH = $(abspath $(mkfile_path)/../xraylib)
INCLUDE_PATH = $(abspath $(mkfile_path)/../Include)
LIB_DIR = $(abspath $(mkfile_path)/../Lib)
SRC_DIR = "."
OBJ_DIR = "../obj"
SRC_FILES = $(wildcard *.cpp)
OBJ_FILES = $(SRC_FILES:.cpp=.o)

ifeq ($(OS),Windows_NT)
	#Windows stuff
	# RM = del
	DLIB = dll
	XRL_LIB = windows
else
	#Linux stuff
	# RM = rm
	DLIB = so
	XRL_LIB = linux
endif

g++ = g++

all: $(OBJ_FILES) cwrapper.o $(LIB_DIR)
	ar crf $(LIB_DIR)/libsim.a $(OBJ_FILES)
	# $(g++) $(OBJ_FILES) cwrapper.o -L$(XRL_PATH)/Lib -lxrl -shared -fPIC -o $(LIB_DIR)/libsim.$(DLIB)
	$(g++) cwrapper.o $(LIB_DIR)/libsim.a -L$(XRL_PATH)/Lib/$(XRL_LIB) -lxrl -shared -fPIC -o $(LIB_DIR)/libsim.$(DLIB) -fopenmp

%.o: %.cpp
	@echo "============="
	@echo "Compiling $<"
	@$(g++) -c $< -fPIC -o $@ -I$(XRL_PATH)/Include -I$(INCLUDE_PATH) -std=c++11 -fopenmp

cwrapper.o: cwrapper.cc
	@echo "============="
	@echo "Compiling $<"
	@$(g++) -c $< -fPIC -o $@ -I$(XRL_PATH)/Include -I$(INCLUDE_PATH) -std=c++11 -fopenmp

$(LIB_DIR):
	- mkdir $(LIB_DIR)

clean:
	$(RM) -rvf *.o $(LIB_DIR)/*.$(DLIB) $(LIB_DIR)/*.a all