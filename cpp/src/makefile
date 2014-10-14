# makefile

XRL_PATH = "../xraylib"
SRC_DIR = "."
OBJ_DIR = "../obj"
LIB_DIR = "../Lib/"
SRC_FILES = $(wildcard *.cpp)
OBJ_FILES = $(SRC_FILES:.cpp=.o)

ifeq ($(OS),Windows_NT)
	#Windows stuff
	# RM = del
	DLIB = dll
else
	#Linux stuff
	# RM = rm
	DLIB = so
endif

g++ = g++

all: $(OBJ_FILES) cwrapper.o ../Lib
	ar crf ../Lib/libsim.a $(OBJ_FILES)
	# $(g++) $(OBJ_FILES) cwrapper.o -L$(XRL_PATH)/Lib -lxrl -shared -fPIC -o ../Lib/libsim.$(DLIB)
	$(g++) cwrapper.o ../Lib/libsim.a -L$(XRL_PATH)/Lib -lxrl -shared -fPIC -o ../Lib/libsim.$(DLIB) -fopenmp
	
%.o: %.cpp
	@echo "============="
	@echo "Compiling $<"
	@$(g++) -c $< -fPIC -o $@ -I$(XRL_PATH)/Include -I../Include -std=c++11 -fopenmp
	
cwrapper.o: cwrapper.cc
	@echo "============="
	@echo "Compiling $<"
	@$(g++) -c $< -fPIC -o $@ -I$(XRL_PATH)/Include -I../Include -std=c++11 -fopenmp

../Lib:
	- mkdir ../Lib

clean:
	$(RM) -rf *.o *.$(DLIB) *.a all

