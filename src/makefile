# makefile

include ../variables.mk

SRC_DIR = "."
OBJ_DIR = "../obj"
SRC_FILES = $(wildcard *.cpp)
OBJ_FILES = $(SRC_FILES:.cpp=.o)


all: $(OBJ_FILES) cwrapper.o $(LIB_PATH)
	ar crf $(LIB_PATH)/libsim.a $(OBJ_FILES)
	$(g++) $(OBJ_FILES) cwrapper.o -L$(XRL_LIB) -lxrl -shared -fPIC -o $(LIB_PATH)/libsim.$(DLIB) -fopenmp
	# $(g++) cwrapper.o $(LIB_PATH)/libsim.a -L$(XRL_LIB) -lxrl -shared -fPIC -o $(LIB_PATH)/libsim.$(DLIB) -fopenmp # works for mingw-w64 on windows, not for linux.

%.o: %.cpp
	@echo "============="
	@echo "Compiling $<"
	@$(g++) -c $< -fPIC -o $@ -I$(XRL_INCLUDE) -I$(INCLUDE_PATH) -std=c++11 -fopenmp

cwrapper.o: cwrapper.cc
	@echo "============="
	@echo "Compiling $<"
	@$(g++) -c $< -fPIC -o $@ -I$(XRL_INCLUDE) -I$(INCLUDE_PATH) -std=c++11 -fopenmp

$(LIB_PATH):
	- mkdir $(LIB_PATH)

clean:
	$(RM) -rvf *.o $(LIB_PATH)/*.$(DLIB) $(LIB_PATH)/*.a all