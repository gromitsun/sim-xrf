# makefile

XRL_PATH = "../xraylib"
SRC_DIR = "."
OBJ_DIR = "../obj"
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


all: $(OBJ_FILES) cwrapper.o
	ar crf ../Lib/libsim.a $(OBJ_FILES) cwrapper.o
	g++ $(OBJ_FILES) cwrapper.o -L$(XRL_PATH)/Lib -lxrl -shared -fPIC -o ../Lib/libsim.$(DLIB)
	
%.o: %.cpp
	@echo "============="
	@echo "Compiling $<"
	@g++ -c $< -fPIC -o $@ -I$(XRL_PATH)/Include -I../Include -std=c++11
	
cwrapper.o: cwrapper.cc
	@echo "============="
	@echo "Compiling $<"
	@g++ -c $< -fPIC -o $@ -I$(XRL_PATH)/Include -I../Include -std=c++11

clean:
	$(RM) -rf *.o *.$(DLIB) *.a all

