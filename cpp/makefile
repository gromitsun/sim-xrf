# makefile

XRL_PATH = "./xraylib"
INCLUDE_PATH = ./Include
LIB_PATH = ./Lib


ifeq ($(OS),Windows_NT)
    #Windows stuff
    # RM = del
	DLIB = dll
	EXE = exe
else
    #Linux stuff
    # RM = rm
	DLIB = so
	EXE = out
endif

g++ = g++

all: main.o
	$(MAKE) -C ./src all
	$(g++) main.o -L./Lib -lsim -L$(XRL_PATH)/Lib -lxrl -o main.$(EXE) -fopenmp #-static
	

main.o:
	$(g++) -c main.cpp -I$(INCLUDE_PATH) -I$(XRL_PATH)/Include -std=c++11 -fopenmp

clean:
	$(MAKE) -C ./src clean
	$(RM) -rf *.o main.exe ./Lib/*.a ./Lib/*.$(DLIB) all

