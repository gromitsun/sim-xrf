# makefile
mkfile_path := $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))
XRL_PATH = $(abspath $(mkfile_path)/xraylib)
INCLUDE_PATH = $(abspath $(mkfile_path)/Include)
LIB_PATH = $(abspath $(mkfile_path)/Lib)


ifeq ($(OS),Windows_NT)
	#Windows stuff
	# RM = del
	DLIB = dll
	EXE = exe
	XRL_LIB = windows
else
	#Linux stuff
	# RM = rm
	DLIB = so
	EXE = out
	XRL_LIB = linux
endif

g++ = g++

all: main.o
	$(MAKE) -C ./src all
	$(g++) main.o -L./Lib -lsim -L$(XRL_PATH)/Lib/$(XRL_LIB) -lxrl -o main.$(EXE) -fopenmp #-static
	echo export LD_LIBRARY_PATH=$(LIB_PATH) > init.sh
	echo export OMP_NUM_THREADS=16 >> init.sh

main.o: main.cpp
	$(g++) -c main.cpp -I$(INCLUDE_PATH) -I$(XRL_PATH)/Include -std=c++11 -fopenmp

clean:
	$(MAKE) -C ./src clean
	$(RM) -rvf *.o main.exe $(LIB_PATH)/*.a $(LIB_PATH)/*.$(DLIB) $(LIB_PATH) all
