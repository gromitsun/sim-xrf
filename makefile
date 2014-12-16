# makefile

include variables.mk

all: main.o
	$(MAKE) -C ./src all
	$(g++) main.o -L./Lib -lsim -L$(XRL_LIB) -lxrl -o main.$(EXE) -fopenmp #-static
	echo export LD_LIBRARY_PATH="$(LIB_PATH)" > init.sh
	echo export OMP_NUM_THREADS=16 >> init.sh
	echo export PYTHONPATH="$(mkfile_path):\$$PYTHONPATH" >> init.sh

main.o: main.cpp
	$(g++) -c main.cpp -I$(INCLUDE_PATH) -I$(XRL_INCLUDE) -std=c++11 -fopenmp

clean:
	$(MAKE) -C ./src clean
	$(RM) -rvf *.o main.exe $(LIB_PATH)/*.a $(LIB_PATH)/*.$(DLIB) $(LIB_PATH) all
