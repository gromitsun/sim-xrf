sim-xrf
=======
This is a python/C++ program for simulating X-Ray Fluorescence (XRF) spectrum acquired at a synchrotron radiation beamline.


Installation
=======
The source files need to be compiled with c++11 (tested with g++ version 4.7 and mingw-w64 x86_64-4.9.1-posix-seh-rt_v3-rev1). 

Xraylib with python bindings is also needed at compiling.

Xraylib is available at https://github.com/tschoonj/xraylib.

Compiling Xraylib
-------
Pre-compiled Xraylib libraries for Ubuntu Linux and Windows are in the /xraylib/ folder. If the libraries do not work, Xraylib has to be recompiled manually.

After compilation of Xraylib, change the XRL_PATH in makefile and src/makefile to the actual path of the installed xraylib.
Or copy the library files to xraylib/Lib/linux or xraylib/Lib/windows.

Compiling sim-xrf
-------
In the top directory do "make". If your c++ compiler is not called with "g++", you must edit the "g++ = g++" line in makefile and src/makefile.

Initialization
-------
Initialize the program by running "source init.sh" before running the program.


Running sim-xrf
=======
Sim-xrf can be run in C++ mode and Python mode. In both modes, an input file defining the parameters needed in the calculation is required.
The output of the C++ mode is a text file; the Python mode can in addition plot the spectrum or read the output files and plot them. 

Input file
-------
The input file is a text file defining the parameters needed in the calculation. See the example "input.txt" file for instructions.

Running in C++ mode
-------
On linux:
    ./main.out [input file name] [output file name]

On windows:
    main.exe [input file name] [output file name]

Running in Python mode
-------
# Calculation + plotting:
    python simpy.py [input file name] [output file name]
    
# Reading output file + plotting:
    python read.py [input file name] [output file name]
    
# Import as a module:
## Calculate spectrum
    from python.pyapi import calc
    
    nout = 3000, 30, 500, 500  # N of channels, N of Z, N of lines, N of thetas
    xlim = 0, 11
    ylim = 1e-14, 1e-4
    
    spec = calc(input_file="input.txt", output_file="output.txt")
    spec.show(xlim=xlim, ylim=ylim)

## Read/plot output files
    from read import read, plot
        
    fname = "output.txt"
    xlim = 0, 11
    ylim = 1e-14, 1e-4
    
    # Read data file
    ev, a, labels = read(fname)
    
    # Plot spectrum
    plt.figure()
    plt.title('Spectrum read from %s' % fname)
    plot(ev, a, labels, xlim=xlim, ylim=ylim)