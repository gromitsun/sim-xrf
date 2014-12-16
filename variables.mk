# # # variables used in makefile

# Makefile Paths
mkfile_path := $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))
INCLUDE_PATH = $(abspath $(mkfile_path)/Include)
LIB_PATH = $(abspath $(mkfile_path)/Lib)

# XrayLib paths
XRL_PATH = "D:\Program Files\xraylib 64-bit"



ifeq ($(OS),Windows_NT)
	#Windows stuff
	# RM = del
	DLIB = dll
	EXE = exe
	XRL_INCLUDE = $(XRL_PATH)/Include
	XRL_LIB = $(XRL_PATH)/Lib
else
	#Linux stuff
	# RM = rm
	DLIB = so
	EXE = out
	XRL_INCLUDE = $(XRL_PATH)/include
	XRL_LIB = $(XRL_PATH)/src/.lib
endif

# C++11 compiler
g++ = g++
