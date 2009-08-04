all: GF2D_test GF2DRad_test GF1D_test GF1DRad_test
	make -f Makefile_GF_test GF2D_test GF2DRad_test GF1D_test GF1DRad_test

GF2D_test : GF2D_test.o FirstPassageGreensFunction2D.o
	libtool --mode=link gcc -o GF2D_test GF2D_test.o findRoot.lo FirstPassageGreensFunction2D.o \
	funcSum.lo freeFunctions.lo -lgsl -lgslcblas -lstdc++

GF2D_test.o : GF2D_test.cpp
	gcc -c GF2D_test.cpp



GF2DRad_test : GF2DRad_test.o FirstPassagePairGreensFunction2D3.o
	libtool --mode=link gcc -o GF2DRad_test GF2DRad_test.o FirstPassagePairGreensFunction2D3.o \
        funcSum.lo freeFunctions.lo findRoot.lo -lgsl -lgslcblas -lstdc++

GF2DRad_test.o : GF2DRad_test.cpp
	gcc -c GF2DRad_test.cpp



GF1DRad_test : GF1DRad_test.o FirstPassageGreensFunction1DRad.o
	libtool --mode=link gcc -o GF1DRad_test GF1DRad_test.o FirstPassageGreensFunction1DRad.o \
	findRoot.lo -lgsl -lgslcblas -lstdc++

GF1DRad_test.o : GF1DRad_test.cpp
	gcc -c GF1DRad_test.cpp



GF1D_test : GF1D_test.o FirstPassageGreensFunction1D.o
	libtool --mode=link gcc -o GF1D_test GF1D_test.o FirstPassageGreensFunction1D.o \
	findRoot.lo -lgsl -lgslcblas -lstdc++

GF1D_test.o : GF1D_test.cpp
	gcc -c GF1D_test.cpp



FirstPassageGreensFunction2D.o :	FirstPassageGreensFunction2D.cpp \
					FirstPassageGreensFunction2D.hpp
	gcc -c FirstPassageGreensFunction2D.cpp

FirstPassageGreensFunction1D.o :	FirstPassageGreensFunction1D.cpp \
					FirstPassageGreensFunction1D.hpp
	gcc -c FirstPassageGreensFunction1D.cpp

FirstPassageGreensFunction1DRad.o :	FirstPassageGreensFunction1DRad.cpp \
					FirstPassageGreensFunction1DRad.hpp
	gcc -c FirstPassageGreensFunction1DRad.cpp

FirstPassagePairGreensFunction2D3.o :	FirstPassagePairGreensFunction2D3.cpp \
					FirstPassagePairGreensFunction2D3.hpp
	gcc -c FirstPassagePairGreensFunction2D3.cpp
