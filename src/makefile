.DEFAULT_GOAL := all

FLAGS = -O3  # -D DEBUG -D THROW 
CC = g++
PYTHON = `pkg-config --cflags python`
MATH = /usr/local/Wolfram/Mathematica/11.1/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/CompilerAdditions/wscc
MATH_INC = -I./ -I/usr/local/Wolfram/Mathematica/11.1/SystemFiles/Links/WSTP/DeveloperKit/Linux/CompilerAdditions/ -I/usr/local/Wolfram/Mathematica/11.1/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/CompilerAdditions/

zeta.o: zeta.h
	$(CC) $(FLAGS) -fPIC -c zeta.cpp -lgsl -lgslcblas -I./  

thermal_funcs.o: thermal_funcs.h
	$(CC) $(FLAGS) -fPIC -c thermal_funcs.cpp -I./
	
thermal_funcs_wrap.cxx:
	swig -c++ -python -py3 thermal_funcs.i
	
thermal_funcs_wrap.o: thermal_funcs_wrap.cxx
	$(CC) -c -fPIC thermal_funcs.cpp thermal_funcs_wrap.cxx -lgsl -lgslcblas $(PYTHON) -I./
	
_thermal_funcs.so: thermal_funcs_wrap.o thermal_funcs.o zeta.o
	$(CC) -shared thermal_funcs.o thermal_funcs_wrap.o zeta.o -lgsl -lgslcblas -I./ -o _thermal_funcs.so
	cd ../thermal_funcs && ln -sf ../src/_thermal_funcs.so .

thermal_funcs.so: thermal_funcs.o zeta.o
	$(CC) $(FLAGS) -fPIC -shared thermal_funcs.o zeta.o -lgsl -lgslcblas -I./ -o ../lib/thermal_funcs.so

math.exe: thermal_funcs.o zeta.o
	$(MATH) -b64 -o math.exe math.tm math.c thermal_funcs.o zeta.o -lgsl -lgslcblas $(MATH_INC)

all: thermal_funcs.so _thermal_funcs.so

clean:
	rm -f ../lib/thermal_funcs.so
	rm -f _thermal_funcs.so
	rm -f thermal_funcs_wrap.o
	rm -f thermal_funcs.o
	rm -f zeta.o
	rm -f thermal_funcs.py
	rm -f thermal_funcs_wrap.cxx
	rm -f thermal_funcs.pyc
	rm -f math.exe
