.DEFAULT_GOAL := all
FLAGS = -O3  # -D DEBUG -D THROW 
CC = g++
PYTHON = `pkg-config --cflags python`

thermal_funcs.o: thermal_funcs.h
	$(CC) $(FLAGS) -c thermal_funcs.cpp -I./
	
thermal_funcs_wrap.cxx:
	swig -c++ -python -py3 thermal_funcs.i
	
thermal_funcs_wrap.o: thermal_funcs_wrap.cxx
	$(CC) -c -fPIC thermal_funcs.cpp thermal_funcs_wrap.cxx -lgsl -lgslcblas $(PYTHON) -I./
	
_thermal_funcs.so: thermal_funcs_wrap.o thermal_funcs.o
	$(CC) -shared thermal_funcs.o thermal_funcs_wrap.o -lgsl -lgslcblas -o _thermal_funcs.so
	cd ../thermal_funcs && ln -sf ../src/_thermal_funcs.so .

thermal_funcs.so:
	$(CC) $(FLAGS) -fPIC -shared thermal_funcs.cpp -lgsl -lgslcblas -I./ -o ../lib/thermal_funcs.so

all: thermal_funcs.so _thermal_funcs.so

clean:
	rm -f ../lib/thermal_funcs.so
	rm -f _thermal_funcs.so
	rm -f thermal_funcs_wrap.o
	rm -f thermal_funcs.o
	rm -f thermal_funcs.py
	rm -f thermal_funcs_wrap.cxx
	rm -f thermal_funcs.pyc
	
