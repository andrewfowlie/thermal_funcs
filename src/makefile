.DEFAULT_GOAL := lib
.PHONY: all clean mathematica python lib example

FLAGS = -O3  # -D DEBUG -D THROW 
CC = g++
PYTHON = $(shell pkg-config --cflags python)
MATH_INC = $(shell echo 'FileNameJoin[{$$InstallationDirectory, "SystemFiles", "Links", "WSTP", "DeveloperKit", $$SystemID, "CompilerAdditions"}]' | math -batchoutput | tr -d '"')
MATH = $(MATH_INC)/wscc

zeta.o: zeta.h
	$(CC) $(FLAGS) -fPIC -c zeta.cpp -lgsl -lgslcblas -I./  

thermal_funcs.o: thermal_funcs.h
	$(CC) $(FLAGS) -fPIC -c thermal_funcs.cpp -I./
	
thermal_funcs_wrap.cxx:
	swig -c++ -python -py3 thermal_funcs.i
	
thermal_funcs_wrap.o: thermal_funcs_wrap.cxx
	$(CC) -c -fPIC derivatives.cpp thermal_funcs.cpp thermal_funcs_wrap.cxx -lgsl -lgslcblas $(PYTHON) -I./
	
_thermal_funcs.so: thermal_funcs_wrap.o thermal_funcs.o zeta.o derivatives.o 
	$(CC) -shared derivatives.o thermal_funcs.o thermal_funcs_wrap.o zeta.o -lgsl -lgslcblas -I./ -o _thermal_funcs.so
	cd ../thermal_funcs && ln -sf ../src/_thermal_funcs.so .

thermal_funcs.so: thermal_funcs.o zeta.o derivatives.o
	$(CC) $(FLAGS) -fPIC -shared derivatives.o thermal_funcs.o zeta.o -lgsl -lgslcblas -I./ -o ../lib/thermal_funcs.so

math.exe: thermal_funcs.o zeta.o
	$(MATH) -b64 -o math.exe math.tm math.c -Wl,thermal_funcs.o,zeta.o -lgsl -lgslcblas -I$(MATH_INC) -I./
	
example: thermal_funcs.so
	# $(CC) $(FLAGS) example.c thermal_funcs.o zeta.o -lgsl -lgslcblas -I./ -o ../bin/example 
	$(CC) $(FLAGS) example.c -L$(shell echo $$PWD)/../lib -I./ -o ../bin/example -l:thermal_funcs.so -lgsl -lgslcblas

derivatives.o:
	$(CC) $(FLAGS) -fPIC -c derivatives.cpp -I./ 

mathematica: math.exe
python: _thermal_funcs.so
lib: thermal_funcs.so
all: lib python mathematica example

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
	rm -f ../bin/example
