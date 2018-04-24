.DEFAULT_GOAL := lib
.PHONY: all clean mathematica python lib example

CXXFLAGS = -O3 -std=c++11 # -D DEBUG -D THROW
CFLAGS = -O3 # -D DEBUG -D THROW
CXX = g++
PYTHON = $(shell pkg-config --cflags python)
MATH_INC = $(shell echo 'FileNameJoin[{$$InstallationDirectory, "SystemFiles", "Links", "WSTP", "DeveloperKit", $$SystemID, "CompilerAdditions"}]' | math -batchoutput | tr -d '"')
MATH = $(MATH_INC)/wscc
OBJECTS = zeta.o thermal_funcs.o derivatives.o

$(OBJECTS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -fPIC -I./ -c $< 
	
thermal_funcs_wrap.cxx: thermal_funcs.i
	swig -c++ -python -py3 $<
	
thermal_funcs_wrap.o: thermal_funcs_wrap.cxx derivatives.cpp thermal_funcs.cpp
	$(CXX) $(CXXFLAGS) -fPIC $(PYTHON) -I./ -c $^
	
../thermal_funcs/_thermal_funcs.so: thermal_funcs_wrap.o thermal_funcs.o zeta.o derivatives.o 
	$(CXX) $(CXXFLAGS) -shared $^ -lgsl -lgslcblas -I./ -o $@

../lib/thermal_funcs.so: thermal_funcs.o zeta.o derivatives.o
	$(CXX) $(CXXFLAGS) -fPIC -shared $^ -lgsl -lgslcblas -I./ -o $@

math.exe: thermal_funcs.o zeta.o derivatives.o math.tm math.c
	$(MATH) $(CFLAGS) -b64 -Wl,thermal_funcs.o,zeta.o,derivatives.o -lgsl -lgslcblas -I$(MATH_INC) -I./ math.tm math.c -o $@
	
../bin/example: example.cpp ../lib/thermal_funcs.so
	# $(CXX) $(CXXFLAGS) $< -lgsl -lgslcblas -I./ thermal_funcs.o zeta.o derivatives.o -o $@
	$(CXX) $(CXXFLAGS) $< -L$(shell echo $$PWD)/../lib -l:thermal_funcs.so -lgsl -lgslcblas -I./ -o $@

example: ../bin/example
mathematica: math.exe
python: ../thermal_funcs/_thermal_funcs.so
lib: ../lib/thermal_funcs.so
all: lib python mathematica example

clean:
	-rm ../lib/thermal_funcs.so
	-rm ../thermal_funcs/_thermal_funcs.so
	-rm thermal_funcs_wrap.o
	-rm thermal_funcs.o
	-rm zeta.o
	-rm derivatives.o
	-rm thermal_funcs.py
	-rm thermal_funcs_wrap.cxx
	-rm math.exe
	-rm ../bin/example
