.DEFAULT_GOAL := all
FLAGS = -O3  # -D FAST -D DEBUG -D THROW 
CC = g++

thermal_funcs.so:
	cd src && $(CC) $(FLAGS) -fPIC -shared thermal_funcs.cpp -lgsl -lgslcblas -I./ -o ../lib/thermal_funcs.so

all: thermal_funcs.so

clean:
	rm -f ./lib/thermal_funcs.so
