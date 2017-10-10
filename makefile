.DEFAULT_GOAL := all

all:
	$(MAKE) -C src

clean:
	$(MAKE) -C src clean
	
math.exe:
	$(MAKE) -C src math.exe
	
thermal_funcs.so:
	$(MAKE) -C src thermal_funcs.so
	
_thermal_funcs.so:
	$(MAKE) -C src _thermal_funcs.so
	
example:
	$(MAKE) -C src example
