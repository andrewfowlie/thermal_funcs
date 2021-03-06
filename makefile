.DEFAULT_GOAL := lib
.PHONY: all clean mathematica python lib example fortran docs

all:
	$(MAKE) -C src all

clean:
	$(MAKE) -C src clean
	
mathematica:
	$(MAKE) -C src mathematica
	
python:
	$(MAKE) -C src python
	
lib:
	$(MAKE) -C src lib
	
example:
	$(MAKE) -C src example

fortran:
	$(MAKE) -C src fortran

docs: doxygen.config
	doxygen $<
	xdg-open ./docs/html/index.html

test: python
	python thermal_funcs/thermal_funcs.py
