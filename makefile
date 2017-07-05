.DEFAULT_GOAL := all

all: 
	$(MAKE) -C src

clean:
	$(MAKE) -C src clean
