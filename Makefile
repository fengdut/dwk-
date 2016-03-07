PRGM = dwk++

CPPFLAGS =  -g -O0
LDFLAGS = 
CFLAGS = -fopenmp 


LIBCONFIG = -I${LIBCONFIG_DIR}/include  -L${LIBCONFIG_DIR}/lib -lconfig -lconfig++
COMPILE = $(CXX) $(CPPFLAGS) $(CFLAGS) $(LIBCONFIG) -c 



all :  
	cd src/; make;
	cp src/dwk++ ./





.PHONY:clean
clean:
	cd src/;make clean; 
	rm dwk++





