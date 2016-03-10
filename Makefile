

all:  
	cd src/; make;
	cp src/dwk++ ./


.PHONY:clean
clean:
	cd src/;make clean; 
	rm dwk++





