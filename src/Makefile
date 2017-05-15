PRGM = dwk++

CPPFLAGS =   
LDFLAGS = 
#CFLAGS = -fopenmp  -framework OpenCL -O2 -ftree-vectorize -fopt-info-vec 
#CFLAGS = -fopenmp   -O0 -ftree-vectorize
CFLAGS =     -fopenmp 
#CFLAGS = -fopenmp   -O0 


#CXX=/usr/local/bin/g++-6

LIBCONFIG = -I${LIBCONFIG_DIR}/include  -L${LIBCONFIG_DIR}/lib -lconfig -lconfig++
NETCDF = -I${NETCDF_DIR}/include -L${NETCDF_DIR}/lib -lnetcdf -L${HDF5_DIR}/lib

MKL_ROOT =/opt/intel/mkl

COMPILE = $(CXX) $(CPPFLAGS) $(CFLAGS) -I${MKL_ROOT}/include  $(LIBCONFIG) $(NETCDF) -c 

LINKCC = $(CXX) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS)  $(LIBCONFIG) $(NETCDF) -lcurl  -lhdf5 -lhdf5_hl

LIBA = 

SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS)) ellpk.o const.o 

all :  $(PRGM)


$(PRGM):$(OBJS)
	$(LINKCC) $(OBJS) $(LIBA) $(HDF5L)  -o $(PRGM)

%.o:%.cpp
	$(COMPILE) $< 

%.o:%.c
	$(COMPILE) $< 

.PHONY:clean
clean:
	rm -f $(OBJS) $(PRGM) 

explain:
	@echo "The information represents in the program:"
	@echo "Final executable name: $(PRGM)"
	@echo "Source files: $(SRCS)"
	@echo "Object files: $(OBJS)"

depend:$(DEPS)
	@echo "Dependencies are now up-to-date"



