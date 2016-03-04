PRGM = dwk++

CPPFLAGS =   -g
LDFLAGS = 
CFLAGS = -fopenmp  -framework OpenCL -O2 -ftree-vectorize -fopt-info-vec 


LIBCONFIG = -I${LIBCONFIG_DIR}/include  -L${LIBCONFIG_DIR}/lib -lconfig -lconfig++
NETCDF = -I${NETCDF_DIR}/include -L${NETCDF_DIR}/lib -lnetcdf
COMPILE = $(CXX) $(CPPFLAGS) $(CFLAGS) $(LIBCONFIG) $(NETCDF) -c 

LINKCC = $(CXX) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS)  $(LIBCONFIG) $(NETCDF) 

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



