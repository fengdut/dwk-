PRGM = dwk++

CC = g++-5

CPPFLAGS =  -g
LDFLAGS = 
CFLAGS = -fopenmp 

LIBCONFIG_DIR=/usr/local

LIBCONFIG = -I${LIBCONFIG_DIR}/include  -L${LIBCONFIG_DIR}/lib -lconfig -lconfig++
COMPILE = $(CC) $(CPPFLAGS) $(CFLAGS) $(LIBCONFIG) -c 

LINKCC = $(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $(BOOST) $(LIBCONFIG)  

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



