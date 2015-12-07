PRGM = dwk++

CC = g++

CPPFLAGS =
LDFLAGS = 
CFLAGS = -Wall -ggdb3 

LIBCONFIG_DIR=/usr/local
LIBCONFIG = -I${LIBCONFIG_DIR}/include  -L${LIBCONFIG_DIR}/lib -lconfig -lconfig++
COMPILE = $(CC) $(CPPFLAGS) $(CFLAGS)  -c 

LINKCC = $(CC) $(LDFLAGS) $(BOOST) $(LIBCONFIG)

LIBA = 

SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

all :  $(PRGM)





$(PRGM):$(OBJS)
	$(LINKCC) $(OBJS) $(LIBA) $(HDF5L)  -o $(PRGM)

%.o:%.cpp
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



