MYTH_PATH = $(HOME)/opt/massivethreads

ifeq ($(origin CC), default)
CC := g++
# CC := icc
endif
CFLAGS  := -std=c++11 -march=native -O3 -g $(CFLAGS)
LDFLAGS := -lm

# CFLAGS += -O0 -g

# MassiveThreads
ifdef SPH_TASK_PARALLEL
CFLAGS  += -DSPH_TASK_PARALLEL -I$(MYTH_PATH)/include
LDFLAGS += -L$(MYTH_PATH)/lib -Wl,-R$(MYTH_PATH)/lib -lmyth
endif

# OpenMP
ifdef SPH_LOOP_PARALLEL
CFLAGS += -DSPH_LOOP_PARALLEL -fopenmp
endif

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.hpp)
EXEFILE = sph.out

all: $(CPPOBJS) $(CPPHDRS)
	$(CC) $(CFLAGS) $(CPPOBJS) -o $(EXEFILE) $(LDFLAGS)

%.o: %.cpp $(CPPHDRS)
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o *.out
