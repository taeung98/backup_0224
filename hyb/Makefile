CC = mpicc 
RM = rm -rf 

GSL_DIR = /home/taeung/GSL
#/opt/gsl/gcc-4.8.5/2.7.1 
#OMP_DIR = /opt/mpi/intel-2021.4.0/openmpi-4.1.0
OMP_DIR = /opt/intel/oneapi/mpi/2021.4.0/lib/release
HDF_DIR = /opt/hdf5/intel-2021.4.0/1.12.0
 
CFLAGS = -g -O2 -Wall -mcmodel=medium -std=c99 -I$(GSL_DIR)/include -I$(HDF_DIR)/include -I$(OMP_DIR)/include  
LDFLAGS = -L$(GSL_DIR)/lib  -L$(HDF_DIR)/lib -L$(OMP_DIR)/lib
LINKS = -lm -lgsl -lgslcblas -lhdf5 -lmpi

OBJS = Lcoeff.o header_hyb.o hyb.o test.o pred_opt.o 

TARGETS = Lcoeff header_hyb hyb test pred_opt 
 
.PHONY: all clean dep 
.SUFFIXES : .c .o

.c .o :
	$(CC) $(CFLAGS) -c $<

all : $(TARGETS)

Lcoeff : Lcoeff.o
	$(CC) -o $@ Lcoeff.o $(LDFLAGS) $(LINKS)

hyb : hyb.o header_hyb.o
	$(CC) -o $@ hyb.o header_hyb.o $(LDFLAGS) $(LINKS)

test: test.o
	$(CC) -o $@ test.o $(LDFLAGS) $(LINKS)

pred_opt: pred_opt.o
	$(CC) -o $@ pred_opt.o $(LDFLAGS) $(LINKS)

clean :
	$(RM) *.o
	$(RM) $(TARGETS)

dep :
	$(CC) $(CFLAGS) -M $(OBJS:.o=c)

#hyb.o: $($@:.o=.c)	hyb.c	Makefile
