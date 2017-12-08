.KEEP_STATE:

all:

MPICC       = mpicc
PNETCDF_DIR = $(HOME)/PnetCDF
SZ_DIR      = $(HOME)/SZ

DFLAGS      = 
OPTFLAGS    = -g -Wall
INCFLAGS    = -I$(PNETCDF_DIR)/include -I$(SZ_DIR)/include
CFLAGS      = $(OPTFLAGS) $(DFLAGS) $(INCFLAGS)

LDFLAGS     = -L$(PNETCDF_DIR)/lib -L$(SZ_DIR)/lib
LIBS        = -lpnetcdf -lSZ -lz -lm


bin2nc: bin2nc.c
	$(MPICC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LIBS)

pnc_sz: pnc_sz.c
	$(MPICC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LIBS)

clean:
	rm -rf core* *.o
	rm -rf bin2nc pnc_sz
