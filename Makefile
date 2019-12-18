.KEEP_STATE:

all: bin2nc pnc_sz


MPICC       = mpicc
PNETCDF_DIR = /home/sdi/Install/parallel-netcdf-1.9.0-install
SZ_DIR      = /home/sdi/Install/sz-2.1.8-install

DFLAGS      =
OPTFLAGS    = -g -Wall
INCFLAGS    = -I$(PNETCDF_DIR)/include -I$(SZ_DIR)/include
CFLAGS      = -O3 $(DFLAGS) $(INCFLAGS)

LDFLAGS     = -L$(PNETCDF_DIR)/lib -L$(SZ_DIR)/lib
LIBS        = -lpnetcdf -lSZ -lz -lzstd -lzlib -lm


bin2nc: bin2nc.c
	$(MPICC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LIBS)

pnc_sz: pnc_sz.c
	$(MPICC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LIBS)

clean:
	rm -rf core* *.o
	rm -rf bin2nc pnc_sz
