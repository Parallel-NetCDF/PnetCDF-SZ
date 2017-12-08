# This project is for developing SZ compression feature into PnetCDF library.

### Build utility program bin2nc. It converts a binary file to a NetCDF file
```
% make bin2nc
```
### Convert a binary file in Little Endian form to a NetCDF file
```
% ./bin2nc testdouble_8_8_128.dat
```
### Rename the file extension to ".nc"
```
% mv testdouble_8_8_128.dat.nc testdouble_8_8_128.nc
```
### Check file header of the newly generated NetCDF file
```
% ncdump -h testdouble_8_8_128.nc
netcdf testdouble_8_8_128 {
dimensions:
	Z = 128 ;
	Y = 8 ;
	X = 8 ;
variables:
	int scalar ;
		scalar:test = "dummy text attribute" ;
	double var3d(Z, Y, X) ;
		var3d:SZ\ test\ dataset = "x86/testdouble_8_8_128.dat" ;
		var3d:hidtory = "2017-12-08 12:36:51" ;

// global attributes:
		:SZ\ URL = "https://collab.cels.anl.gov/display/ESR/SZ" ;
		:Purpose = "This file is created to test integrated SZ compression and decompression in PnetCDF" ;
}
```
### Build utility program pnc_sz, a prototyped PnetCDF program with SZ compression and decompression feature.
```
% make pnc_sz
```
### The usage of pnc_sz can be obtained from command "./pnc_sz -h"
```
Usage: ./pnc_sz [-h] | [-z sz.conf] [-v var1[,...]] input_file
       [-h]            Print help
       [-v var1[,...]] Output for variable(s) <var1>,... only
                       Without this option, all variables are compressed
       [-z sz.conf]    Input SZ configure file
       input_file      Input netCDF file name
```

### Test run on 4 MPI processes to compress variables in input NetCDF file. Note scalar variables are not compressed.
```
% mpiexec -n 4 ./pnc_sz -z sz.config testdouble_8_8_128.nc
```
### A new file with file name extension ".sz" was created.  Check file header of the new NetCDF file.
```
% ncdump -h testdouble_8_8_128.nc.sz
netcdf testdouble_8_8_128.nc {
dimensions:
	Z = 128 ;
	Y = 8 ;
	X = 8 ;
	SZ.var3d = 18744 ;
variables:
	byte var3d(SZ.var3d) ;
		var3d:SZ\ test\ dataset = "x86/testdouble_8_8_128.dat" ;
		var3d:hidtory = "2017-12-08 12:36:51" ;
		var3d:SZ.nc_type = 6 ;
		var3d:SZ.ndims = 3 ;
		var3d:SZ.dimids = 0, 1, 2 ;
		var3d:SZ.nblocks = 4 ;
		var3d:SZ.block_lens = 4686, 4686, 4686, 4686 ;
		var3d:SZ.starts = 0LL, 0LL, 0LL, 32LL, 0LL, 0LL, 64LL, 0LL, 0LL, 96LL, 0LL, 0LL ;
		var3d:SZ.counts = 32LL, 8LL, 8LL, 32LL, 8LL, 8LL, 32LL, 8LL, 8LL, 32LL, 8LL, 8LL ;

// global attributes:
		:SZ\ URL = "https://collab.cels.anl.gov/display/ESR/SZ" ;
		:Purpose = "This file is created to test integrated SZ compression and decompression in PnetCDF" ;
}
```
### Test run on 3 MPI processes to decompress variables in input NetCDF file.
```
% mpiexec -n 3 ./pnc_sz -d -z sz.config testdouble_8_8_128.nc.sz
```
### A new file with file name extension ".unsz" was created.  Check file header of the new NetCDF file.
```
% ncdump -h testdouble_8_8_128.nc.sz.unsz
netcdf testdouble_8_8_128.nc.sz {
dimensions:
	Z = 128 ;
	Y = 8 ;
	X = 8 ;
variables:
	double var3d(Z, Y, X) ;
		var3d:SZ\ test\ dataset = "x86/testdouble_8_8_128.dat" ;
		var3d:hidtory = "2017-12-08 12:36:51" ;

// global attributes:
		:SZ\ URL = "https://collab.cels.anl.gov/display/ESR/SZ" ;
		:Purpose = "This file is created to test integrated SZ compression and decompression in PnetCDF" ;
}
```
Contact: @wkliao
