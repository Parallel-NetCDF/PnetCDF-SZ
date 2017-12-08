#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>

#include <mpi.h>
#include <pnetcdf.h>

#define ERR {if(err!=NC_NOERR){printf("Error at line=%d: %s\n", __LINE__, ncmpi_strerror(err));nerrs++;}}

#define NZ 128
#define NY   8
#define NX   8

int main(int argc, char** argv)
{
    int err, nerrs=0, fd, ncid, cmode, varid, dimid[3];
    char filename[256];
    double *buf;
    ssize_t len, expect;

    MPI_Init(&argc, &argv);
    if (argc != 2) {
        printf("Usage: %s input_file\n",argv[0]);
        MPI_Finalize();
        return 1;
    }

    time_t timer;
    char str_time[26];
    struct tm* tm_info;

    time(&timer);
    tm_info = localtime(&timer);
    strftime(str_time, 26, "%Y-%m-%d %H:%M:%S", tm_info);

    fd = open(argv[1], O_RDONLY, 0600);
    if (fd < 0) {
        printf("Error at line %d : open file %s (%s)\n", __LINE__, argv[1], strerror(errno));
        MPI_Finalize();
        return 1;
    }
    expect = NZ*NY*NX*sizeof(double);
    buf = (double*) malloc(expect);
    len = read(fd, buf, expect);
    if (len < 0) {
        printf("Error at line %d : read %s (%s)\n", __LINE__, argv[1], strerror(errno));
        MPI_Finalize();
        return 1;
    }
    else if (len != expect) {
        printf("Error at line %d : read amount expect %zd but got %zd\n",__LINE__,expect,len);
        MPI_Finalize();
        return 1;
    }
    if (0 != close(fd)) {
        printf("Error at line %d : close (%s)\n",__LINE__, strerror(errno));
    }

    sprintf(filename, "%s.nc", argv[1]);
    cmode = NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, &ncid); ERR
    char *attr="https://collab.cels.anl.gov/display/ESR/SZ";
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "SZ URL", strlen(attr), attr); ERR
    attr="This file is created to test integrated SZ compression and decompression in PnetCDF";
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "Purpose", strlen(attr), attr); ERR

    err = ncmpi_def_dim(ncid, "Z", NZ, &dimid[0]); ERR
    err = ncmpi_def_dim(ncid, "Y", NY, &dimid[1]); ERR
    err = ncmpi_def_dim(ncid, "X", NX, &dimid[2]); ERR

    err = ncmpi_def_var(ncid, "scalar", NC_INT, 0, NULL, &varid); ERR
    attr="dummy text attribute";
    err = ncmpi_put_att_text(ncid, varid, "test", strlen(attr), attr); ERR

    err = ncmpi_def_var(ncid, "var3d", NC_DOUBLE, 3, dimid, &varid); ERR
    attr="x86/testdouble_8_8_128.dat";
    err = ncmpi_put_att_text(ncid, varid, "SZ test dataset", strlen(attr), attr); ERR

    err = ncmpi_put_att_text(ncid, varid, "hidtory", strlen(str_time), str_time); ERR
    err = ncmpi_enddef(ncid); ERR
    err = ncmpi_put_var_double_all(ncid, varid, buf); ERR
    err = ncmpi_close(ncid); ERR
    free(buf);

    MPI_Finalize();
    return nerrs;
}

