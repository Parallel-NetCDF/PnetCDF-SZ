/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <sz.h>
#include <mpi.h>
#include <pnetcdf.h>

#define ERR {if(err!=NC_NOERR){printf("Error at line=%d: %s\n", __LINE__, ncmpi_strerror(err));goto fn_exit;}}

static int
nc2SZtype(nc_type xtype)
{
    switch(xtype){
        case NC_CHAR :   return SZ_INT8;
        case NC_BYTE :   return SZ_INT8;
        case NC_SHORT :  return SZ_INT16;
        case NC_INT :    return SZ_INT32;
        case NC_FLOAT :  return SZ_FLOAT;
        case NC_DOUBLE : return SZ_DOUBLE;
        case NC_UBYTE :  return SZ_UINT8;
        case NC_USHORT : return SZ_UINT16;
        case NC_UINT :   return SZ_UINT32;
        case NC_INT64 :  return SZ_INT64;
        case NC_UINT64 : return SZ_UINT64;
        default:         return -1;
    }
}

static MPI_Datatype
nc2mpitype(nc_type xtype)
{
    switch(xtype){
        case NC_CHAR :   return MPI_CHAR;
        case NC_BYTE :   return MPI_SIGNED_CHAR;
        case NC_SHORT :  return MPI_SHORT;
        case NC_INT :    return MPI_INT;
        case NC_FLOAT :  return MPI_FLOAT;
        case NC_DOUBLE : return MPI_DOUBLE;
        case NC_UBYTE :  return MPI_UNSIGNED_CHAR;
        case NC_USHORT : return MPI_UNSIGNED_SHORT;
        case NC_UINT :   return MPI_UNSIGNED;
        case NC_INT64 :  return MPI_LONG_LONG_INT;
        case NC_UINT64 : return MPI_UNSIGNED_LONG_LONG;
        default:         return MPI_DATATYPE_NULL;
    }
}

static int
xlen_nc_type(nc_type xtype, int *size)
{
    switch(xtype) {
        case NC_BYTE:
        case NC_CHAR:
        case NC_UBYTE:  *size = 1; return NC_NOERR;
        case NC_SHORT:
        case NC_USHORT: *size = 2; return NC_NOERR;
        case NC_INT:
        case NC_UINT:
        case NC_FLOAT:  *size = 4; return NC_NOERR;
        case NC_DOUBLE:
        case NC_INT64:
        case NC_UINT64: *size = 8; return NC_NOERR;
        default: return NC_EBADTYPE;
    }
}

struct fspec {
    int nlvars;    /* Number of variables specified with -v
                    * option on command line */
    char** lvars;  /* list of variable names specified with -v
                    * option on command line */
};

static int 
var_decompress(MPI_Comm comm,
               int      in_ncid,  /* ID of input file */
               int      in_varid, /* ID of input variable */
               int      out_ncid) /* ID of output file */
{
    int j, k, err, nprocs, rank;
    int ndims, *dimids, dimid, out_varid, nblocks, *block_lens;
    nc_type xtype;
    MPI_Offset var_size, *starts, *counts;

    MPI_Comm_size(comm, &nprocs); /* number of MPI processes */
    MPI_Comm_rank(comm, &rank);   /* MPI process rank */

    /* obtain original variable size */
    char var_name[1024], dim_name[1024];
    err = ncmpi_inq_varname(in_ncid, in_varid, var_name); ERR
    sprintf(dim_name, "SZ.%s", var_name);
    err = ncmpi_inq_dimid(in_ncid, dim_name, &dimid); ERR
    if (err == NC_EBADDIM) {
        printf("Error: missing compressed size attribute for variable %s\n",var_name);
        return err;
    } ERR
    err = ncmpi_inq_dimlen(in_ncid, dimid, &var_size); ERR

    /* obtain original variable's nc_type */
    err = ncmpi_get_att(in_ncid, in_varid, "SZ.nc_type", &xtype);
    if (err == NC_ENOTATT) {
        printf("Error: missing original nc_type attribute for variable %s\n",var_name);
        return err;
    } ERR

    /* obtain original number of dimensions */
    err = ncmpi_get_att(in_ncid, in_varid, "SZ.ndims", &ndims);
    if (err == NC_ENOTATT) {
        printf("Error: missing original no. dimensions attribute for variable %s\n",var_name);
        return err;
    } ERR

    /* obtain original dimension IDs */
    dimids = (int*) malloc(ndims * sizeof(int));
    err = ncmpi_get_att(in_ncid, in_varid, "SZ.dimids", dimids);
    if (err == NC_ENOTATT) {
        printf("Error: missing original dimension IDs attribute for variable %s\n",var_name);
        return err;
    } ERR

    /* obtain number of compressed blocks */
    err = ncmpi_get_att(in_ncid, in_varid, "SZ.nblocks", &nblocks);
    if (err == NC_ENOTATT) {
        printf("Error: missing no. compressed blocks attribute for variable %s\n",var_name);
        return err;
    } ERR

    /* obtain sizes of compressed blocks */
    block_lens = (int*) malloc(nblocks * sizeof(int));
    err = ncmpi_get_att(in_ncid, in_varid, "SZ.block_lens", block_lens);
    if (err == NC_ENOTATT) {
        printf("Error: missing compressed block sizes attribute for variable %s\n",var_name);
        return err;
    } ERR

    /* obtain original put starts of the compressed blocks */
    starts = (MPI_Offset*) malloc(2 * nblocks * ndims * sizeof(MPI_Offset));
    err = ncmpi_get_att(in_ncid, in_varid, "SZ.starts", starts);
    if (err == NC_ENOTATT) {
        printf("Error: missing compressed block starts attribute for variable %s\n",var_name);
        return err;
    } ERR

    /* obtain original put counts of the compressed blocks */
    counts = starts + nblocks * ndims;
    err = ncmpi_get_att(in_ncid, in_varid, "SZ.counts", counts);
    if (err == NC_ENOTATT) {
        printf("Error: missing compressed block counts attribute for variable %s\n",var_name);
        return err;
    } ERR

    err = ncmpi_redef(out_ncid); ERR

    /* define a new variable in decompressed form in out_ncid */
    err = ncmpi_def_var(out_ncid, var_name, xtype, ndims, dimids, &out_varid); ERR

    /* copy variable's attributes, exclude SZ attributes */
    int natts;
    err = ncmpi_inq_varnatts(in_ncid, in_varid, &natts); ERR
    for (j=0; j<natts; j++) {
        char att_name[256];
        err = ncmpi_inq_attname(in_ncid, in_varid, j, att_name); ERR
        if (!strncmp("SZ.", att_name, 3)) continue;
        err = ncmpi_copy_att(in_ncid, in_varid, att_name, out_ncid, out_varid); ERR
    }

    /* done with updating metadata for out_ncid, exit define mode */
    err = ncmpi_enddef(out_ncid); ERR

    /* prepare to read variable from in_ncid. First, calculate the read amount
     * and starting indices for each process */
    int local_nblocks = nblocks / nprocs;
    int local_blockID = local_nblocks * rank;
    if (rank < nblocks % nprocs) {
        local_blockID += rank;
        local_nblocks++;
    }
    else {
        local_blockID += nblocks % nprocs;
    }

    /* allocate read buffer */
    signed char **inBuf = (signed char**) malloc(local_nblocks * sizeof(signed char*));

    /* calculate subarray start indices and counts */
    MPI_Offset start, count;
    start = 0;
    for (j=0; j<nblocks; j++) {
        if (j == local_blockID) break;
        start += block_lens[j];
    }

    /* read variable using nonblocking APIs, because number of subarrays to
     * read may be different among processes */
    int start_index = j;
    for (k=0; k<local_nblocks; k++) {
        count = block_lens[j];
        inBuf[k] = (signed char*) malloc(count);
        err = ncmpi_iget_vara_schar(in_ncid, in_varid, &start, &count, inBuf[k], NULL); ERR
        start += block_lens[j++];
    }
    err = ncmpi_wait_all(in_ncid, NC_REQ_ALL, NULL, NULL); ERR

    /* decompress read buffer */
    signed char **outBuf = (signed char**) malloc(local_nblocks * sizeof(signed char*));
    for (k=0; k<local_nblocks; k++) {
        size_t r[5];

        for (j=0; j<ndims; j++) r[j] = counts[(start_index+k+1)*ndims - j - 1];
        for (; j<5; j++) r[j] = 0;
        outBuf[k] = SZ_decompress(nc2SZtype(xtype), (unsigned char*)inBuf[k], block_lens[start_index+k], r[4], r[3], r[2], r[1], r[0]);
    }

    /* write decompressed buffer to file out_ncid */
    for (k=0; k<local_nblocks; k++) {
        int index = (start_index + k) * ndims;
        for (count=1,j=0; j<ndims; j++) count *= counts[index+j];
        err = ncmpi_iput_vara(out_ncid, out_varid, starts+index, counts+index, outBuf[k], count, nc2mpitype(xtype), NULL); ERR
    }
    err = ncmpi_wait_all(out_ncid, NC_REQ_ALL, NULL, NULL); ERR

    free(block_lens);
    free(dimids);
    free(starts);
    for (k=0; k<local_nblocks; k++) free(inBuf[k]);
    free(inBuf);
    free(outBuf);
fn_exit:
    return err;
}

static int
var_compress(MPI_Comm comm,
             int      in_ncid,  /* ID of input file */
             int      in_varid, /* ID of input variable */
             int      out_ncid) /* ID of output file */
{
    void *buf, *outbuf;
    int j, err, nprocs, rank, ndims, *dimids, dimid, out_varid, el_size;
    nc_type xtype;
    MPI_Offset len, *dimlen, *start, *count;

    MPI_Comm_size(comm, &nprocs); /* number of MPI processes */
    MPI_Comm_rank(comm, &rank);   /* MPI process rank */

    /* obtain variable size from in_ncid */
    err = ncmpi_inq_vartype(in_ncid, in_varid, &xtype); ERR
    err = ncmpi_inq_varndims(in_ncid, in_varid, &ndims); ERR
    if (ndims > 5) {
        char var_name[1024];
        err = ncmpi_inq_varname(in_ncid, in_varid, var_name); ERR
        printf("Error: no. dimensions of variable %s exceeds SZ limit (5)\n",var_name);
        printf("Error: skip variable %s\n", var_name);
        return err;
    }

    /* obtain variable dimension IDs and lengths */
    dimids = (int*) malloc(ndims * sizeof(int));
    dimlen = (MPI_Offset*) malloc(ndims * sizeof(MPI_Offset));
    err = ncmpi_inq_vardimid(in_ncid, in_varid, dimids); ERR
    err = xlen_nc_type(xtype, &el_size); ERR
    for (j=0; j<ndims; j++)
        err = ncmpi_inq_dimlen(in_ncid, dimids[j], &dimlen[j]); ERR

    /* partition the variable along the most significant dimension */
    start = (MPI_Offset*) malloc(2 * ndims * sizeof(MPI_Offset));
    count = start + ndims;
    count[0] = dimlen[0] / nprocs;
    start[0] = count[0] * rank;
    if (rank < dimlen[0] % nprocs) {
        start[0] += rank;
        count[0]++;
    }
    else {
        start[0] += dimlen[0] % nprocs;
    }
    len = count[0];
    for (j=1; j<ndims; j++) {
        start[j] = 0;
        count[j] = dimlen[j];
        len *= count[j];
    }
    free(dimlen);

    /* allocate read buffer */
    buf = malloc(len*el_size); if (buf == NULL) err = NC_ENOMEM; ERR

    /* read the entire variable in parallel */
    err = ncmpi_get_vara_all(in_ncid, in_varid, start, count, buf, len, nc2mpitype(xtype)); ERR

    /* compress read buffer */
    size_t outSize=0, r[5];
    for (j=0; j<ndims; j++) r[j] = count[ndims-j-1];
    for (; j<5; j++) r[j] = 0;
    if (len > 0)
        outbuf = SZ_compress(nc2SZtype(xtype), buf, &outSize, r[4], r[3], r[2], r[1], r[0]);

    /* gather compressed sizes from all processes */
    int *block_lens, local_size=outSize;
    block_lens = (int*) malloc(nprocs * sizeof(int));
    MPI_Allgather(&local_size, 1, MPI_INT, block_lens, 1, MPI_INT, comm);

    /* calculate the concatenated variable size and starting write offset
     * for this process */
    MPI_Offset offset=0, var_size=0;
    for (j=0; j<nprocs; j++) var_size += block_lens[j];
    for (j=0; j<nprocs; j++) {
        if (j == rank) break;
        offset += block_lens[j];
    }

    /* reenter define mode for out_ncid to add new metadata */
    err = ncmpi_redef(out_ncid); ERR

    /* define a new dimension, also size of new variable */
    char var_name[256], dim_name[256];
    err = ncmpi_inq_varname(in_ncid, in_varid, var_name); ERR
    sprintf(dim_name, "SZ.%s", var_name);
    err = ncmpi_def_dim(out_ncid, dim_name, var_size, &dimid); ERR

    /* define 1D variable of type NC_BYTE */
    err = ncmpi_def_var(out_ncid, var_name, NC_BYTE, 1, &dimid, &out_varid); ERR

    /* copy variable's attributes from in_ncid to out_ncid */
    int natts;
    err = ncmpi_inq_varnatts(in_ncid, in_varid, &natts); ERR
    for (j=0; j<natts; j++) {
        char att_name[256];
        err = ncmpi_inq_attname(in_ncid, in_varid, j, att_name); ERR
        err = ncmpi_copy_att(in_ncid, in_varid, att_name, out_ncid, out_varid); ERR
    }

    /* save original data type as special attributes in out_ncid */
    err = ncmpi_put_att_int(out_ncid, out_varid, "SZ.nc_type", NC_INT, 1, &xtype); ERR
    /* save original number of dimensions */
    err = ncmpi_put_att_int(out_ncid, out_varid, "SZ.ndims", NC_INT, 1, &ndims); ERR
    /* save original dimension IDs */
    err = ncmpi_put_att_int(out_ncid, out_varid, "SZ.dimids", NC_INT, ndims, dimids); ERR
    /* save number of compressed blocks */
    err = ncmpi_put_att_int(out_ncid, out_varid, "SZ.nblocks", NC_INT, 1, &nprocs); ERR
    /* save sizes of compressed blocks */
    err = ncmpi_put_att_int(out_ncid, out_varid, "SZ.block_lens", NC_INT, nprocs, block_lens); ERR

    /* collect starts[] and counts[] from all processes */
    MPI_Offset *starts, *counts;
    starts = (MPI_Offset*) malloc(2 * nprocs * ndims * sizeof(MPI_Offset));
    counts = starts + nprocs * ndims;
    MPI_Allgather(start, ndims, MPI_OFFSET, starts, ndims, MPI_OFFSET, comm);
    MPI_Allgather(count, ndims, MPI_OFFSET, counts, ndims, MPI_OFFSET, comm);

    /* save original subarray start[] for each block */
    err = ncmpi_put_att_longlong(out_ncid, out_varid, "SZ.starts", NC_INT64, ndims*nprocs, starts); ERR
    /* save original subarray count[] for each block */
    err = ncmpi_put_att_longlong(out_ncid, out_varid, "SZ.counts", NC_INT64, ndims*nprocs, counts); ERR
    err = ncmpi_enddef(out_ncid); ERR

    /* write variable in parallel */
    start[0] = offset;
    count[0] = outSize;
    err = ncmpi_put_vara_schar_all(out_ncid, out_varid, start, count, outbuf); ERR

    free(starts);
    free(block_lens);
    free(dimids);
    free(start);
    free(outbuf);
    free(buf);

fn_exit:
    return err;
}

static void
make_lvars(char *optarg, struct fspec* fspecp)
{
    char *cp = optarg;
    int nvars = 1;
    char ** cpp;

    /* compute number of variable names in comma-delimited list */
    fspecp->nlvars = 1;
    while (*cp++)
        if (*cp == ',')
            nvars++;

    fspecp->lvars = (char **) malloc(nvars * sizeof(char*));

    cpp = fspecp->lvars;
    /* copy variable names into list */
    for (cp = strtok(optarg, ",");
         cp != NULL;
         cp = strtok((char *) NULL, ",")) {

        *cpp = (char *) malloc(strlen(cp) + 1);
        strcpy(*cpp, cp);
        cpp++;
    }
    fspecp->nlvars = nvars;
}

static void
usage(char *cmd)
{
    char *help =
"Usage: %s [-h] | [-z sz.conf] [-v var1[,...]] input_file\n"
"       [-h]            Print help\n"
"       [-d]            Perform data decompression\n"
"       [-z sz.conf]    Input SZ configure file\n"
"       [-v var1[,...]] Output for variable(s) <var1>,... only\n"
"                       Without this option, all variables are compressed\n"
"       input_file      Input netCDF file name\n"
"*Parallel netCDF library version PNETCDF_RELEASE_VERSION of PNETCDF_RELEASE_DATE\n";
    fprintf(stderr, help, cmd);
}

int main(int argc, char** argv)
{
    extern int optind;
    int i, j, err=NC_NOERR, nerrs=0, opt, rank, nprocs, decompress=0;
    int in_ncid, out_ncid, cmode, ngatts, ndims, nvars, varid;
    int *in_varids=NULL;
    char outfile[1024], *cmd, cfgFile[1024], *infile;
    struct fspec *fspecp=NULL;
    MPI_Info info=MPI_INFO_NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    fspecp = (struct fspec*) calloc(1, sizeof(struct fspec));
    cmd = (char*) malloc(strlen(argv[0])+1);
    strcpy(cmd, argv[0]);
    cfgFile[0] = '\0';

    /* get command-line arguments */
    while ((opt = getopt(argc, argv, "dv:z:h")) != EOF) {
        switch(opt) {
            case 'd': decompress = 1;
                      break;
            case 'v': make_lvars(optarg, fspecp);
                      break;
            case 'z': strcpy(cfgFile, optarg);
                      break;
            case 'h':
            default:  usage(cmd);
                      free(fspecp);
                      goto fn_exit;
        }
    }
    argc -= optind;
    argv += optind;
    if (argc != 1) { /* input NetCDF file is required argument */
        fprintf(stderr, "Error: %s missing input file name\n", cmd);
        usage(cmd);
        for (i=0; i<fspecp->nlvars; i++)
            free(fspecp->lvars[i]);
        if (fspecp->lvars != NULL) free(fspecp->lvars);
        free(fspecp);
        free(cmd);
        nerrs++;
        goto fn_exit;
    }
    if (cfgFile[0] == '\0') { /* input cfg file is required argument */
        fprintf(stderr, "Error: %s missing SZ configure file\n", cmd);
        usage(cmd);
        for (i=0; i<fspecp->nlvars; i++)
            free(fspecp->lvars[i]);
        if (fspecp->lvars != NULL) free(fspecp->lvars);
        free(fspecp);
        free(cmd);
        nerrs++;
        goto fn_exit;
    }
    free(cmd);
    infile = argv[0]; /* required argument */

    /* initialize SZ */
    err = SZ_Init(cfgFile);
    if (err != SZ_SCES) {
        printf("Error in file %s line %d SZ_Init failed\n", __FILE__,__LINE__);
        nerrs++;
        goto fn_exit;
    }

    /* open input file */
    err = ncmpi_open(MPI_COMM_WORLD, infile, NC_NOWRITE, info, &in_ncid); ERR
    if (err != NC_NOERR) goto fn_exit;

    /* First check if all selected variables can be found in input file. */
    nvars = 0;
    if (fspecp->nlvars > 0) {
        in_varids = (int*) malloc(fspecp->nlvars * sizeof(int));
        for (i=0; i<fspecp->nlvars; i++) {
            err = ncmpi_inq_varid(in_ncid, fspecp->lvars[i], &varid);
            if (err == NC_ENOTVAR) {
                printf("Error: variable %s not found in %s\n",
                       fspecp->lvars[i], infile);
                continue;
            }
            else ERR

            if (decompress) { /* check if already compressed (SZ dimension) */
                char dim_name[1204];
                sprintf(dim_name, "SZ.%s",fspecp->lvars[i]);
                err = ncmpi_inq_dimid(in_ncid, dim_name, NULL);
                if (err == NC_EBADDIM) {
                    printf("Error: variable %s not compressed in %s\n",
                           fspecp->lvars[i], infile);
                    continue;
                }
                else ERR
            }
            in_varids[nvars++] = varid;
        }
    }
    else { /* -v is not given, compress/decompress all variables */
        err = ncmpi_inq_nvars(in_ncid, &nvars); ERR
        if (nvars > 0) {
            in_varids = (int*) malloc(nvars * sizeof(int));
            for (i=0; i<nvars; i++) in_varids[i] = i;
        }
    }

    if (decompress) { /* exclude variables that are not compressed */
        for (j=0,i=0; i<nvars; i++) {
            char var_name[1024], dim_name[1204];
            if (j < i) in_varids[j] = in_varids[i];
            err = ncmpi_inq_varname(in_ncid, in_varids[i], var_name); ERR
            sprintf(dim_name, "SZ.%s", var_name);
            err = ncmpi_inq_dimid(in_ncid, dim_name, NULL);
            if (err == NC_EBADDIM) {
                if (fspecp->nlvars > 0) /* variable explicitly requested */
                    printf("Warn: variable %s not compressed in %s\n",
                           var_name, infile);
                /* skip variable in_varids[i] */
                continue;
            }
            else if (err != NC_NOERR)
                goto fn_exit;
            j++;
        }
        nvars = j;
    }
    else { /* exclude variables that are already compressed */
        for (j=0,i=0; i<nvars; i++) {
            int ndims, dimid, SZ_dimid;
            char var_name[1024], dim_name[1204];
            if (j < i) in_varids[j] = in_varids[i];

            /* already compressed variables are 1D */
            err = ncmpi_inq_varndims(in_ncid, in_varids[i], &ndims); ERR
            if (ndims > 1) {j++; continue;}

            /* skip scalar variables */
            if (ndims == 0) continue;

            /* check variable dimension name against prefix "SZ." */
            err = ncmpi_inq_varname(in_ncid, in_varids[i], var_name); ERR
            sprintf(dim_name, "SZ.%s", var_name);
            err = ncmpi_inq_vardimid(in_ncid, in_varids[i], &dimid); ERR
            err = ncmpi_inq_dimid(in_ncid, dim_name, &SZ_dimid);
            if (err != NC_EBADDIM && dimid == SZ_dimid) {
                if (fspecp->nlvars > 0) /* variable explicitly requested */
                    printf("Warn: variable %s not compressed in %s\n",
                           var_name, infile);
                /* skip variable in_varids[i] */
                continue;
            }
            else if (err != NC_EBADDIM && err != NC_NOERR)
                goto fn_exit;
            j++;
        }
        nvars = j;
    }
    free(fspecp);
    if (nvars == 0) { /* no variables fit the request */
        if (decompress)
            printf("Error: no variable are compressed\n");
        else
            printf("Error: all variable are compressed\n");
        err = ncmpi_close(in_ncid); ERR
        nerrs++;
        goto fn_exit;
    }

    /* add extention to create file name for output file */
    if (decompress)
        sprintf(outfile, "%s.unsz", infile);
    else
        sprintf(outfile, "%s.sz", infile);

    /* create a new output file */
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(MPI_COMM_WORLD, outfile, cmode, info, &out_ncid); ERR

    /* copy the global attributes */
    err = ncmpi_inq_natts(in_ncid, &ngatts); ERR
    for (i=0; i<ngatts; i++) {
        char att_name[1024];
        err = ncmpi_inq_attname(in_ncid, NC_GLOBAL, i, att_name); ERR
        err = ncmpi_copy_att(in_ncid, NC_GLOBAL, att_name, out_ncid, NC_GLOBAL); ERR
    }
    /* copy all dimensions, except those with SZ prefix */
    err = ncmpi_inq_ndims(in_ncid, &ndims); ERR
    for (i=0; i<ndims; i++) {
        int dimid;
        char dimname[1024];
        MPI_Offset dimlen;

        err = ncmpi_inq_dimname(in_ncid, i, dimname); ERR
        if (decompress && !strncmp(dimname, "SZ.", 3)) continue;
        err = ncmpi_inq_dimlen(in_ncid, i, &dimlen); ERR
        err = ncmpi_def_dim(out_ncid, dimname, dimlen, &dimid); ERR
    }
    err = ncmpi__enddef(out_ncid, 8192, 1, 1, 1); ERR

    for (i=0; i<nvars; i++) { /* loop thru all variables */
        if (decompress)
            err = var_decompress(MPI_COMM_WORLD, in_ncid, in_varids[i], out_ncid);
        else
            err = var_compress(MPI_COMM_WORLD, in_ncid, in_varids[i], out_ncid);
        ERR
    }
    free(in_varids);

    err = ncmpi_close(in_ncid); ERR
    err = ncmpi_close(out_ncid); ERR

    SZ_Finalize();

fn_exit:
    MPI_Finalize();
    return (err != NC_NOERR || nerrs > 0);
}

