#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <cstring>

#include "stokes_mob.h"
#include "stokes_util.h"

#include "mob_debug.h"

using namespace stokesdt;

const struct option long_options[] = {
    {"npts",    1, NULL, 10},
    {"xyz",     1, NULL, 20},
    {"xi",      1, NULL, 30},
    {"rmax",    1, NULL, 40},
    {"dim",     1, NULL, 50},
    {"porder",  1, NULL, 60},
    {"ref",     1, NULL, 70},
    {NULL, 0, NULL, 0}
};

const char *const short_options = ":h";

static void usage (char *call)
{
    fprintf(stderr, "Usage: %s [OPTIONS]\n", call);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\t-h or --help         Display this information\n");
    fprintf(stderr, "\t--npts               Model file\n");
    fprintf(stderr, "\t--xyz                XYZ file\n");
    fprintf(stderr, "\t--xi                 Ewald paramter\n");    
    fprintf(stderr, "\t--rmax               Real-space cutoff\n");
    fprintf(stderr, "\t--dim                Dimension of FFT grid\n");
    fprintf(stderr, "\t--porder             Interpolation order\n");  
    fprintf(stderr, "\t--ref                Referenece file\n");  
}


int main(int argc, char **argv)
{
    char *xyz_file = NULL;
    char *ref_file = NULL;
    double xi = 0.5;
    double rmax = 4.0;
    int dim = 64;
    int porder = 4;

    /* parse arguments */
    int npos;
    int c = 0;
    while ((c = getopt_long(argc, argv, short_options,
                            long_options, NULL)) != -1) {
        switch (c) {
            case 'h':
                usage(argv[0]);
                return 0;
            case 10:
                npos = atoi(optarg);
                break;           
            case 20:
                xyz_file = strdup(optarg);
                break;             
            case 70:
                ref_file = strdup(optarg);
                break;             
            case 30:
                xi = atof(optarg);
                break;                
            case 40:
                rmax = atof(optarg);
                break;                
            case 50:
                dim = atoi(optarg);
                break;                
            case 60:
                porder = atoi(optarg);
                break;
            case ':':
                fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                return -1;
            case '?':
                fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                return -1;
            default:
                usage (argv[0]);
                return -1;
        }
    }
    if (xyz_file == NULL || ref_file == NULL) {
        usage(argv[0]);
        return -1;
    }
    
    fprintf(stdout, "\nTest SPME\n");
    fprintf(stdout, "----------\n");
    
    // input particles
    double *pos = (double *)malloc(sizeof(double) * npos * 3);
    double *rdi = (double *)malloc(sizeof(double) * npos);
    double r_min[3];
    double box_len[3];

    FILE *inf = fopen(xyz_file, "r");

    //  1. box size
    double dummy;
    for (int i = 0; i < 3; i++)
        fscanf(inf, "%lf,", &r_min[i]);
    fscanf(inf, "%lf\n", &dummy);
    for (int i = 0; i < 3; i++)
        fscanf(inf, "%lf,", &box_len[i]);                                                                      
    fscanf(inf, "%lf\n", &dummy);

    //  2. coordinates and radii
    for (int i = 0; i < npos; i++)
    {
        for (int j = 0; j < 3; j++)
            fscanf(inf, "%lf,", &pos[i*3+j]);
        fscanf(inf, "%lf\n", &rdi[i]);
    }
    fclose(inf);

    //  Assuming a cubic box
    double box_size = box_len[0];
    for (int i = 0; i < npos; i++)
    {
        pos[i*3]   -= r_min[0];
        pos[i*3+1] -= r_min[1];
        pos[i*3+2] -= r_min[2]; 
    }
    //printf("14th pts: (%lf, %lf, %lf, %lf)\n", pos[13*3+0], pos[13*3+1], pos[13*3+2], rdi[13]);
    //printf("512th pts: (%lf, %lf, %lf, %lf)\n", pos[511*3+0], pos[511*3+1], pos[511*3+2], rdi[511]);

    //  Matvec check
    int num_rhs  = 1;

    //  compute spme
    printf("SPME engine construction\n");
    MobSpme *mob_spme = new MobSpme(npos, rdi, box_size, xi, rmax, dim, porder);
    assert(mob_spme->Init());
    mob_spme->Update(pos, rdi);

    // read multiplicand and matvec ref results 
    int nm  = npos * 3;
    int ldm = nm;
    double *x   = (double *)malloc(sizeof(double) * ldm * num_rhs);
    double *ref = (double *)malloc(sizeof(double) * ldm * num_rhs);
    assert(x != NULL && ref != NULL);

    inf = fopen(ref_file, "r");
    for (int i = 0; i < num_rhs; i++)
    {
        for (int j = 0; j < 3*npos; j++)
            fscanf(inf, "%lf\n", &x[i*3*npos+j]);

        for (int j = 0; j < 3*npos; j++)
            fscanf(inf, "%lf\n", &ref[i*3*npos+j]);
    }
    fclose(inf);


    // Multiply vectors   
    double *v_spme = (double *)malloc(sizeof(double) * ldm * num_rhs);
    assert(v_spme != NULL);
    memset(v_spme, 0, sizeof(double) * ldm * num_rhs);
    double alpha = 1.0;
    double beta  = 0.0;
    fprintf(stdout, "x10 = %lf\n", x[10]);
    mob_spme->MulVector(num_rhs, alpha, ldm, x, beta, ldm, v_spme);
    fprintf(stdout, "x10 = %lf\n", x[10]);
   
    printf("SPME:  %lf, %lf, %lf\n%lf, %lf, %lf\n", v_spme[0], ref[0], v_spme[0]/ref[0], v_spme[1], ref[1], v_spme[1]/ref[1]); 

    // check results
    double error = 0.0;
    double vnorm = 0.0;
    for (int i = 0; i < 3 * npos * num_rhs; i++)
    {
        double diff = v_spme[i] - ref[i];
        error += diff * diff;
        vnorm += ref[i]*ref[i];
    }
    error = sqrt(error);
    vnorm = sqrt(vnorm);
    error = error / vnorm;
    fprintf(stdout, "SPME error = %le\n", error);

    if (1)
    {
        printf("\nEwald matrix construction\n");
        fprintf(stdout, "x10 = %lf\n", x[10]);
        MobEwald *mob_ewald = new MobEwald(npos, rdi, box_size, 1.0e-12);
        assert(mob_ewald->Init());
        mob_ewald->Update(pos, rdi);

        double *v_ewald = (double *)malloc(sizeof(double) * ldm * num_rhs);
        assert(v_ewald != NULL);
        memset(v_ewald, 0, sizeof(double) * ldm * num_rhs);
        mob_ewald->MulVector(num_rhs, alpha, ldm, x, beta, ldm, v_ewald);


        MobDebug *mob_debug = new MobDebug(npos, rdi, box_size, 1.0e-12, MobDebug::EWALD);
        assert(mob_debug->Init());
        double *v_ewald1 = (double *)malloc(sizeof(double) * ldm * num_rhs);
        assert(v_ewald1 != NULL);
        memset(v_ewald1, 0, sizeof(double) * ldm * num_rhs);
        mob_debug->MulVector(pos, rdi, num_rhs, alpha, ldm, x, beta, ldm, v_ewald1);

        fprintf(stdout, "x10 = %lf\n", x[10]);

        printf("Ewald: %lf, %lf, %lf\n%lf, %lf, %lf\n", v_ewald[0], ref[0], v_ewald[0]/ref[0], v_ewald[1], ref[1], v_ewald[1]/ref[1]); 

        // check results
        error = 0.0;
        vnorm = 0.0;
        double error1 = 0.0;
        double vnorm1 = 0.0;
        for (int i = 0; i < 3 * npos * num_rhs; i++)
        {
            double diff = v_ewald1[i] - ref[i];
            double diff1 = v_ewald[i] - v_ewald1[i];
            error += diff * diff;
            vnorm += ref[i]*ref[i];

            error1 += diff1 * diff1;
            vnorm1 += v_ewald1[i] * v_ewald1[i];
        }
        error = sqrt(error);
        vnorm = sqrt(vnorm);
        error = error / vnorm;

        error1 = sqrt(error1);
        vnorm1 = sqrt(vnorm1);
        error1 = error1 / vnorm1;

        fprintf(stdout, "Ewald-Ref error = %le\n", error);
        fprintf(stdout, "Ewald-SPME error = %le\n", error1);

        delete mob_ewald;
        delete mob_debug;
        free(v_ewald);
        free(v_ewald1);
    }

    /*
    FILE *log = fopen("./result.dat", "w");
    for (int i =0; i < 3 * npos * num_rhs; i++)
        fprintf(log, "%.12f\n", ref[i]);
    for (int i =0; i < 3 * npos * num_rhs; i++)
        fprintf(log, "%.12f\n", v_spme[i]);
    for (int i =0; i < 3 * npos * num_rhs; i++)
        fprintf(log, "%.12f\n", v_ewald[i]);
    fclose(log);
    */
    
    // clean up
    delete mob_spme;
    free(pos);
    free(rdi);
    free(x);
    free(ref);
    free(v_spme);
    free(xyz_file);
    free(ref_file);
    
    return 0;
}
