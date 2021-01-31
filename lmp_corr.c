#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <errno.h>
#include <glob.h>
#include <zlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "lmp_corr.h"

# define LINE_BUF_SIZE 1024 // maximum line length from input

static glob_t G_FILE_NAMES; // List of files to load
static input G_IN; // Input

// -------------------------------------------------------------------
// CALLER FOR THE OTHER FUNCTIONS IN THE FILE
// -------------------------------------------------------------------

// ------ Caller for the other functions in the file  ------
void analyze_lmp(input in) {

  // Set global variable with input from user
  G_IN = in;

  // Set global variable with file names
  if (strcmp(G_IN.fluct_file, "NO-INPUT")==0) get_file_names(in.config_file);
  
  // Set number of threads for parallel calculations
  omp_set_num_threads(G_IN.num_threads);

  // Compute intermediate scattering function
  if (in.isf) isf();

  // Compute longitudinal velocity correlation function
  if (in.lvcf) lvcf();
  
  // Free memory
  globfree(&G_FILE_NAMES);

}

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE INTERMEDIATE SCATTERING FUNCTION
// -------------------------------------------------------------------

// ------ Compute intermediate scattering function ------
void isf(){

  // Variables to store the information read from the file
  int n_files = 0; // Number of files
  int n_atoms = 0; // Number of atoms
  double *sim_box = NULL; // Simulation box
  double *xx = NULL, *yy = NULL, *zz = NULL; // Coordinates
  double *vx = NULL, *vy = NULL, *vz = NULL; // Velocities

  // Variables to compute the density fluctuations
  int nq = 0, nq_dir = 0; // Number of points in the wave-vector grid
  double dq = 0; // Resolution of the wave-vector grid
  double LL = 0; // Smallest simulation box size 
  double *qq = NULL; // Wave-vector grid
  double complex *drhok = NULL; // Density fluctuations
  double complex *drhomk = NULL; // Density fluctuations (complex conjugate)

  // Variables to compute the intermediate scattering function
  double complex *fkt = NULL;
  double norm_fact;

  if (strcmp(G_IN.fluct_file, "NO-INPUT")==0){ // Analyze LAMMPS configurations

    // Set number of files
    n_files = G_FILE_NAMES.gl_pathc;

    // Initialize
    init_corr(&n_atoms, &LL, &dq, &nq,
	      &nq_dir, &sim_box,
	      &xx, &yy, &zz,
	      &vx, &vy, &vz, &qq,
	      &drhok, &drhomk, &fkt);
    
    // Compute density fluctuations
    compute_fluct(n_atoms, LL, dq, nq, nq_dir, sim_box,
		  xx, yy, zz, vx, vy, vz, qq, 
		  drhok, drhomk, false);

  }
  else { // Load pre-computed fluctuations
    fluct_input(&n_atoms, &dq, &nq, &nq_dir, &n_files, 
		&qq, &drhok, &drhomk, &fkt, false); 
  }


  // Compute intermediate scattering function (averaged over directions)
  int arr_idx;
  for (int ii=0; ii<nq_dir; ii++){
    for (int jj=0; jj<nq; jj++){
      for (int kk=0; kk<n_files; kk++){
  	arr_idx = idx2(jj,kk,nq);
  	for (int ll=0; ll<n_files-kk; ll++){
  	  fkt[arr_idx] += drhomk[idx3(jj,ii,ll,nq,nq_dir)]*
  	                  drhok[idx3(jj,ii,kk+ll,nq,nq_dir)]/(nq_dir*n_atoms);
  	}
      }
    }
  }

  // Normalize intermediate scattering function
  for (int jj=0; jj<nq; jj++){
    for (int kk=0; kk<n_files; kk++){
      norm_fact = 0.0;
      arr_idx = idx2(jj,kk,nq);
      for (int ll=0; ll<n_files-kk; ll++){
  	norm_fact+= 1.0;
      }
      fkt[arr_idx] /= norm_fact;
    }
  }
  
  // Write output
  isf_output(fkt, dq, nq, n_files);

  // Free memory
  free_corr(sim_box, xx, yy, zz,
	    vx, vy, vz,
	    qq, drhok, drhomk,fkt);

}

// ------ Write output for isf calculations ------
void isf_output(double complex *fkt, double dq, int nq, int n_files){
  
  
  FILE *fid;
  
  // Real part
  fid = fopen("isf_real.dat", "w");
  fprintf(fid, "ITEM: NUMBER OF WAVE-VECTORS (k)\n");
  fprintf(fid, "%d\n", nq);
  fprintf(fid, "ITEM: NUMBER OF TIME-STEPS (t)\n");
  fprintf(fid, "%d\n", n_files);
  fprintf(fid, "ITEM: INTERMEDIATE SCATTERING FUNCTION: k, t=0, t=dt, t=2*dt, ...\n"); 
  for (int ii=0; ii<nq; ii++){
    fprintf(fid, "%.8f ", (ii+1)*dq);
    for (int jj=0; jj<n_files; jj++){
      fprintf(fid, "%.8f ", creal(fkt[idx2(ii,jj,nq)]));
    }
    fprintf(fid, "\n");
  }
  fclose(fid);

  // Imaginary part (should be close to zero)
  fid = fopen("isf_imag.dat", "w");
  fprintf(fid, "ITEM: NUMBER OF WAVE-VECTORS (k)\n");
  fprintf(fid, "%d\n", nq);
  fprintf(fid, "ITEM: NUMBER OF TIME-STEPS (t)\n");
  fprintf(fid, "%d\n", n_files);
  fprintf(fid, "ITEM: TIME-STEP (dt)\n");
  fprintf(fid, "%f\n", G_IN.dt);
  fprintf(fid, "ITEM: INTERMEDIATE SCATTERING FUNCTION: k, t=0, t=dt, t=2*dt, ...\n");
  for (int ii=0; ii<nq; ii++){
    fprintf(fid, "%.8f ", (ii+1)*dq);
    for (int jj=0; jj<n_files; jj++){
      fprintf(fid, "%.8f ", cimag(fkt[idx2(ii,jj,nq)]));
    }
    fprintf(fid, "\n");
  }
  fclose(fid);
  
}

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE LONGITUDINAL VELOCITY CORRELATIONS
// -------------------------------------------------------------------

// ------ Compute intermediate scattering function ------
void lvcf(){

  // Variables to store the information read from the file
  int n_files = 0; // Number of files
  int n_atoms = 0; // Number of atoms
  double *sim_box = NULL; // Simulation box
  double *xx = NULL, *yy = NULL, *zz = NULL; // Coordinates
  double *vx = NULL, *vy = NULL, *vz = NULL; // Velocities

  // Variables to compute the density fluctuations
  int nq = 0, nq_dir = 0; // Number of points in the wave-vector grid
  double dq = 0; // Resolution of the wave-vector grid
  double LL = 0; // Smallest simulation box size 
  double *qq = NULL; // Wave-vector grid
  double complex *dvk = NULL; // Velocity fluctuations
  double complex *dvmk = NULL; // Velocity fluctuations (complex conjugate)

  // Variables to compute the intermediate scattering function
  double complex *ckt = NULL;
  double norm_fact;


  if (strcmp(G_IN.fluct_file, "NO-INPUT")==0){ // Analyze LAMMPS configurations

    // Set number of files
    n_files = G_FILE_NAMES.gl_pathc;

    // Initialize
    init_corr(&n_atoms, &LL, &dq, &nq,
	      &nq_dir, &sim_box,
	      &xx, &yy, &zz,
	      &vx, &vy, &vz, &qq,
	      &dvk, &dvmk, &ckt);
    
    // Compute velocity fluctuations
    compute_fluct(n_atoms, LL, dq, nq, nq_dir, sim_box,
		  xx, yy, zz, vx, vy, vz, qq, 
		  dvk, dvmk, true);
  }
  else { // Load pre-computed fluctuations
    fluct_input(&n_atoms, &dq, &nq, &nq_dir, &n_files, 
		&qq, &dvk, &dvmk, &ckt, true); 
  }

  // Compute longitudinal velocity correlation function (averaged over directions)
  int arr_idx;
  for (int ii=0; ii<nq_dir; ii++){
    for (int jj=0; jj<nq; jj++){
      for (int kk=0; kk<n_files; kk++){
  	arr_idx = idx2(jj,kk,nq);
  	for (int ll=0; ll<n_files-kk; ll++){
  	  ckt[arr_idx] += dvmk[idx3(jj,ii,ll,nq,nq_dir)]*
  	                  dvk[idx3(jj,ii,kk+ll,nq,nq_dir)]/(nq_dir*n_atoms);
  	}
      }
    }
  }

  // Normalize longitudinal velocity correlation function
  for (int jj=0; jj<nq; jj++){
    for (int kk=0; kk<n_files; kk++){
      norm_fact = 0.0;
      arr_idx = idx2(jj,kk,nq);
      for (int ll=0; ll<n_files-kk; ll++){
  	norm_fact+= 1.0;
      }
      ckt[arr_idx] /= norm_fact;
    }
  }
  
  // Write output
  lvcf_output(ckt, dq, nq, n_files);

  // Free memory
  free_corr(sim_box, xx, yy, zz,
	    vx, vy, vz,
	    qq, dvk, dvmk, ckt);

}


// ------ Write output for lvcf calculations ------
void lvcf_output(double complex *ckt, double dq, int nq, int n_files){
  
  
  FILE *fid;
  
  // Real part
  fid = fopen("lvcf_real.dat", "w");
  fprintf(fid, "ITEM: NUMBER OF WAVE-VECTORS (k)\n");
  fprintf(fid, "%d\n", nq);
  fprintf(fid, "ITEM: NUMBER OF TIME-STEPS (t)\n");
  fprintf(fid, "%d\n", n_files);
  fprintf(fid, "ITEM: LONGITUDINAL VELOCITY CORRELATION FUNCTION: k, t=0, t=dt, t=2*dt, ...\n"); 
  for (int ii=0; ii<nq; ii++){
    fprintf(fid, "%.8f ", (ii+1)*dq);
    for (int jj=0; jj<n_files; jj++){
      fprintf(fid, "%.8f ", creal(ckt[idx2(ii,jj,nq)]));
    }
    fprintf(fid, "\n");
  }
  fclose(fid);

  // Imaginary part (should be zero)
  fid = fopen("lvcf_imag.dat", "w");
  fprintf(fid, "ITEM: NUMBER OF WAVE-VECTORS (k)\n");
  fprintf(fid, "%d\n", nq);
  fprintf(fid, "ITEM: NUMBER OF TIME-STEPS (t)\n");
  fprintf(fid, "%d\n", n_files);
  fprintf(fid, "ITEM: TIME-STEP (dt)\n");
  fprintf(fid, "%f\n", G_IN.dt);
  fprintf(fid, "ITEM: LONGITUDINAL VELOCITY CORRELATION FUNCTION: k, t=0, t=dt, t=2*dt, ...\n"); 
  for (int ii=0; ii<nq; ii++){
    fprintf(fid, "%.8f ", (ii+1)*dq);
    for (int jj=0; jj<n_files; jj++){
      fprintf(fid, "%.8f ", cimag(ckt[idx2(ii,jj,nq)]));
    }
    fprintf(fid, "\n");
  }
  fclose(fid);
  
}

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE FLUCTUATIONS 
// -------------------------------------------------------------------

// ------ Compute fluctuations ------
void compute_fluct(int n_atoms, double LL, double dq, int nq,
		   int nq_dir, double *sim_box,
		   double *xx, double *yy, double *zz,
		   double *vx, double *vy, double *vz, 
		   double *qq, double complex *dfk, double complex *dfmk,
		   bool vel){

  // Variables 
  int time_step = 0; 
  int n_files = G_FILE_NAMES.gl_pathc;
  int n_atoms_tmp;

  // Loop through the configuration files
  for (int ii=0; ii<n_files; ii++){
    
    // Read file content (it is assumed that one file contains one configuration)
    printf("%s\n",G_FILE_NAMES.gl_pathv[ii]);
    fflush(stdout);
    read_file(ii, &time_step, &n_atoms_tmp, sim_box,
  	      xx, yy, zz, vx, vy, vz);
    
    // check input consistency
    check_consistency(n_atoms, n_atoms_tmp, LL, sim_box);


    //Loop through the wave-vector magnitudes
    #pragma omp parallel for
    for (int jj=0; jj<nq; jj++){
      int arr_idx;
      double qx, qy, qz, qmag, qr, cosqr, sinqr, flong;
      // Loop through the wave-vector directions
      for (int kk=0; kk<nq_dir; kk++){
    	arr_idx = idx3(jj,kk,ii,nq,nq_dir);
	qx = qq[idx3(jj,kk,0,nq,nq_dir)];
	qy = qq[idx3(jj,kk,1,nq,nq_dir)];
	qz = qq[idx3(jj,kk,2,nq,nq_dir)];
	qmag = sqrt(qx*qx + qy*qy + qz*qz);
    	dfk[arr_idx] = 0.0 + I*0.0;
    	dfmk[arr_idx] = 0.0 + I*0.0;
    	// Loop through the atom positions
    	for (int mm=0; mm<n_atoms; mm++){
    	  qr = qx*xx[mm] + qy*yy[mm] + qz*zz[mm];
    	  cosqr = cos(qr);
    	  sinqr = sin(qr);
	  if (vel) flong = (qx*vx[mm] + qy*vy[mm] + qz*vz[mm])/qmag; // Velocity fluctuations
	  else flong = 1.0; // Density fluctuations
	  dfk[arr_idx] += flong*(cosqr + I * sinqr);
	  dfmk[arr_idx] += flong*(cosqr - I * sinqr);
    	}
      }
    }
  }

  // Write output
  fluct_output(dfk, dfmk, qq, dq, nq, n_files, nq_dir, n_atoms, vel);


}

// ------ Write output for fluctuations ------
void fluct_output(double complex *dfk, double complex *dfmk, 
		  double *qq, double dq, int nq, 
		  int n_files, int nq_dir, int n_atoms, bool vel){
  
  // Open binary file
  FILE *fid = NULL;
  if (vel) fid = fopen("vel_fluct.bin", "wb");
  else fid = fopen("dens_fluct.bin", "wb");
  if (fid == NULL) {
    fprintf(stderr,"Error while creating fluctuation file");
    exit(EXIT_FAILURE);
  }

  // Fluctuation identifier
  fwrite(&vel, sizeof(bool), 1, fid);

  // Number of wave-vectors magnitudes
  fwrite(&nq, sizeof(int), 1, fid);

  // Magnitude step
  fwrite(&dq, sizeof(double), 1, fid);

  // Number of wave-vectors directions
  fwrite(&nq_dir, sizeof(int), 1, fid);

  // Number of time-steps
  fwrite(&n_files, sizeof(int), 1, fid);

  // Time-step
  fwrite(&G_IN.dt, sizeof(double), 1, fid);

  // Number of atoms
  fwrite(&n_atoms, sizeof(int), 1, fid);

  // Wave-vectors
  fwrite(qq, sizeof(double complex), nq * nq_dir * 3, fid);

  // Fluctuations
  fwrite(dfk, sizeof(double complex), nq * nq_dir * n_files, fid);

  // Fluctuations (complex conjugate)
  fwrite(dfmk, sizeof(double complex), nq * nq_dir * n_files, fid);

  // Close binary file
  fclose(fid);

}


// ------ Read fluctuations from file ------
void fluct_input(int *out_n_atoms, double *out_dq, int *out_nq,  int *out_nq_dir, 
		 int *out_n_files, double **out_qq, 
		 double complex **out_dfk, double complex **out_dfmk, 
		 double complex **out_cf, bool vel){
  

  // Variables
  bool vel_r;
  int n_atoms, nq, n_files, nq_dir;
  double dq, *qq;
  double complex *dfk, *dfmk, *cf;

  // Open binary file
  FILE *fid = NULL;
  if (vel) fid = fopen("vel_fluct.bin", "rb");
  else fid = fopen("dens_fluct.bin", "rb");
  if (fid == NULL) {
    fprintf(stderr,"Error while opening fluctuation file\n");
    exit(EXIT_FAILURE);
  }

  // Fluctuation identifier
  fread(&vel_r, sizeof(bool), 1, fid);
  if (vel_r != vel){
    if (vel) fprintf(stderr,"ERROR: Attempt to load velocity fluctuations from file %s, but file %s contains only density fluctuations\n",
		    G_IN.fluct_file, G_IN.fluct_file);
    else fprintf(stderr,"ERROR: Attempt to load density fluctuations from file %s, but file %s contains only velocity fluctuations\n",
		    G_IN.fluct_file, G_IN.fluct_file);
    exit(EXIT_FAILURE);
  }
  
  // Number of wave-vectors magnitudes
  fread(&nq, sizeof(int), 1, fid);

  // Magnitude step
  fread(&dq, sizeof(double), 1, fid);

  // Number of wave-vectors directions
  fread(&nq_dir, sizeof(int), 1, fid);

  // Number of time-steps
  fread(&n_files, sizeof(int), 1, fid);

  // Time-step
  fread(&G_IN.dt, sizeof(double), 1, fid);

  // Number of atoms
  fread(&n_atoms, sizeof(int), 1, fid);

  // Allocate array for wave-vector grid
  qq = malloc(sizeof(double) * nq * nq_dir * 3);
  if (qq == NULL){
    fprintf(stderr,"ERROR: Failed allocation for wave-vector grid array (qq)\n");
    exit(EXIT_FAILURE);
  }

  // Allocate arrays for fluctuations
  dfk = malloc(sizeof(double complex) * nq * nq_dir * n_files);
  if (dfk == NULL){
    fprintf(stderr,"ERROR: Failed allocation for fluctuations array (dfk)\n");
    exit(EXIT_FAILURE);
  }
  dfmk = malloc(sizeof(double complex) * nq * nq_dir * n_files);
  if (dfmk == NULL){
    fprintf(stderr,"ERROR: Failed allocation for fluctuations array (dfmk) \n");
    exit(EXIT_FAILURE);
  }

  // Allocate array for correlation function
  cf = malloc(sizeof(double complex) * nq * n_files);
  if (cf == NULL){
    fprintf(stderr,"ERROR: Failed allocation for correlation function\n");
    exit(EXIT_FAILURE);
  }


  // Wave-vectors
  fread(qq, sizeof(double complex), nq * nq_dir * 3, fid);

  // Fluctuations
  fread(dfk, sizeof(double complex), nq * nq_dir * n_files, fid);

  // Fluctuations (complex conjugate)
  fread(dfmk, sizeof(double complex), nq * nq_dir * n_files, fid);

  // Close binary file
  fclose(fid);

  // Output
  *out_n_atoms = n_atoms;
  *out_nq = nq;
  *out_dq = dq;
  *out_nq_dir = nq_dir;
  *out_n_files = n_files;
  *out_qq = qq;
  *out_dfk = dfk;
  *out_dfmk = dfmk;
  *out_cf = cf;

}


// -------------------------------------------------------------------
// FUNCTIONS USED TO INITIALIZE AND FREE MEMORY
// -------------------------------------------------------------------

// ------ Initialize calculation for correlation functions ------
void init_corr(int *out_n_atoms, double *out_LL, double *out_dq, int *out_nq,
	       int *out_nq_dir, double **out_sim_box, 
	       double **out_xx, double **out_yy, double **out_zz,
	       double **out_vx, double **out_vy, double **out_vz,
	       double **out_qq, double complex **out_dfk,
	       double complex **out_dfmk, double complex **out_cf){

  // Variables 
  int n_atoms = 0, nq = 0, n_files = G_FILE_NAMES.gl_pathc, 
    nq_dir = 0, idx_dir = 0;
  double LL, dq, dtheta, dphi, qtmp, qsint, qcost;
  double *sim_box = NULL, *xx = NULL, *yy = NULL, *zz = NULL,
    *vx = NULL, *vy = NULL, *vz = NULL, *qq = NULL;
  double complex *dfk = NULL, *dfmk = NULL, *cf = NULL;
  
  // Read the first configuration file and allocate the necessary arrays
  // Note: it is assumed that all files refer to the same NVT 
  // simulation. i.e. the number of atoms and the simulation box, 
  // remain constant while reading the files. If this is 
  // not the case (which could happen for NpT or for grand-canonical
  // simulations) the code will stop and produce an error
  read_file_init(&n_atoms, &LL);

  //Wave-vector grid
  dq = 2.0*M_PI/LL;
  nq = (int)G_IN.q_max/dq;
  dtheta = M_PI/(G_IN.num_theta-1);
  dphi = 2.0*M_PI/G_IN.num_phi;
  if (G_IN.num_theta==1)
    nq_dir = 1;
  else
    nq_dir = G_IN.num_phi*(G_IN.num_theta-2) + 2;

  // Allocate array to store the simulation box information
  sim_box = malloc(sizeof(double) * 3);
  if (sim_box == NULL){
    fprintf(stderr,"ERROR: Failed allocation for simulation box array\n");
    exit(EXIT_FAILURE);
  }

  // Allocate arrays to store the configuration
  xx = malloc(sizeof(double) * n_atoms);
  if (xx == NULL){
    fprintf(stderr,"ERROR: Failed allocation for configuration array (xx)\n");
    exit(EXIT_FAILURE);
  }
  yy = malloc(sizeof(double) * n_atoms);
  if (yy == NULL){
    fprintf(stderr,"ERROR: Failed allocation for configuration array (yy)\n");
    exit(EXIT_FAILURE);
  }
  zz = malloc(sizeof(double) * n_atoms);
  if (zz == NULL){
    fprintf(stderr,"ERROR: Failed allocation for configuration array (zz)\n");
    exit(EXIT_FAILURE);
  }
  vx = malloc(sizeof(double) * n_atoms);
  if (vx == NULL){
    fprintf(stderr,"ERROR: Failed allocation for configuration array (vx)\n");
    exit(EXIT_FAILURE);
  }
  vy = malloc(sizeof(double) * n_atoms);
  if (vy == NULL){
    fprintf(stderr,"ERROR: Failed allocation for configuration array (vy)\n");
    exit(EXIT_FAILURE);
  }
  vz = malloc(sizeof(double) * n_atoms);
  if (vz == NULL){
    fprintf(stderr,"ERROR: Failed allocation for configuration array (vz)\n");
    exit(EXIT_FAILURE);
  }

  // Allocate array for wave-vector grid
  qq = malloc(sizeof(double) * nq * nq_dir * 3);
  if (qq == NULL){
    fprintf(stderr,"ERROR: Failed allocation for wave-vector grid array (qq)\n");
    exit(EXIT_FAILURE);
  }

  // Allocate arrays for fluctuations
  dfk = malloc(sizeof(double complex) * nq * nq_dir * n_files);
  if (dfk == NULL){
    fprintf(stderr,"ERROR: Failed allocation for fluctuations array (dfk)\n");
    exit(EXIT_FAILURE);
  }
  dfmk = malloc(sizeof(double complex) * nq * nq_dir * n_files);
  if (dfmk == NULL){
    fprintf(stderr,"ERROR: Failed allocation for fluctuations array (dfmk) \n");
    exit(EXIT_FAILURE);
  }

  // Allocate array for correlation function
  cf = malloc(sizeof(double complex) * nq * n_files);
  if (cf == NULL){
    fprintf(stderr,"ERROR: Failed allocation for correlation function\n");
    exit(EXIT_FAILURE);
  }

  // Initialize the wave-vector grid
  for (int ii=0; ii<nq; ii++){
    qtmp = (ii+1)*dq;
    // theta = 0 
    qq[idx3(ii, 0, 0, nq, nq_dir)] = 0;
    qq[idx3(ii, 0, 1, nq, nq_dir)] = 0;
    qq[idx3(ii, 0, 2, nq, nq_dir)] = qtmp;
    /* qq[idx3(ii, 0, 0, nq, nq_dir)] = 0; */
    /* qq[idx3(ii, 0, 1, nq, nq_dir)] = qtmp/sqrt(2.0); */
    /* qq[idx3(ii, 0, 2, nq, nq_dir)] = qtmp/sqrt(2.0); */
    /* if (ii==0) printf("%f %f %f %f %f\n",qq[idx3(ii, 0, 0, nq, nq_dir)],qq[idx3(ii, 0, 1, nq, nq_dir)],qq[idx3(ii, 0, 2, nq, nq_dir)], */
    /* 		      qtmp,sqrt(qq[idx3(ii, 0, 0, nq, nq_dir)]*qq[idx3(ii, 0, 0, nq, nq_dir)] +  */
    /* 				qq[idx3(ii, 0, 1, nq, nq_dir)]*qq[idx3(ii, 0, 1, nq, nq_dir)] +  */
    /* 				qq[idx3(ii, 0, 2, nq, nq_dir)]*qq[idx3(ii, 0, 2, nq, nq_dir)])); */
    // theta > 0 (if necessary)
    if (nq_dir > 1) {
      for (int jj=0; jj<G_IN.num_theta-2; jj++){
	qsint = qtmp * sin(dtheta*(jj+1));
	qcost = qtmp * cos(dtheta*(jj+1));
	for (int kk=0; kk<G_IN.num_phi; kk++){
	  idx_dir = jj*G_IN.num_phi + kk + 1;
	  //if (ii==0) printf("%d\n",idx_dir);
	  qq[idx3(ii, idx_dir, 0, nq, nq_dir)] = qsint * cos(dphi*kk);
	  qq[idx3(ii, idx_dir, 1, nq, nq_dir)] = qsint * sin(dphi*kk);
	  qq[idx3(ii, idx_dir, 2, nq, nq_dir)] = qcost;
	  //if (ii==0) printf("%f %f %f\n",qq[idx3(ii, idx_dir, 0, nq, nq_dir)],qq[idx3(ii, idx_dir, 1, nq, nq_dir)],qq[idx3(ii, idx_dir, 2, nq, nq_dir)]);
	}
      }
      // theta = pi
      qq[idx3(ii, nq_dir-1, 0, nq, nq_dir)] = 0;
      qq[idx3(ii, nq_dir-1, 1, nq, nq_dir)] = 0;
      qq[idx3(ii, nq_dir-1, 2, nq, nq_dir)] = -qtmp;
      //if (ii==0) printf("%f %f %f\n",qq[idx3(ii, nq_dir-1, 0, nq, nq_dir)],qq[idx3(ii, nq_dir-1, 1, nq, nq_dir)],qq[idx3(ii, nq_dir-1, 2, nq, nq_dir)]);
    }
  }

  // Initialize intermediate scattering function
  for (int ii=0; ii<nq; ii++){
      for (int jj=0; jj<n_files; jj++){
  	cf[idx2(ii,jj,nq)] = 0.0 + I*0.0;
      }
  }

  // Output
  *out_n_atoms = n_atoms;
  *out_LL = LL;
  *out_dq = dq;
  *out_nq = nq;
  *out_nq_dir = nq_dir;
  *out_sim_box = sim_box;
  *out_xx = xx;
  *out_yy = yy;
  *out_zz = zz;
  *out_vx = vx;
  *out_vy = vy;
  *out_vz = vz;
  *out_qq = qq;
  *out_dfk = dfk;
  *out_dfmk = dfmk;
  *out_cf = cf;

}

// ------ Free arrays associated with correlation function calculation ------
void free_corr(double *sim_box, double *xx, double *yy, double *zz,
	       double *vx, double *vy, double *vz,
	       double *qq, double complex *dfk, double complex *dfmk,
	       double complex *cf){

  free(sim_box);
  free(xx);
  free(yy);
  free(zz);
  free(vx);
  free(vy);
  free(vz);
  free(qq);
  free(dfk);
  free(dfmk);
  free(cf);

}


// -------------------------------------------------------------------
// FUNCTIONS USED TO ACCESS THE MULTIDIMENSIONAL ARRAYS
// -------------------------------------------------------------------

// ------ Access element of two dimensional array ------
int idx2(int xx, int yy, int x_size) { 
  return (yy * x_size) + xx; 
}

// ------ Access element of three dimensional array ------
int idx3(int xx, int yy, int zz, 
	 int x_size, int y_size) { 
  return (zz * x_size * y_size) + (yy * x_size) + xx; 
}


// -------------------------------------------------------------------
// FUNCTIONS USED TO READ LAMMPS OUTPUT
// -------------------------------------------------------------------

// ------ Obtain file names with the pattern specified in input ------
void get_file_names(char *file_pattern){

  int glob_check;
  glob_check = glob(file_pattern, GLOB_MARK, NULL, &G_FILE_NAMES);
  if (glob_check == GLOB_NOSPACE){
    fprintf(stderr, "ERROR: Glob, out of memory\n");
    exit(EXIT_FAILURE);
  }
  else if (glob_check == GLOB_ABORTED){
    fprintf(stderr, "ERROR: Glob, read error\n");
    exit(EXIT_FAILURE);
  }
  else if (glob_check == GLOB_NOMATCH){
    fprintf(stderr, "ERROR: Glob, no match\n");
    exit(EXIT_FAILURE);
  }

}

// ------ Read one file for initialization purposes ------
void read_file_init(int *out_n_atoms, double *out_min_box_size){

  // Variables
  int n_atoms = -1, foi; 
  double *sim_box = NULL; 
  bool n_atoms_flag = false, box_flag = false,
       input_read = false, fob;
  gzFile fid;
  char line[LINE_BUF_SIZE];

  // Open file
  fid = gzopen(G_FILE_NAMES.gl_pathv[0], "r");
  read_line(fid, line);

  // Read file
  while (!gzeof(fid) && !input_read)
    {      

      // Get what type of input to expect
      get_headers_info(line, &fob, &n_atoms_flag,
		       &box_flag, &fob,
		       &foi, &foi, &foi,
		       &foi, &foi, &foi);

      // Read number of atoms
      if (n_atoms_flag){
	read_n_atoms(fid, line, &n_atoms);	
  	n_atoms_flag = false;
      }
      
      // Read simulation box info
      if (box_flag){
	sim_box = malloc(sizeof(double) * 3);
	if (sim_box == NULL){
	  fprintf(stderr, "ERROR: Failed simulation box allocation\n");
	  exit(EXIT_FAILURE);
	}
	read_sim_box(fid, line, sim_box);	
	box_flag = false;
      }

      // Exit if all the necessary input was read
      if (n_atoms != -1 && sim_box != NULL)
	input_read = true;
      
      // Read new line
      read_line(fid, line);

    }

  // Close file
  gzclose(fid);
  
  // Output
  *out_n_atoms = n_atoms;
  *out_min_box_size = sim_box[0];
  if (sim_box[1] < *out_min_box_size) *out_min_box_size = sim_box[1];
  if (sim_box[2] < *out_min_box_size) *out_min_box_size = sim_box[2];

  // Free memory
  free(sim_box);

}

// ------ Read one file and store the coordinates ------
void read_file(int file_id, int *time_step, int *n_atoms, 
	       double *sim_box, double *xx, double *yy, double *zz,
	       double *vx, double *vy, double *vz){
  
  // Variables
  bool time_flag = false;
  bool n_atoms_flag = false;
  bool box_flag = false;
  bool config_flag = false;
  int x_idx = -1, y_idx = -1, z_idx = -1;
  int vx_idx = -1, vy_idx, vz_idx = -1;
  gzFile fid;
  char line[LINE_BUF_SIZE];

  // Open file
  fid = gzopen(G_FILE_NAMES.gl_pathv[file_id], "r");
  read_line(fid, line);

  // Read line by line until enf of file
  while (!gzeof(fid))
    {      

      // Get what type of input to expect
      get_headers_info(line, &time_flag, &n_atoms_flag,
		       &box_flag, &config_flag,
		       &x_idx, &y_idx, &z_idx,
		       &vx_idx, &vy_idx, &vz_idx);

      // Read time-step
      if (time_flag){
	read_time_step(fid, line, time_step);	
  	time_flag = false;
      }

      // Read number of atoms
      if (n_atoms_flag){
	read_n_atoms(fid, line, n_atoms);	
  	n_atoms_flag = false;
      }
      
      // Read simulation box info
      if (box_flag){
	read_sim_box(fid, line, sim_box);	
	box_flag = false;
      }


      /* // Read configuration */
      if (config_flag){
      	read_config(fid, line, *n_atoms,
      		    x_idx, y_idx, z_idx,
      		    vx_idx, vy_idx, vz_idx,
      		    xx, yy, zz,
      		    vx, vy, vz);
      	config_flag = false;
      }
      
      // Read new line
      read_line(fid, line);

    }

  // Close file
  gzclose(fid);

}

// ------ Read information in the headers ------
void get_headers_info(char *line, bool *time_flag, bool *n_atoms_flag,
		      bool *box_flag, bool *config_flag,
		      int *x_idx, int *y_idx, int *z_idx,
		      int *vx_idx, int *vy_idx, int *vz_idx){

  // The headers contain information about what type on input
  // can be expected in the next line

  // Check what type of input to expect
  if (strncmp(line,"ITEM: TIMESTEP", 14) == 0){
    *time_flag = true; // Next line contains the time step
  }
  else if (strncmp(line,"ITEM: NUMBER OF ATOMS", 21) == 0){
    *n_atoms_flag = true; // Next line contains the number of atoms
  }
  else if (strncmp(line,"ITEM: BOX BOUNDS", 16) == 0){
    *box_flag = true; // Next lines contain the simulation box
  }
  else if (strncmp(line,"ITEM: ATOMS", 11) == 0){
    *config_flag = true; // Next lines contains  the configuration
    // Check in which column are the positions and the velocities
    char *key = strtok(line, " ");
    int key_counter = -2;
    while (key != NULL){
      if (strcmp(key, "x") == 0) *x_idx = key_counter;
      else if (strcmp(key, "y") == 0) *y_idx = key_counter;
      else if (strcmp(key, "z") == 0) *z_idx = key_counter;
      else if (strcmp(key, "vx") == 0) *vx_idx = key_counter;
      else if (strcmp(key, "vy") == 0) *vy_idx = key_counter;
      else if (strcmp(key, "vz") == 0) *vz_idx = key_counter;
      key = strtok(NULL, " ");
      key_counter++;
    }
    // Exit if the configuration is missing necessary information
    if (*x_idx == -1 || *y_idx == -1 || *z_idx == -1){
      fprintf(stderr,"ERROR: the configuration file does not contain"
	      " full information on the particle's positions\n");
      exit(EXIT_FAILURE);
    }
    if (*vx_idx == -1 || *vy_idx == -1 || *vz_idx == -1){
      fprintf(stderr,"ERROR: the configuration file does not contain"
	      " full information on the particle's velocity\n");
      exit(EXIT_FAILURE);
    }
  }
  
      // Exit if the file does not have the correct format
  if ( (*time_flag == false) && (*n_atoms_flag == false) &&
       (*box_flag == false) && (*config_flag) == false){
    fprintf(stderr,"ERROR: the configuration file has an unexpected format\n");
    exit(EXIT_FAILURE);
  }
  
}

// ------ Read time step ------
void read_time_step(gzFile fid, char *line, int *tstep){

  read_line(fid, line);
  *tstep = read_int(line);

}

// ------ Read number of atoms ------
void read_n_atoms(gzFile fid, char *line, int *nat){

  read_line(fid, line);
  *nat = read_int(line);

}

// ------ Read box information ------
void read_sim_box(gzFile fid, char *line, double *sim_box){

  double lmin = -1, lmax = -1;

  for (int ii=0; ii<3; ii++){
    
    // Read line
    read_line(fid,line);

    // Extract simulation box size from line
    char *tok = strtok(line, " ");
    int tok_counter = 0;
    while (tok != NULL){
      if (tok_counter == 0) lmin = read_double(tok);
      else if (tok_counter == 1) lmax = read_double(tok);
      tok = strtok(NULL, " ");
      tok_counter++;
    }

    // Check that the simulation box was read correctly
    if (lmin == -1 || lmax == -1){
      fprintf(stderr,"ERROR: the simulation box information could"
	      " not be read from the configuration file\n");
      exit(EXIT_FAILURE);
    }

    // Set simulation box size
    sim_box[ii] = lmax - lmin;

  }

}

// ------ Read configuration ------
void read_config(gzFile fid, char *line, int n_atoms,
		 int x_idx, int y_idx, int z_idx,
		 int vx_idx, int vy_idx, int vz_idx,
		 double *xx, double *yy, double *zz,
		 double *vx, double *vy, double *vz){

  for (int ii=0; ii<n_atoms; ii++){
    
    // Initialize the positions and velocities 
    xx[ii] = -1;
    yy[ii] = -1;
    zz[ii] = -1;
    vx[ii] = -1;
    vy[ii] = -1;
    vz[ii] = -1;

    // Read line 
    read_line(fid,line);

    // Extract coordinates and velocities from line
    char *tok = strtok(line, " ");
    int tok_counter = 0;
    while (tok != NULL){
      if (tok_counter == x_idx) xx[ii] = read_double(tok);
      if (tok_counter == y_idx) yy[ii] = read_double(tok);
      if (tok_counter == z_idx) zz[ii] = read_double(tok);
      if (tok_counter == vx_idx) vx[ii] = read_double(tok);
      if (tok_counter == vy_idx) vy[ii] = read_double(tok);
      if (tok_counter == vz_idx) vz[ii] = read_double(tok);
      tok = strtok(NULL, " ");
      tok_counter++;
    }

  }

}


// ------ Read one line from file ------
void read_line(gzFile fid, char *line){

  // Read line
  gzgets( fid, line, LINE_BUF_SIZE );
  
  // Check if the eof has been reached
  if (!gzeof(fid)){

    // Check for buffer overflow
    if(line[strlen(line)-1] != '\n'){
      fprintf(stderr,"ERROR: buffer overflow while reading input"
	      " (increase  LINE_BUF_SIZE if it is necessary"
	      " to read more than %d characters per line)\n", 
	      LINE_BUF_SIZE);
      exit(EXIT_FAILURE);
    }

    // Drop the new line character
    line[strlen(line)-1] = '\0';
  

  }

}


// ------ Read one integer from string ------
int read_int(char *str){

  long num;
  char *end;

  errno = 0;

  // Read number from string
  num = strtol(str, &end, 10);        //10 specifies base-10

  // Check if a number was read
  if (end == str){    
    fprintf(stderr, "ERROR: can't convert string to number\n");
    exit(EXIT_FAILURE);
  }

  // If sizeof(int) == sizeof(long), we have to explicitly check for overflows
  if ((num == LONG_MAX || num == LONG_MIN) && errno == ERANGE){  
    fprintf(stderr, "ERROR: number out of range for LONG\n");
    exit(EXIT_FAILURE);
  }
  // Because strtol produces a long, check for overflow
  if ( (num > INT_MAX) || (num < INT_MIN) ){
    fprintf(stderr, "ERROR: number out of range for INT\n");
    exit(EXIT_FAILURE);
  }

  // Return
  return (int)num;

}

// ------ Read one double from string ------
double read_double(char *str){

  double num;
  char *end;

  errno = 0;

  // Read number from string
  num = strtod(str, &end);        

  // Check if a number was read
  if (end == str){     
    fprintf(stderr, "ERROR: can't convert string to number\n");
    exit(EXIT_FAILURE);
  }

  // Check for overflows
  if ((num == DBL_MAX || num == DBL_MIN) && errno == ERANGE){  
    fprintf(stderr, "ERROR: number out of range for LONG\n");
    exit(EXIT_FAILURE);
  }

  // Return
  return num;

}

// ------ Check that the input from one file is consistent with the initialization
void check_consistency(int n_atoms_ref, int n_atoms_check, double LL_ref, double *sim_box){

  // Check if the number of atoms is consistent
  if (n_atoms_check != n_atoms_ref){
    fprintf(stderr, "ERROR: Number of atoms not consistent with initizialization\n");
    exit(EXIT_FAILURE);
  }

  // Check if the simulation box is consistent
  double LL_check = sim_box[0];
  if (sim_box[1] < LL_check) LL_check = sim_box[1];
  if (sim_box[2] < LL_check) LL_check = sim_box[2];
  if ( abs(LL_ref-LL_check) > 1e-10 ) {
      fprintf(stderr,"ERROR: Simulation box not consistent with initialization\n");
      exit(EXIT_FAILURE);
  }  
  

}

