#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <glob.h>
#include <zlib.h>
#include <math.h>
#include <complex.h>
#include "lmp_corr.h"

# define LINE_BUF_SIZE 1024 // maximum line length from input

static glob_t G_FILE_NAMES; // List of files to load
static double G_QMAX; // wave-vector cutoff


// ------ Main ------
int main(int argc, char *argv[]) {

  // Default values
  // Input file
  char *file_pattern;
  file_pattern = (char*) malloc(sizeof(char) * 21);
  strcpy(file_pattern, "trajectories*.dat");
  // Wave-vector cutoff
  G_QMAX = 10.0;

  // Parse command line
  static const char *usage = "Description coming soon...";
  int opt;
  while ((opt = getopt(argc, argv, "hi:q:")) != -1) {
    switch (opt) {
    case 'i':
      file_pattern = (char*)realloc(file_pattern, sizeof(char) * (strlen(optarg)+1));
      strcpy(file_pattern, optarg);
      break;
    case 'q':
      G_QMAX = read_double(optarg);
      break;
    case 'h':
      printf(usage, argv[0]);
      exit(EXIT_SUCCESS);
    case '?':
      printf("%s: unknown option: -%c\n", argv[0], optopt);
      printf(usage, argv[0]);
      exit(EXIT_FAILURE);
      break;
    case ':':
      printf("%s: missing option argument for option: -%c\n", argv[0], optopt);
      printf(usage, argv[0]);
      exit(EXIT_FAILURE);
      break;
    default:
      printf(usage, argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  
  
  // Get file names
  get_file_names(file_pattern);

  // Compute intermediate scattering function
  clock_t start = clock();
  isf();
  clock_t end = clock();
  printf("Elapsed time: %f seconds\n",
         (double)(end - start) / CLOCKS_PER_SEC);

  
  // Free memory
  free(file_pattern);
  globfree(&G_FILE_NAMES);

  return 0;

}


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

// ------ Compute intermediate scattering function ------
void isf(){

  // Variables to store the information read from the file
  int time_step = 0; // Timestep
  int n_atoms = 0; // Number of atoms
  double *sim_box = NULL; // Simulation box
  double *xx = NULL, *yy = NULL, *zz = NULL; // Coordinates
  double *vx = NULL, *vy = NULL, *vz = NULL; // Velocities

  // Variables to compute the density fluctuations
  bool init = true;
  double LL, LL_tmp, dq, *qq = NULL;
  int nq = 0;
  double complex (*drhok)[G_FILE_NAMES.gl_pathc] = NULL;
  double complex (*drhomk)[G_FILE_NAMES.gl_pathc] = NULL;
  
  
  // Loop through the configuration files
  for (int ii=0; ii<G_FILE_NAMES.gl_pathc; ii++){
    
    // Read file content (it is assumed that one file contains one configuration)
    printf("%s\n",G_FILE_NAMES.gl_pathv[ii]);
    fflush(stdout); 
    read_file(ii, &time_step, &n_atoms, &sim_box,
	      &xx, &yy, &zz, &vx, &vy, &vz);

    // Initialize calculations (if necessary)
    if (init) {
      // Wave-vector resolution
      LL = sim_box[0];
      if (sim_box[1] < LL) LL = sim_box[1];
      if (sim_box[2] < LL) LL = sim_box[2];
      dq = 2.0*M_PI/LL;
      // Allocate wave-vector grid
      nq = (int)G_QMAX/dq;
      qq = (double*)malloc(sizeof(double) * nq);
      if (qq == NULL){
	printf("ERROR: Failed allocation for wave-vector grid\n");
	exit(EXIT_FAILURE);
      }
      // Initialize wave-vector grid
      for (int jj=0; jj<nq; jj++){
      	qq[jj] = (jj+1)*dq;
      }
      // Allocate matrix to store density fluctuations
      drhok = malloc(sizeof(*drhok) * nq);
      drhomk = malloc(sizeof(*drhomk) * nq);
      if (drhok == NULL || drhomk == NULL){
	printf("ERROR: Failed allocation for density fluctuations\n");
	exit(EXIT_FAILURE);
      }      
      // De-activate initialization
      init = false;
    }
    else {
      // check that was read is consistent with the previous files
      LL_tmp = sim_box[0];
      if (sim_box[1] < LL_tmp) LL_tmp = sim_box[1];
      if (sim_box[2] < LL_tmp) LL_tmp = sim_box[2];
      if ( abs(LL_tmp-LL) > 1e-10 ) {
	printf("ERROR: The simulation box changed for file %s\n", 
	       (G_FILE_NAMES.gl_pathv[ii]));
	exit(EXIT_FAILURE);
      }    
    }
 

    // Compute density fluctuations
    for (int jj=0; jj<nq; jj++){
      drhok[jj][ii] = 0.0 + I*0.0;
      drhomk[jj][ii] = 0.0 + I*0.0;
      for (int kk=0; kk<n_atoms; kk++){
	drhok[jj][ii] += cexp(I * qq[jj] * xx[kk]);
	drhomk[jj][ii] += cexp(-I * qq[jj] * xx[kk]);
      }
    }
        
    // Free memory associated to file
    free(sim_box);
    free(xx);
    free(yy);
    free(zz);
    free(vx);
    free(vy);
    free(vz);
  
  }

  // Compute intermediate scattering function   
  int lag = (int)G_FILE_NAMES.gl_pathc;
  double complex (*Fkt)[lag] = malloc(sizeof(*Fkt) * nq);
  if (Fkt == NULL){
    printf("ERROR: Failed allocation for intermediate scattering function \n");
    exit(EXIT_FAILURE);
  }
  int norm_fact;
  for (int jj=0; jj<nq; jj++){
    for (int kk=0; kk<lag; kk++){
      Fkt[jj][kk] = 0.0;
      norm_fact = 0;
      for (int ll=0; ll<lag-kk; ll++){
      	Fkt[jj][kk] += drhomk[jj][ll]*drhok[jj][kk+ll];
	norm_fact++;
      }
      Fkt[jj][kk] /= norm_fact*n_atoms;
    }
  }
  
  // Write output
  FILE *fid;
  fid = fopen("isf_real.dat", "w");
  /* fprintf(fid, "#########"); */
  /* for (int ii=0; ii<lag; ii++) fprintf(fid,"%d        ",ii); */
  /* fprintf(fid, "\n"); */
  for (int jj=0; jj<nq; jj++){
    fprintf(fid, "%.8f ", qq[jj]);
    for (int kk=0; kk<lag; kk++){
      fprintf(fid, "%.8f ", creal(Fkt[jj][kk]));
    }
    fprintf(fid, "\n");
  }
  fclose(fid);

  fid = fopen("isf_imag.dat", "w");
  /* fprintf(fid, "#########"); */
  /* for (int ii=0; ii<lag; ii++) fprintf(fid,"%d        ",ii); */
  /* fprintf(fid, "\n"); */
  for (int jj=0; jj<nq; jj++){
    fprintf(fid, "%.8f ", qq[jj]);
    for (int kk=0; kk<lag; kk++){
      fprintf(fid, "%.8f ", cimag(Fkt[jj][kk]));
    }
    fprintf(fid, "\n");
  }
  fclose(fid);  

  
  // Free memory associated to density fluctuations
  free(qq);
  free(drhok);
  free(drhomk);

  // Free memory associated to intermediate scattering function
  free(Fkt);

}

// ------ Read one file ------
void read_file(int file_id, int *out_time_step, int *out_n_atoms, 
	       double **out_sim_box,
	       double **out_xx, double **out_yy, double **out_zz,
	       double **out_vx, double **out_vy, double **out_vz){
  
  // Variables
  int time_step = -1; 
  int n_atoms = -1; 
  double *sim_box = NULL; 
  double *xx = NULL, *yy = NULL, *zz = NULL;
  double *vx = NULL, *vy = NULL, *vz = NULL;
  bool time_flag = false;
  bool n_atoms_flag = false;
  bool box_flag = false;
  bool config_flag = false;
  int x_idx = -1, y_idx = -1, z_idx = -1;
  int vx_idx = -1, vy_idx = -1, vz_idx = -1;
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
	read_time_step(fid, line, &time_step);	
  	time_flag = false;
      }

      // Read number of atoms
      if (n_atoms_flag){
	read_n_atoms(fid, line, &n_atoms);	
  	n_atoms_flag = false;
      }
      
      // Read simulation box info
      if (box_flag){
	sim_box = (double*)malloc(sizeof(double) * 3);
	if (sim_box == NULL){
	  fprintf(stderr, "ERROR: Failed simulation box allocation\n");
	  exit(EXIT_FAILURE);
	}
	read_sim_box(fid, line, sim_box);	
	box_flag = false;
      }


      // Read configuration
      if (config_flag){
	xx = (double*)malloc(sizeof(double) * n_atoms);
	yy = (double*)malloc(sizeof(double) * n_atoms);
	zz = (double*)malloc(sizeof(double) * n_atoms);
	vx = (double*)malloc(sizeof(double) * n_atoms);
	vy = (double*)malloc(sizeof(double) * n_atoms);
	vz = (double*)malloc(sizeof(double) * n_atoms);
	if ((xx == NULL) || (xx == NULL) || (xx == NULL) ||
	    (vx == NULL) || (vy == NULL) || (vz == NULL) ){
          fprintf(stderr, "ERROR: Failed configuration allocation\n");
          exit(EXIT_FAILURE);
        }
	read_config(fid, line, n_atoms,
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

  // Check that all the necessary information was read from the file
  if ( time_step == -1 || n_atoms == -1 || sim_box == NULL ||
       xx == NULL || yy == NULL || zz == NULL || 
       vx == NULL || vy == NULL || vz == NULL ||
       x_idx == -1 || y_idx == -1 || z_idx == -1 ||
       vx_idx == -1 || vy_idx == -1 || vz_idx == -1){
    fprintf(stderr, "ERROR: It was not possible to read all necessary information"
	    " from the configuration file\n");
    exit(EXIT_FAILURE);
  }

  // Output
  *out_time_step = time_step;
  *out_n_atoms = n_atoms;
  *out_sim_box = sim_box;
  *out_xx = xx;
  *out_yy = yy;
  *out_zz = zz;
  *out_vx = vx;
  *out_vy = vy;
  *out_vz = vz;

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
  if (end == str)    
    fprintf(stderr, "ERROR: can't convert string to number\n");

  // If sizeof(int) == sizeof(long), we have to explicitly check for overflows
  if ((num == LONG_MAX || num == LONG_MIN) && errno == ERANGE)  
    fprintf(stderr, "ERROR: number out of range for LONG\n");

  // Because strtol produces a long, check for overflow
  if ( (num > INT_MAX) || (num < INT_MIN) )
    fprintf(stderr, "ERROR: number out of range for INT\n");

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
  if (end == str)     
    fprintf(stderr, "ERROR: can't convert string to number\n");

  // Check for overflows
  if ((num == DBL_MAX || num == DBL_MIN) && errno == ERANGE)  
    fprintf(stderr, "ERROR: number out of range for LONG\n");

  // Return
  return num;

}
