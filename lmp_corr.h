#ifndef LMP_CORR_H
#define LMP_CORR_H

#include <zlib.h>
#include <complex.h>

typedef struct {

  bool isf;
  bool lvcf;
  char *config_file;
  double q_max;
  double num_theta;
  double num_phi;
  int num_threads;

} input;

// -------------------------------------------------------------------
// CALLER FOR THE OTHER FUNCTIONS IN THE FILE
// -------------------------------------------------------------------

void analyze_lmp(input in);

// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE INTERMEDIATE SCATTERING FUNCTION
// -------------------------------------------------------------------

void isf();


void isf_init(int *out_n_atoms, double *out_LL, double *out_dq, int *out_nq,
              int *out_nq_dir, double **out_sim_box,
              double **out_xx, double **out_yy, double **out_zz,
              double **out_vx, double **out_vy, double **out_vz,
              double **out_qq, double complex **out_drhok,
              double complex **out_drhomk, double complex **out_fkt);

void isf_output(double complex *fkt, double dq, int nq, int n_files);

void isf_free(double *out_sim_box, double *xx, double *yy, double *zz,
	      double *vx, double *vy, double *vz,
	      double *qq, double complex *drhok, double complex *drhomk,
	      double complex *fkt);


// -------------------------------------------------------------------
// FUNCTIONS USED TO COMPUTE THE LONGITUDINAL VELOCITY CORRELATIONS
// -------------------------------------------------------------------

void lvcf();

void lvcf_init(int *out_n_atoms, double *out_LL, double *out_dq, int *out_nq,
              int *out_nq_dir, double **out_sim_box,
              double **out_xx, double **out_yy, double **out_zz,
              double **out_vx, double **out_vy, double **out_vz,
              double **out_qq, double complex **out_dvk,
              double complex **out_dvmk, double complex **out_ckt);

void lvcf_output(double complex *ckt, double dq, int nq, int n_files);

void lvcf_free(double *out_sim_box, double *xx, double *yy, double *zz,
	      double *vx, double *vy, double *vz,
	      double *qq, double complex *dvk, double complex *dvmk,
	      double complex *ckt);


// -------------------------------------------------------------------
// FUNCTIONS USED TO ACCESS THE MULTIDIMENSIONAL ARRAYS
// -------------------------------------------------------------------


int idx2(int xx, int yy, int x_size);

int idx3(int xx, int yy, int zz, int x_size, int y_size);


// -------------------------------------------------------------------
// FUNCTIONS USED TO READ LAMMPS OUTPUT
// -------------------------------------------------------------------

void get_file_names(char *file_pattern);

void read_file_init(int *num_atoms, double *min_box_size);

void read_file(int file_id, int *out_time_step, int *out_n_atoms,
	       double *sim_box, double *xx, double *yy, double *zz,
	       double *vx, double *vy, double *vz);

void get_headers_info(char *line, bool *time_flag, bool *n_atoms_flag,
                      bool *box_flag, bool *config_flag,
                      int *x_idx, int *y_idx, int *z_idx,
                      int *vx_idx, int *vy_idx, int *vz_idx);

void read_time_step(gzFile fid, char *line, int *tstep);

void read_n_atoms(gzFile fid, char *line, int *nat);

void read_sim_box(gzFile fid, char *line, double *box_info);

void read_config(gzFile fid, char *line, int n_atoms,
		 int x_idx, int y_idx, int z_idx,
		 int vx_idx, int vy_idx, int vz_idx,
		 double *xx, double *yy, double *zz,
		 double *vx, double *vy, double *vz);

void read_line(gzFile fid, char *line);

int read_int(char *str);

double read_double(char *str);

void check_consistency(int n_atoms_ref, int n_atoms_check, 
		       double LL_ref, double *sim_box);

#endif

