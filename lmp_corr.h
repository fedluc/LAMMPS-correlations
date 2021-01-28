#ifndef LMP_CORR_H
#define LMP_CORR_H

#include <zlib.h>
#include <complex.h>

typedef struct {
  char *config_file;
  double q_max;
  bool isf;
  bool lvcf;

} input;

void analyze_lmp(input in);

void get_file_names(char *file_pattern);

void isf();


void isf_init(int *out_n_atoms, double *out_LL, double *out_dq, int *out_nq,
	      double **out_sim_box, double **out_xx, double **out_yy, double **out_zz,
              double **out_vx, double **out_vy, double **out_vz,
              double **out_qq, double complex **out_drhok, double complex **out_drhomk,
	      double complex **out_fkt);

void isf_output(double complex *fkt, double *qq, int nq, int n_files);

void isf_free(double *out_sim_box, double *xx, double *yy, double *zz,
	      double *vx, double *vy, double *vz,
	      double *qq, double complex *drhok, double complex *drhomk,
	      double complex *fkt);

void lvcf();

int idx2(int xx, int yy, int x_size);

int idx3(int xx, int yy, int zz, int x_size, int y_size);

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

