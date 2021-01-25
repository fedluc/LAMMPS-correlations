#ifndef ISF_H
#define ISF_H

void get_file_names(char *file_pattern);

void compute_isf(int file_id);

void read_file(int file_id, int *out_time_step, int *out_n_atoms,
               double **out_sim_box,
               double **out_xx, double **out_yy, double **out_zz,
               double **out_vx, double **out_vy, double **out_vz);

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

#endif

