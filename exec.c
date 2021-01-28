#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <argp.h>
#include <time.h>
#include "lmp_corr.h"

// ----------------------------------------
// ------ Set-up command line parser ------
// ----------------------------------------

// Documentation
static char doc[] =
  "lmp_corr analyses the configurations generated with "
  "by LAMMPS (see dump command) and it computes the Fourier "
  "space representation of the density and velocity "
  "correlation functions";

// Optional arguments
static struct argp_option options[] = {
  {"conf",    'c', "CONFIG_FILE", 0,
   "File(s) with the configuration (wild-cards are accepted)"}, 
  {"dq",  'q', "DQ_MAX", 0,
  "Wave-vector cutoff"},
  {"isf",   'd', 0, 0,
   "Use this flag to compute the intermediate scattering function" },
  {"lvcf",   'v', 0, 0,
   "Use this flag to compute the longitudinal velocity correlation function" },
  { 0 }
};

// Structure to communicate between main and parser
struct arguments
{
  char *config_file;
  double q_max;
  bool isf;
  bool lvcf;

};


// Single option parser
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{

  struct arguments *arguments = state->input;

  switch (key)
    {

    case 'c':
      arguments->config_file = arg;
      break;
    case 'q':
      arguments->q_max = atof(arg);
      break;
    case 'd':
      arguments->isf = true;
      break;
    case 'v':
      arguments->lvcf = true;
      break;

    case ARGP_KEY_ARG:
      if (state->arg_num > 0) // Too many necessary arguments
        argp_usage (state);

      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

static struct argp argp = { options, parse_opt, 0, doc };

// ----------------------------------------
// ---------------- Main ------------------
// ----------------------------------------


int main (int argc, char **argv){

  struct arguments arguments;

  // Default values for optional arguments
  arguments.config_file  = "trajectories*.dat.gz";
  arguments.q_max = 10.0;
  arguments.isf = false;
  arguments.lvcf = false;

  // Parse command line
  argp_parse (&argp, argc, argv, 0, 0, &arguments);

  // Prepare input for LAMMPS analyzer
  input in;
  in.config_file = arguments.config_file;
  in.q_max = arguments.q_max;
  in.isf = arguments.isf;
  in.lvcf = arguments.lvcf;

  // Start timing
  clock_t t_start = clock();

  // Analyze LAMMPS output
  analyze_lmp(in);

  // Stop timing
  clock_t t_end = clock();
  printf("Analysis complete. Elapsed time: %f seconds\n",
         (double)(t_end - t_start) / CLOCKS_PER_SEC);


  return 0;

}

