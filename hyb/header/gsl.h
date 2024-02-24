#ifndef GSL_H
#define GSL_H
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>

typedef struct
{
  int iter;
  double step;
  double max_step;
  double tol;
  gsl_vector *x1;
  gsl_vector *dx1;
  gsl_vector *x2;
  double pnorm;
  gsl_vector *p; 
  double g0norm;
  gsl_vector *g0;
}
conjugate_state_t; //fr, pr common

typedef struct
{
  int iter;
  double step;
  double max_step;
  double tol;
  gsl_vector *x1;
  gsl_vector *dx1;
  gsl_vector *x2;
  double g0norm;
  double pnorm;
  gsl_vector *p; 
  gsl_vector *x0;
  gsl_vector *g0;
  gsl_vector *dx0;
  gsl_vector *dg0;
}
vector_bfgs_state_t;

typedef struct
{
	gsl_function_fdf fdf_linear;
	gsl_multimin_function_fdf *fdf;
	/* fixed values */
	const gsl_vector *x;
	const gsl_vector *g;
	const gsl_vector *p;

	/* cached values, for x(alpha) = x + alpha * p */
	double f_alpha;
	double df_alpha;
	gsl_vector *x_alpha;
	gsl_vector *g_alpha;

	/* cache "keys" */
	double f_cache_key;
	double df_cache_key;
	double x_cache_key;
	double g_cache_key;
}
wrapper_t;

typedef struct
{
  int iter;
  double step;
  double g0norm;
  double pnorm;
  double delta_f;
  double fp0;                   /* f'(0) for f(x-alpha*p) */
  gsl_vector *x0;
  gsl_vector *g0;
  gsl_vector *p; 
  /* work space */
  gsl_vector *dx0;
  gsl_vector *dg0;
  gsl_vector *x_alpha;
  gsl_vector *g_alpha;
  /* wrapper function */
  wrapper_t wrap;
  /* minimization parameters */
  double rho;
  double sigma;
  double tau1;
  double tau2;
  double tau3;
  int order;
}
vector_bfgs2_state_t;

#endif
