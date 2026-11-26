/* C interface of the Rust sparse direct solver (rust/febio_rs_solver). */
#pragma once
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef struct FrsSolver FrsSolver;

typedef struct {
	int64_t n, nnz_a, nnz_l, n_supernodes, n_supervars;
	double  flops;
	double  t_analyze, t_factor, t_solve;
	int64_t n_factor, n_solve, n_perturbed, n_refine;
} FrsStats;

enum {
	FRS_OPT_ORDERING  = 1,	/* 0 = METIS nested dissection, 1 = AMD, 2 = natural */
	FRS_OPT_PIVOT_EPS = 2,	/* relative pivot perturbation (default 1e-8, as PARDISO) */
	FRS_OPT_REFINE    = 3,	/* 0 = refine only if pivots were perturbed, n>0 = always n steps */
	FRS_OPT_MAX_REFINE = 4,	/* max refinement steps in auto mode (default 2) */
	FRS_OPT_PAR_FLOPS = 5	/* flop threshold for multithreaded dense kernels */
};

FrsSolver* frs_create(int nthreads);	/* nthreads = 0 -> all cores */
void frs_destroy(FrsSolver* s);
int  frs_set_option(FrsSolver* s, int key, double val);
/* symmetric: lower triangle in CSC == upper triangle in CSR (PARDISO mtype -2 layout) */
int  frs_analyze_symmetric(FrsSolver* s, int n, const int* ptr, const int* idx, int offset);
/* unsymmetric: full CSR */
int  frs_analyze_unsymmetric(FrsSolver* s, int n, const int* ptr, const int* idx, int offset);
int  frs_factor(FrsSolver* s, const double* vals, int64_t nnz);
int  frs_solve(FrsSolver* s, double* x, const double* b);
void frs_get_stats(const FrsSolver* s, FrsStats* out);

#ifdef __cplusplus
}
#endif
