/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#include "stdafx.h"
#include "RustSolver.h"
#include <FECore/log.h>
#include <FECore/fecore_threads.h>
#include <stdio.h>

#ifdef RUST_SOLVER
#include <febio_rs_solver.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

#if defined(RUST_SOLVER) && defined(__GNUC__)
// Intel OpenMP runtime (libiomp5, used with MKL): after a parallel region, idle OpenMP
// threads spin for KMP_BLOCKTIME (default 200 ms) before sleeping. Those spinning threads
// compete for the cores the Rust solver's own thread pool needs. Weak reference so that
// builds using GNU libgomp (no such function) still link.
extern "C" void kmp_set_blocktime(int) __attribute__((weak));
#endif

BEGIN_FECORE_CLASS(RustSolver, LinearSolver)
	ADD_PARAMETER(m_nthreads, "nthreads");
	ADD_PARAMETER(m_ordering, "ordering");
	ADD_PARAMETER(m_refine  , "refine");
	ADD_PARAMETER(m_pivotEps, "pivot_eps");
	ADD_PARAMETER(m_printLevel, "print_level");
	ADD_PARAMETER(m_reuseOrdering, "reuse_ordering");
	ADD_PARAMETER(m_ompBlocktime, "omp_blocktime");
END_FECORE_CLASS();

RustSolver::RustSolver(FEModel* fem) : LinearSolver(fem), m_pA(nullptr), m_symmetric(true), m_solver(nullptr)
{
	m_nthreads = 0;
	m_ordering = 0;
	m_refine = 0;
	m_pivotEps = 0.0;
	m_printLevel = 0;
	m_reuseOrdering = true;
	m_ompBlocktime = 5;
}

SparseMatrix* RustSolver::CreateSparseMatrix(Matrix_Type ntype)
{
	switch (ntype)
	{
	case REAL_SYMMETRIC     : m_symmetric = true ; m_pA = new CompactSymmMatrix(0); break;
	case REAL_UNSYMMETRIC   :
	case REAL_SYMM_STRUCTURE: m_symmetric = false; m_pA = new CRSSparseMatrix(0); break;
	default:
		assert(false);
		m_pA = nullptr;
	}
	return m_pA;
}

bool RustSolver::SetSparseMatrix(SparseMatrix* pA)
{
	m_pA = dynamic_cast<CompactMatrix*>(pA);
	if (m_pA == nullptr) return false;
	m_symmetric = (dynamic_cast<CompactSymmMatrix*>(pA) != nullptr);
	if (!m_symmetric && (dynamic_cast<CRSSparseMatrix*>(pA) == nullptr)) return false;
	return true;
}

bool RustSolver::PreProcess()
{
#ifdef RUST_SOLVER
	if (m_pA == nullptr) return false;

	// The Rust solver object survives Destroy() (FEBio calls Destroy + PreProcess whenever it
	// rebuilds the matrix profile), so it can detect identical / similar sparsity patterns and
	// reuse its symbolic analysis or ordering.
	// use the number of threads requested by the user (the Rust solver is bitwise
	// deterministic for any thread count)
	int nthreads = m_nthreads;
	if (nthreads <= 0) nthreads = FECoreThreads::RequestedThreads();
	if (m_solver == nullptr)
	{
#if defined(__GNUC__)
		// let idle OpenMP threads sleep soon (5 ms) instead of spinning while the Rust
		// thread pool runs; a few ms keep them awake between the back-to-back element
		// loops of one assembly (-1 = leave the OpenMP runtime setting unchanged)
		if ((m_ompBlocktime >= 0) && kmp_set_blocktime) kmp_set_blocktime(m_ompBlocktime);
#endif
		FrsSolver* s = frs_create(nthreads);
		m_solver = s;
		frs_set_option(s, FRS_OPT_ORDERING, m_ordering);
		if (m_pivotEps > 0) frs_set_option(s, FRS_OPT_PIVOT_EPS, m_pivotEps);	// 0 = solver default (1e-13 of max|Aij| for LDLT and LU)
		frs_set_option(s, FRS_OPT_REFINE, m_refine);
		frs_set_option(s, FRS_OPT_REUSE_ORDERING, m_reuseOrdering ? 1.0 : 0.0);
	}
	FrsSolver* s = (FrsSolver*)m_solver;

	int n = m_pA->Rows();
	if (n == 0) return LinearSolver::PreProcess();
	int ret = (m_symmetric ?
		frs_analyze_symmetric  (s, n, m_pA->Pointers(), m_pA->Indices(), m_pA->Offset()) :
		frs_analyze_unsymmetric(s, n, m_pA->Pointers(), m_pA->Indices(), m_pA->Offset()));
	if (ret != 0) { feLogError("Rust solver: symbolic analysis failed."); return false; }

	if (m_printLevel > 0)
	{
		FrsStats st; frs_get_stats(s, &st);
		feLog("\tRust solver: n = %lld, nnz(A) = %lld, nnz(L) = %lld, supernodes = %lld, threads = %d, analyses = %lld, reused = %lld\n",
			(long long)st.n, (long long)st.nnz_a, (long long)st.nnz_l, (long long)st.n_supernodes, nthreads,
			(long long)st.n_analyze, (long long)st.n_reused);
	}
	return LinearSolver::PreProcess();
#else
	feLogError("FEBio was not built with the Rust solver (USE_RUST_SOLVER).");
	return false;
#endif
}

bool RustSolver::Factor()
{
#ifdef RUST_SOLVER
	if (m_pA->Rows() == 0) return true;
	int ret = frs_factor((FrsSolver*)m_solver, m_pA->Values(), (int64_t)m_pA->NonZeroes());
	if (ret != 0) { feLogError("Rust solver: numerical factorization failed."); return false; }
	if (m_printLevel > 1)
	{
		FrsStats st; frs_get_stats((FrsSolver*)m_solver, &st);
		if (st.n_perturbed > 0) feLog("\tRust solver: %lld perturbed pivots\n", (long long)st.n_perturbed);
	}
	return true;
#else
	return false;
#endif
}

bool RustSolver::BackSolve(double* x, double* b)
{
#ifdef RUST_SOLVER
	if (m_pA->Rows() == 0) return true;
	int ret = frs_solve((FrsSolver*)m_solver, x, b);
	if (ret != 0) { feLogError("Rust solver: solve failed."); return false; }
	UpdateStats(1);
	return true;
#else
	return false;
#endif
}

void RustSolver::Destroy()
{
	// free the factors only; the symbolic analysis is kept for pattern reuse
#ifdef RUST_SOLVER
	if (m_solver) frs_release_numeric((FrsSolver*)m_solver);
#endif
	LinearSolver::Destroy();
}

RustSolver::~RustSolver()
{
#ifdef RUST_SOLVER
	if (m_solver)
	{
		if (m_printLevel > 0)
		{
			FrsStats st; frs_get_stats((FrsSolver*)m_solver, &st);
			feLog("\tRust solver: %lld analyses (%lld patterns reused) %.3fs, %lld factorizations %.3fs, %lld solves %.3fs, %lld fallbacks\n",
				(long long)st.n_analyze, (long long)st.n_reused, st.t_analyze, (long long)st.n_factor, st.t_factor,
				(long long)st.n_solve, st.t_solve, (long long)st.n_fallback);
		}
		frs_destroy((FrsSolver*)m_solver);
		m_solver = nullptr;
	}
#endif
}
