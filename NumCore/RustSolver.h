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
#pragma once
#include <FECore/LinearSolver.h>
#include <FECore/CompactUnSymmMatrix.h>
#include <FECore/CompactSymmMatrix.h>

//! Parallel sparse direct solver implemented in Rust (rust/febio_rs_solver).
//! Symmetric matrices: supernodal LDL^T with METIS nested-dissection ordering,
//! elimination-tree + dense-kernel parallelism and PARDISO-style pivot perturbation.
//! Unsymmetric matrices: supernodal sparse LU with partial pivoting.
//! Select with <linear_solver type="rust"/>.
class RustSolver : public LinearSolver
{
public:
	RustSolver(FEModel* fem);
	~RustSolver();

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override;
	bool SetSparseMatrix(SparseMatrix* pA) override;

	bool PreProcess() override;
	bool Factor() override;
	bool BackSolve(double* x, double* y) override;
	void Destroy() override;

	void SetPrintLevel(int n) override { m_printLevel = n; }

private:
	CompactMatrix*	m_pA;
	bool			m_symmetric;
	void*			m_solver;	// FrsSolver*

	int		m_nthreads;		// 0 = use omp_get_max_threads()
	int		m_ordering;		// 0 = metis, 1 = amd
	int		m_refine;		// 0 = auto, n>0 = always n steps
	double	m_pivotEps;
	int		m_printLevel;
	bool	m_reuseOrdering;	// reuse ordering for slightly changed patterns (contact)
	int		m_ompBlocktime;		// OpenMP (libiomp5) blocktime in ms while using this solver; -1 = unchanged

	DECLARE_FECORE_CLASS();
};
