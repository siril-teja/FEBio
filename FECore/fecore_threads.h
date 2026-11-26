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
#include "fecore_api.h"
#include <vector>
#include <map>
#include <string>

//! Thread configuration and run-to-run determinism.
//!
//! FEBio's element loops assemble the global matrix and vectors from several OpenMP
//! threads with atomic additions, so the order of floating-point additions (and hence the
//! last bits of the result) changes from run to run. In deterministic mode (the default)
//! FEBio's own OpenMP loops run on one thread, which makes every run bitwise reproducible.
//! Linear solvers that are deterministic themselves (the Rust solver, MKL PARDISO in
//! conditional-numerical-reproducibility mode) still use all requested threads.
//!
//! Deterministic mode can be switched off with the -nondeterministic command line option
//! or the environment variable FEBIO_DETERMINISTIC=0.
namespace FECoreThreads {

	//! Records the number of threads requested by the user (OMP_NUM_THREADS / default) and
	//! applies the deterministic mode. Called by FECore::InitModule; safe to call repeatedly.
	FECORE_API void Initialize();

	//! The number of threads requested by the user. Parallel linear solvers use this.
	FECORE_API int RequestedThreads();

	//! Switch deterministic mode on/off.
	FECORE_API void SetDeterministic(bool b);

	//! Is deterministic mode on?
	FECORE_API bool IsDeterministic();

	//! Deterministic parallel element loop with ordered assembly.
	//!
	//! compute(i, item) runs in parallel (all requested threads, also in deterministic mode)
	//! and must only write into its own item. assemble(i, item) is then called on the
	//! calling thread strictly in index order. All floating-point accumulation into shared
	//! data (global matrix, residual) therefore happens in exactly the order of a serial loop:
	//! the result is bitwise identical to a single-threaded run, for any thread count.
	template <class T, class Compute, class Assemble>
	void ComputeThenAssembleInOrder(int n, Compute compute, Assemble assemble, int chunk = 1024)
	{
		const int nt = RequestedThreads();
		std::vector<T> buf;
		for (int c0 = 0; c0 < n; c0 += chunk)
		{
			const int c1 = (c0 + chunk < n ? c0 + chunk : n);
			buf.resize(c1 - c0);
			#pragma omp parallel for num_threads(nt) schedule(static)
			for (int i = c0; i < c1; ++i) compute(i, buf[i - c0]);
			for (int i = c0; i < c1; ++i) assemble(i, buf[i - c0]);
		}
	}

	//! Collects error messages from a parallel loop and reports them in index order, so that
	//! the log output does not depend on thread timing either.
	class OrderedMessages
	{
	public:
		void add(int i, const std::string& msg);	//!< thread safe
		template <class F> void flush(F f) { for (auto& m : m_msg) f(m.second); m_msg.clear(); }
		bool empty() const { return m_msg.empty(); }
	private:
		std::map<int, std::string> m_msg;
	};
}
