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
	//! and must only write into its own item. assemble(i, item) is called strictly in index
	//! order, by one thread at a time. All floating-point accumulation into shared data (global
	//! matrix, residual) therefore happens in exactly the order of a serial loop: the result is
	//! bitwise identical to a single-threaded run, for any thread count.
	//!
	//! The loop is pipelined over chunks: while one thread assembles chunk k-1, the other
	//! threads compute chunk k. Two chunk buffers are reused across calls (per item type), so
	//! the items keep their allocations.
	template <class T, class Compute, class Assemble>
	void ComputeThenAssembleInOrder(int n, Compute compute, Assemble assemble, int chunk = 256)
	{
		if (n <= 0) return;
		const int nt = RequestedThreads();
		static thread_local std::vector<T> tls[2];
		std::vector<T>* buf = tls;	// the calling thread's buffers, shared with the team below
		if (nt <= 1)
		{
			if (buf[0].empty()) buf[0].resize(1);
			for (int i = 0; i < n; ++i) { compute(i, buf[0][0]); assemble(i, buf[0][0]); }
			return;
		}
		const int nchunks = (n + chunk - 1) / chunk;
		for (int k = 0; k < 2; ++k) if ((int)buf[k].size() < chunk) buf[k].resize(chunk);
		#pragma omp parallel num_threads(nt)
		{
			for (int c = 0; c <= nchunks; ++c)
			{
				// one thread assembles the previous chunk, in index order ...
				if (c > 0)
				{
					#pragma omp single nowait
					{
						const int a0 = (c - 1) * chunk, a1 = (a0 + chunk < n ? a0 + chunk : n);
						std::vector<T>& b = buf[(c - 1) & 1];
						for (int i = a0; i < a1; ++i) assemble(i, b[i - a0]);
					}
				}
				// ... while the others compute this one (the assembling thread joins when done).
				// The implicit barrier of the loop ends both before the buffers are swapped.
				if (c < nchunks)
				{
					const int c0 = c * chunk, c1 = (c0 + chunk < n ? c0 + chunk : n);
					std::vector<T>& b = buf[c & 1];
					#pragma omp for schedule(dynamic, 4)
					for (int i = c0; i < c1; ++i) compute(i, b[i - c0]);
				}
			}
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
