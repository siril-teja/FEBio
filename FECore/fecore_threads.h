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
#include <atomic>
#include <memory>
#include <algorithm>

//! Thread configuration and run-to-run determinism.
//!
//! FEBio's element loops assemble the global matrix and vectors from several OpenMP
//! threads with atomic additions, so the order of floating-point additions (and hence the
//! last bits of the result) changes from run to run. In deterministic mode (the default)
//! FEBio's own OpenMP loops run on one thread, which makes every run bitwise reproducible.
//! Linear solvers that are deterministic themselves (the Rust solver, MKL PARDISO in
//! conditional-numerical-reproducibility mode) still use all requested threads.
//! The element loops of the solid, three-field and biphasic domains also use all threads:
//! they assemble color by color (ComputeAndAssembleColored), which is bitwise reproducible
//! for any thread count.
//!
//! Deterministic mode can be switched off with the -nondeterministic command line option
//! or the environment variable FEBIO_DETERMINISTIC=0.
namespace FECoreThreads {

	//! Records the number of threads requested by the user (OMP_NUM_THREADS / default) and
	//! applies the deterministic mode. Called by FECore::InitModule; safe to call repeatedly.
	FECORE_API void Initialize();

	//! The number of threads requested by the user. Parallel linear solvers use this.
	FECORE_API int RequestedThreads();

	//! Thread number and team size of the calling thread (0 and 1 outside a parallel region).
	FECORE_API int ThreadNum();
	FECORE_API int TeamSize();

	//! Temporary experiment switch for the colored assembly (FEBIO_ASM_STATIC=1).
	FECORE_API bool StaticAssembly();

	//! Switch deterministic mode on/off.
	FECORE_API void SetDeterministic(bool b);

	//! Is deterministic mode on?
	FECORE_API bool IsDeterministic();

	//! CPU hint for spin-wait loops.
	inline void SpinPause()
	{
#if defined(__x86_64__) || defined(__i386__)
		__builtin_ia32_pause();
#endif
	}

	//! Deterministic greedy coloring of a list of items (elements) for parallel assembly.
	//!
	//! The items are grouped into blocks of consecutive items, and the blocks are colored by
	//! the keys (nodes) their items write to: two blocks that share a key never get the same
	//! color. Blocks are colored first-fit in index order, so the coloring depends only on the
	//! connectivity (not on the number of threads or on timing). Consecutive items keep the
	//! memory locality of a plain loop.
	class Coloring
	{
	public:
		//! keys(i, out) appends the keys (0 <= key < nkeys) that item i writes to.
		//! The block size depends on n only (never on the thread count): 16, or less for
		//! small loops so that they still have enough blocks to run in parallel.
		template <class Keys> void Build(int n, int nkeys, Keys keys)
		{
			m_n = (n > 0 ? n : 0);
			m_bs = m_n / 64;
			if (m_bs < 1) m_bs = 1;
			if (m_bs > 16) m_bs = 16;
			m_ptr.assign(1, 0);
			m_block.clear();
			if (m_n == 0) return;
			const int nb = (m_n + m_bs - 1) / m_bs;
			std::vector<int> color(nb), k;
			for (int W = 1; ; W *= 2)	// colors in use, in units of 64
			{
				std::vector<unsigned long long> mask((size_t)nkeys * W, 0ull);
				int ncol = 0;
				bool ok = true;
				for (int b = 0; b < nb; ++b)
				{
					k.clear();
					for (int i = b * m_bs, i1 = (i + m_bs < m_n ? i + m_bs : m_n); i < i1; ++i) keys(i, k);
					int c = -1;
					for (int w = 0; (w < W) && (c < 0); ++w)
					{
						unsigned long long used = 0ull;
						for (int key : k) used |= mask[(size_t)key * W + w];
						if (~used) c = 64 * w + lowestZero(used);
					}
					if (c < 0) { ok = false; break; }	// more than 64*W colors: retry wider
					color[b] = c;
					for (int key : k) mask[(size_t)key * W + c / 64] |= (1ull << (c % 64));
					if (c + 1 > ncol) ncol = c + 1;
				}
				if (!ok) continue;
				m_ptr.assign(ncol + 1, 0);
				for (int b = 0; b < nb; ++b) m_ptr[color[b] + 1]++;
				for (int c = 0; c < ncol; ++c) m_ptr[c + 1] += m_ptr[c];
				m_block.resize(nb);
				std::vector<int> pos(m_ptr.begin(), m_ptr.end() - 1);
				for (int b = 0; b < nb; ++b) m_block[pos[color[b]]++] = b;

				// Dependencies, in processing order (color, then index): for each key, a block
				// waits for the previous block that writes to the key. Blocks that share a key
				// have different colors, so this orders the writes to every key by color, as
				// processing the colors one after another would.
				m_predPtr.assign(1, 0);
				m_pred.clear();
				std::vector<int> last(nkeys, -1);
				for (int p = 0; p < nb; ++p)
				{
					const int blk = m_block[p];
					k.clear();
					for (int i = blk * m_bs, i1 = (i + m_bs < m_n ? i + m_bs : m_n); i < i1; ++i) keys(i, k);
					const size_t first = m_pred.size();
					for (int key : k) if (last[key] >= 0) m_pred.push_back(last[key]);
					std::sort(m_pred.begin() + first, m_pred.end());
					m_pred.erase(std::unique(m_pred.begin() + first, m_pred.end()), m_pred.end());
					for (int key : k) last[key] = p;
					m_predPtr.push_back((int)m_pred.size());
				}
				return;
			}
		}

		int Items() const { return m_n; }
		int Colors() const { return (int)m_ptr.size() - 1; }
		int Blocks(int c) const { return m_ptr[c + 1] - m_ptr[c]; }
		int Blocks() const { return (int)m_block.size(); }
		//! the items [i0, i1) of the block at position p of the processing order
		void OrderedBlockRange(int p, int& i0, int& i1) const
		{
			i0 = m_block[p] * m_bs;
			i1 = (i0 + m_bs < m_n ? i0 + m_bs : m_n);
		}
		//! the positions (in processing order) of the blocks that block p must wait for
		const int* Preds(int p, int& n) const { n = m_predPtr[p + 1] - m_predPtr[p]; return m_pred.data() + m_predPtr[p]; }
		int BlockSize() const { return m_bs; }
		//! the block ids of color c (ascending)
		const int* BlockIds(int c) const { return m_block.data() + m_ptr[c]; }
		//! the items [i0, i1) of the j-th block of color c
		void BlockRange(int c, int j, int& i0, int& i1) const
		{
			i0 = m_block[m_ptr[c] + j] * m_bs;
			i1 = (i0 + m_bs < m_n ? i0 + m_bs : m_n);
		}

	private:
		static int lowestZero(unsigned long long used)
		{
			int b = 0;
			while (used & 1ull) { used >>= 1; ++b; }
			return b;
		}

	private:
		int m_n = 0;				//!< number of items
		int m_bs = 1;				//!< block size
		std::vector<int> m_ptr;		//!< color c holds the blocks m_block[m_ptr[c] .. m_ptr[c+1])
		std::vector<int> m_block;	//!< block indices, ordered by color and index (processing order)
		std::vector<int> m_predPtr;	//!< block p waits for m_pred[m_predPtr[p] .. m_predPtr[p+1])
		std::vector<int> m_pred;
	};

	//! Deterministic parallel element loop with colored assembly.
	//!
	//! The colors are processed one after another. Within a color the blocks run in parallel;
	//! each block runs compute(i, item) and assemble(i, item) for its items in index order.
	//! Blocks of one color write to disjoint global entries, so every global entry receives its
	//! contributions in a fixed order (color, then index). The result is therefore bitwise
	//! identical for any thread count and from run to run. With serial = true (shared writes
	//! the coloring does not cover, e.g. linear constraints), or for a small loop, the blocks
	//! are processed on one thread in the same order.
	template <class T, class Compute, class Assemble>
	void ComputeAndAssembleColored(const Coloring& col, Compute compute, Assemble assemble, bool serial = false)
	{
		const int n = col.Items();
		if (n <= 0) return;
		const int nt = RequestedThreads();
		const int ncol = col.Colors();
		if (serial || (nt <= 1) || (n < 8))
		{
			static thread_local T it;
			for (int c = 0; c < ncol; ++c)
				for (int j = 0, m = col.Blocks(c); j < m; ++j)
				{
					int i0, i1; col.BlockRange(c, j, i0, i1);
					for (int i = i0; i < i1; ++i) { compute(i, it); assemble(i, it); }
				}
			return;
		}
		// The threads take the blocks in processing order (color, then index). A block starts
		// once the blocks it depends on (earlier colors, shared keys) are done; those were
		// taken earlier, so the wait always ends. Every global entry thus receives its
		// contributions in color order, as in the serial loop above.
		const int nb = col.Blocks();
		std::unique_ptr<std::atomic<int>[]> done(new std::atomic<int>[nb]);
		for (int p = 0; p < nb; ++p) done[p].store(0, std::memory_order_relaxed);
		std::atomic<int> next(0);
		const bool bstatic = StaticAssembly();
		#pragma omp parallel num_threads(nt)
		{
			static thread_local T it;
			if (bstatic)
			{
				const int team = TeamSize(), tid = ThreadNum();
				const int slice = (n + team - 1) / team;
				for (int p = 0; p < nb; ++p)
				{
					int i0, i1; col.OrderedBlockRange(p, i0, i1);
					if (i0 / slice != tid) continue;
					int np = 0;
					const int* pred = col.Preds(p, np);
					for (int q = 0; q < np; ++q)
						while (done[pred[q]].load(std::memory_order_acquire) == 0) SpinPause();
					for (int i = i0; i < i1; ++i) { compute(i, it); assemble(i, it); }
					done[p].store(1, std::memory_order_release);
				}
			}
			else
			for (int p; (p = next.fetch_add(1, std::memory_order_relaxed)) < nb; )
			{
				int np = 0;
				const int* pred = col.Preds(p, np);
				for (int q = 0; q < np; ++q)
					while (done[pred[q]].load(std::memory_order_acquire) == 0) SpinPause();
				int i0, i1; col.OrderedBlockRange(p, i0, i1);
				for (int i = i0; i < i1; ++i) { compute(i, it); assemble(i, it); }
				done[p].store(1, std::memory_order_release);
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
