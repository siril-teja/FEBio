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
}
