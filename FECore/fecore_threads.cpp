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
#include "fecore_threads.h"
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
	bool s_initialized = false;
	bool s_deterministic = true;
	int  s_requested = 1;

	void initOnce()
	{
		if (s_initialized) return;
		s_initialized = true;
#ifdef _OPENMP
		s_requested = omp_get_max_threads();
#endif
		if (s_requested < 1) s_requested = 1;
		const char* sz = getenv("FEBIO_DETERMINISTIC");
		if (sz && ((strcmp(sz, "0") == 0) || (strcmp(sz, "false") == 0) || (strcmp(sz, "off") == 0)))
			s_deterministic = false;
	}

	void apply()
	{
#ifdef _OPENMP
		omp_set_num_threads(s_deterministic ? 1 : s_requested);
#endif
	}
}

void FECoreThreads::Initialize()
{
	if (s_initialized) return;
	initOnce();
	apply();
}

int FECoreThreads::RequestedThreads()
{
	initOnce();
	return s_requested;
}

void FECoreThreads::SetDeterministic(bool b)
{
	initOnce();
	s_deterministic = b;
	apply();
}

int FECoreThreads::ThreadNum()
{
#ifdef _OPENMP
	return omp_get_thread_num();
#else
	return 0;
#endif
}

int FECoreThreads::TeamSize()
{
#ifdef _OPENMP
	return omp_get_num_threads();
#else
	return 1;
#endif
}

bool FECoreThreads::StaticAssembly()
{
	static const bool b = [] { const char* s = getenv("FEBIO_ASM_STATIC"); return s && (s[0] == '1'); }();
	return b;
}

bool FECoreThreads::IsDeterministic()
{
	initOnce();
	return s_deterministic;
}

void FECoreThreads::OrderedMessages::add(int i, const std::string& msg)
{
#pragma omp critical (FECoreThreads_OrderedMessages)
	m_msg[i] = msg;
}
