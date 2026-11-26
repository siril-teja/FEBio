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
#include "FEElasticFiberMaterial.h"
#include "FEFiberMaterialPoint.h"
#include <FECore/FEConstValueVec3.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEElasticFiberMaterial, FEElasticMaterial)
	ADD_PROPERTY(m_fiber, "fiber");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMaterialPointData* FEElasticFiberMaterial::CreateMaterialPointData()
{
    FEFiberMaterialPoint* mp = new FEFiberMaterialPoint(new FEElasticMaterialPoint);
    return mp;
}

//-----------------------------------------------------------------------------
FEElasticFiberMaterial::FEElasticFiberMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_fiber = nullptr;

    m_Us = mat3dd(1);
    m_bUs = false;
}

//-----------------------------------------------------------------------------
// get the fiber vector (in global coordinates)
bool FEElasticFiberMaterial::Init()
{
	if (FEElasticMaterial::Init() == false) return false;
	m_cacheFiber = (m_fiber != nullptr) && m_fiber->isReferenceInvariant() && LocalCSIsReferenceInvariant();
	return true;
}

//-----------------------------------------------------------------------------
vec3d FEElasticFiberMaterial::FiberVector(FEMaterialPoint& mp)
{
	// The reference fiber direction Q*fiber is constant for reference-invariant
	// material axes and fiber generators: compute it once per material point.
	// (memo keyed by this material and the point's reference position)
	if (m_cacheFiber && (mp.m_cacheOwner == this) && (mp.m_cacheR0 == mp.m_r0))
		return FiberPreStretch(mp.m_cacheVec);

	// get the material coordinate system
	mat3d Q = GetLocalCS(mp);

	// get the fiber vector in local coordinates
	vec3d fiber = m_fiber->unitVector(mp);

	// convert to global coordinates
	vec3d a0 = Q*fiber;
	if (m_cacheFiber) { mp.m_cacheOwner = this; mp.m_cacheR0 = mp.m_r0; mp.m_cacheVec = a0; }

	// account for prior deformation in multigenerational formulation
	return FiberPreStretch(a0);
}

//-----------------------------------------------------------------------------
vec3d FEElasticFiberMaterial::FiberPreStretch(const vec3d a0)
{
    // account for prior deformation in multigenerational formulation
    if (m_bUs) {
        vec3d a = (m_Us*a0);
        a.unit();
        return a;
    }
    else
        return a0;
}
