//  Copyright (C) 2016, Ludwig Jens Papenfort <papenfort@th.physik.uni-frankfurt.de>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "types.hh"

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <limits>


extern "C" void am_computePointwise(CCTK_ARGUMENTS) {
  /*
   */
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!compute_geodesic && !compute_bernoulli)
    return;

  // Anything to compute this iteration?
  if(cctk_iteration % integrate_every == 0)
    if(verbose)
      CCTK_INFO("Not integrating anything this iteration, skipping computation.");
    return;

  gf_size = cctkGH->cctk_ash[0]*cctkGH->cctk_ash[1]*cctkGH->cctk_ash[2];
  CCTK_REAL* velx = &vel[0*gf_size];
  CCTK_REAL* vely = &vel[1*gf_size];
  CCTK_REAL* velz = &vel[2*gf_size];


  #pragma omp parallel for schedule(static)
  for(int k = cctk_nghostzones[2]; k < cctk_lsh[2]-cctk_nghostzones[2]; ++k) {
    for(int j = cctk_nghostzones[1]; j < cctk_lsh[1]-cctk_nghostzones[1]; ++j) {
      for(int i = cctk_nghostzones[0]; i < cctk_lsh[0]-cctk_nghostzones[0]; ++i) {
        const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

        const CCTK_REAL v_x = gxx[ijk]*velx[ijk] + gxy[ijk]*vely[ijk] + gxz[ijk]*velz[ijk];
        const CCTK_REAL v_y = gxy[ijk]*velx[ijk] + gyy[ijk]*vely[ijk] + gyz[ijk]*velz[ijk];
        const CCTK_REAL v_z = gxz[ijk]*velx[ijk] + gyz[ijk]*vely[ijk] + gzz[ijk]*velz[ijk];

        const CCTK_REAL u_t = w_lorentz[ijk]*(v_x*betax[ijk] + v_y*betay[ijk] + v_z*betaz[ijk] - alp[ijk]);
        const CCTK_REAL h = 1.0 + eps[ijk] + press[ijk]/rho[ijk];

        if(std:abs(u_t) > 1) {
          // mass
          pointwise_terms_geodesic[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = dens[ijk];
          // total energy
          pointwise_terms_geodesic[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = tau[ijk];
          // internal energy
          pointwise_terms_geodesic[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = dens[ijk]*eps[ijk];
        }
        else {
          pointwise_terms_geodesic[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0;
          pointwise_terms_geodesic[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = 0;
          pointwise_terms_geodesic[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = 0;
        }

        if(std:abs(h*u_t) > 1) {
          // mass
          pointwise_terms_bernoulli[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = dens[ijk];
          // total energy
          pointwise_terms_bernoulli[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = tau[ijk];
          // internal energy
          pointwise_terms_bernoulli[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = dens[ijk]*eps[ijk];
        }
        else {
          pointwise_terms_bernoulli[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0;
          pointwise_terms_bernoulli[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = 0;
          pointwise_terms_bernoulli[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = 0;
        }
      }
    }
  }
}
