#//  Copyright (C) 2016, Ludwig Jens Papenfort <papenfort@th.physik.uni-frankfurt.de>
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


extern "C" void outint_computePointwise(CCTK_ARGUMENTS) {
  /*
   */
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // stop here if there's nothing to compute
  if((!compute_geodesic && !compute_bernoulli) || compute_every == 0) {
    if(verbose)
      CCTK_VInfo(CCTK_THORNSTRING,"Always skipping computation.");
    return;
  }

  // Anything to compute this iteration?
  if(cctk_iteration % compute_every != 0) {
    if(verbose) {
      CCTK_VInfo(CCTK_THORNSTRING,"Skipping computation. (it=%i)",cctk_iteration);
    }
    return;
  }

  if(verbose)
    CCTK_VInfo(CCTK_THORNSTRING,"Updating volume integral GFs... (it=%i)",cctk_iteration);

  // GF size and three velocity pointers
  CCTK_INT gf_size = cctkGH->cctk_ash[0]*cctkGH->cctk_ash[1]*cctkGH->cctk_ash[2];
  CCTK_REAL* velx = &vel[0*gf_size];
  CCTK_REAL* vely = &vel[1*gf_size];
  CCTK_REAL* velz = &vel[2*gf_size];

  bool const do_mhd = CCTK_IsThornActive("MHD_Analysis");
  
  CCTK_REAL * magnetic_energy_temp = nullptr;
  if(do_mhd) {
      magnetic_energy_temp = 
        static_cast<CCTK_REAL*>(CCTK_VarDataPtr(cctkGH,0,"MHD_Analysis::magnetic_energy_temp"));
      assert(magnetic_energy_temp != nullptr);
  }


  // loop over local grid extend
  #pragma omp parallel for schedule(static)
  for(int k = 0; k < cctk_lsh[2]; ++k) {
    for(int j = 0; j < cctk_lsh[1]; ++j) {
      for(int i = 0; i < cctk_lsh[0]; ++i) {
        // get compressed index
        const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

        // compute Jacobian for non-flat volume element
        CCTK_REAL dV = std::sqrt(utils::metric::spatial_det(gxx[ijk],
                                                            gxy[ijk],
                                                            gxz[ijk],
                                                            gyy[ijk],
                                                            gyz[ijk],
                                                            gzz[ijk]));

        // transform to covariant three velocity components
        const CCTK_REAL v_x = gxx[ijk]*velx[ijk] + gxy[ijk]*vely[ijk] + gxz[ijk]*velz[ijk];
        const CCTK_REAL v_y = gxy[ijk]*velx[ijk] + gyy[ijk]*vely[ijk] + gyz[ijk]*velz[ijk];
        const CCTK_REAL v_z = gxz[ijk]*velx[ijk] + gyz[ijk]*vely[ijk] + gzz[ijk]*velz[ijk];

        // get u_t
        const CCTK_REAL u_t = w_lorentz[ijk]*(v_x*betax[ijk] + v_y*betay[ijk] + v_z*betaz[ijk] - alp[ijk]);
        outint_ut[ijk] = u_t;

        // get the enthalpy
        const CCTK_REAL h = 1.0 + eps[ijk] + press[ijk]/rho[ijk];
        outint_h[ijk] = h;

        // precompute conserved density and energy
        const CCTK_REAL D = w_lorentz[ijk]*rho[ijk];
        const CCTK_REAL tau = D*h*w_lorentz[ijk] - press[ijk] - D;

        // geodesic criterion
        if(std::abs(u_t) > 1) {
          //CCTK_VInfo(CCTK_THORNSTRING,"Found unbound matter (geodesic) (it=%i,x=%e,y=%e,z=%e,u_t=%e)",cctk_iteration,x[ijk],y[ijk],z[ijk],u_t);

          // mass
          outint_terms_geo[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = dV*D;
          // total energy
          outint_terms_geo[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = dV*tau;
          // internal energy
          outint_terms_geo[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = dV*D*eps[ijk];
          // Ye_star
          outint_terms_geo[CCTK_GFINDEX4D(cctkGH,i,j,k,3)] = dV*D*Y_e[ijk];
          // magnetic_energy
          outint_terms_geo[CCTK_GFINDEX4D(cctkGH,i,j,k,4)] = (do_mhd) ? magnetic_energy_temp[ijk] : 0.;

        }
        else {
          outint_terms_geo[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0;
          outint_terms_geo[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = 0;
          outint_terms_geo[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = 0;
          outint_terms_geo[CCTK_GFINDEX4D(cctkGH,i,j,k,3)] = 0;
          outint_terms_geo[CCTK_GFINDEX4D(cctkGH,i,j,k,4)] = 0;
        }

        // bernoulli criterion
        if(std::abs(h*u_t) > 1) {
          //CCTK_VInfo(CCTK_THORNSTRING,"Found unbound matter (bernoulli) (it=%i,x=%e,y=%e,z=%e,h=%e,u_t=%e)",cctk_iteration,x[ijk],y[ijk],z[ijk],h,u_t);

          // mass
          outint_terms_bern[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = dV*D;
          // total energy
          outint_terms_bern[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = dV*tau;
          // internal energy
          outint_terms_bern[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = dV*D*eps[ijk];
          // Ye_star
          outint_terms_bern[CCTK_GFINDEX4D(cctkGH,i,j,k,3)] = dV*D*Y_e[ijk];
          // magnetic_energy
          outint_terms_bern[CCTK_GFINDEX4D(cctkGH,i,j,k,4)] = (do_mhd) ? magnetic_energy_temp[ijk] : 0.;
        }
        else {
          outint_terms_bern[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0;
          outint_terms_bern[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = 0;
          outint_terms_bern[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] = 0;
          outint_terms_bern[CCTK_GFINDEX4D(cctkGH,i,j,k,3)] = 0;
          outint_terms_bern[CCTK_GFINDEX4D(cctkGH,i,j,k,4)] = 0;
        }
      }
    }
  }
}
