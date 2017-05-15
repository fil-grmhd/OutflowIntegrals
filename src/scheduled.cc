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

  if(!compute_global_quantities)
    return;

  // Anything to compute this iteration?
  if(cctk_iteration % integrate_every == 0)
    if(verbose)
      CCTK_INFO("Not integrating anything this iteration, skipping computation.");
    return;

  #pragma omp parallel for schedule(static)
  for(int k = cctk_nghostzones[2]; k < cctk_lsh[2]-cctk_nghostzones[2]; ++k) {
    for(int j = cctk_nghostzones[1]; j < cctk_lsh[1]-cctk_nghostzones[1]; ++j) {
      for(int i = cctk_nghostzones[0]; i < cctk_lsh[0]-cctk_nghostzones[0]; ++i) {
        const int ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
        if(abs(u_t)>1) {
          outint_vol_mass[ijk] = 1;
          outint_vol_energy_tot[ijk] = 1;
          outint_vol_energy_int[ijk] = 1;
        }
        else {
          outint_vol_mass[ijk] = 0;
          outint_vol_energy_tot[ijk] = 0;
          outint_vol_energy_int[ijk] = 0;
        }
      }
    }
  }
}
