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

extern "C" void am_registerIntegrals(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(compute_global_quantities) {
    SphericalIntegrator_RegisterVolumeIntegral("OutflowIntegrals::outint_terms_threshold[0]",
                                               "OutflowIntegrals::outint_vol_mass",
                                               sphere_id,
                                               integrate_every);
    SphericalIntegrator_RegisterVolumeIntegral("OutflowIntegrals::outint_terms_threshold[1]",
                                               "OutflowIntegrals::outint_vol_energy_tot",
                                               sphere_id,
                                               integrate_every);
    SphericalIntegrator_RegisterVolumeIntegral("OutflowIntegrals::outint_terms_threshold[2]",
                                               "OutflowIntegrals::outint_vol_energy_int",
                                               sphere_id,
                                               integrate_every);
  }
}

