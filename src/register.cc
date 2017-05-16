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

extern "C" void outint_registerIntegrals(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(compute_bernoulli) {
    SphericalIntegrator_RegisterVolumeIntegral("OutflowIntegrals::outint_terms_bern[0]",
                                               "OutflowIntegrals::bernoulli_vol_mass",
                                               sphere_id,
                                               compute_every);
    SphericalIntegrator_RegisterVolumeIntegral("OutflowIntegrals::outint_terms_bern[1]",
                                               "OutflowIntegrals::bernoulli_vol_energy_tot",
                                               sphere_id,
                                               compute_every);
    SphericalIntegrator_RegisterVolumeIntegral("OutflowIntegrals::outint_terms_bern[2]",
                                               "OutflowIntegrals::bernoulli_vol_energy_int",
                                               sphere_id,
                                               compute_every);
  }
  if(compute_geodesic) {
    SphericalIntegrator_RegisterVolumeIntegral("OutflowIntegrals::outint_terms_geo[0]",
                                               "OutflowIntegrals::geodesic_vol_mass",
                                               sphere_id,
                                               compute_every);
    SphericalIntegrator_RegisterVolumeIntegral("OutflowIntegrals::outint_terms_geo[1]",
                                               "OutflowIntegrals::geodesic_vol_energy_tot",
                                               sphere_id,
                                               compute_every);
    SphericalIntegrator_RegisterVolumeIntegral("OutflowIntegrals::outint_terms_geo[2]",
                                               "OutflowIntegrals::geodesic_vol_energy_int",
                                               sphere_id,
                                               compute_every);
  }
}

