// Copyright (C) 2003-2005 CEREA
//     Author: Vivien Mallet
//
// CEREA (http://www.enpc.fr/cerea) is a joint laboratory of
// ENPC (http://www.enpc.fr) and EDF R&D (http://www.edf.fr).
//
// This file is part of AtmoData library.
// AtmoData library is a tool for data processing in atmospheric
// sciences.
//
// AtmoData is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// AtmoData is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License (file "license") for more details.
//
// For more information, please see the AtmoData home page:
//     http://www.enpc.fr/cerea/atmodata/


#ifndef ATMODATA_FILE_AEROSOL_HXX

// Fortran functions.
#if defined(__GNUG__) && __GNUG__ < 4 && !defined (__INTEL_COMPILER)

#define _compute_gas_diffusivity compute_gas_diffusivity__
#define _gerber_wet_diameter gerber_wet_diameter__
#define _compute_collision_integral compute_collision_integral__
#define _compute_condensation_transfer_rate \
compute_condensation_transfer_rate__
#define _compute_quadratic_mean_velocity compute_quadratic_mean_velocity__
#define _compute_saturation_concentration compute_saturation_concentration__
#define _compute_kelvin_coefficient compute_kelvin_coefficient__

#else

#define _compute_gas_diffusivity compute_gas_diffusivity_
#define _gerber_wet_diameter gerber_wet_diameter_
#define _compute_collision_integral compute_collision_integral_
#define _compute_condensation_transfer_rate \
compute_condensation_transfer_rate_
#define _compute_quadratic_mean_velocity compute_quadratic_mean_velocity_
#define _compute_saturation_concentration compute_saturation_concentration_
#define _compute_kelvin_coefficient compute_kelvin_coefficient_

#endif

extern "C"
{
  void _compute_gas_diffusivity(const double*, const double*, const double*,
				const double*, const double*, double*);
  void _gerber_wet_diameter(const double*, const double*, const double*,
			    const double*);
  void _compute_collision_integral(const double*, const double*);
  void _compute_condensation_transfer_rate(const double*, const double*,
					   const double*, const double*,
					   double*);
  void _compute_quadratic_mean_velocity(const double*, const double*,
					double*);
  void _compute_saturation_concentration(const double*, const double*,
					 const double*, const double*,
					 double*);
  void _compute_kelvin_coefficient(const double*, const double*,
				   const double*, const double*,
				   const double*, double*);
}

#define ATMODATA_FILE_AEROSOL_HXX
#endif
