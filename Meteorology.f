c     Copyright (C) 2006 CEREA (ENPC)
c     
c     Author: Vivien Mallet
c     
c     This file is part of AtmoData library. AtmoData library is a tool
c     for data processing in atmospheric sciences.
c     
c     AtmoData is free software; you can redistribute it and/or modify
c     it under the terms of the GNU General Public License as published
c     by the Free Software Foundation; either version 2 of the License,
c     or (at your option) any later version.
c     
c     AtmoData is distributed in the hope that it will be useful, but
c     WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c     General Public License (file "LICENSE") for more details.
c     
c     For more information, please see the AtmoData home page:
c     http://www.enpc.fr/cerea/atmodata/


c     Function: compute_relative_humidity
c     Computes relative humidity from specific humidity.
c
c     Parameters:
c     specific_humidity - Specific humidity (kg/kg).
c     temperature - Temperature (K).
c     pressure - Pressure (Pa).
c
c     Returns:
c     relative_humidity - relative humidity.
      subroutine compute_relative_humidity(specific_humidity,
     $     temperature, pressure, relative_humidity)

      double precision specific_humidity
      double precision temperature
      double precision pressure

      double precision relative_humidity

      double precision pressure_sat

      pressure_sat = 611.2d0 * dexp(17.67d0 * (temperature - 273.15d0) /
     $     (temperature - 29.65d0))

      relative_humidity = specific_humidity * pressure / ( (0.62197d0
     $     * (1.d0 - specific_humidity) + specific_humidity)
     $     * pressure_sat)

      end
