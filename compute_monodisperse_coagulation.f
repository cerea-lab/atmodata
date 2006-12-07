c     Copyright (C) 2006 CEREA (ENPC)
c
c     Authors: Vivien Mallet, Edouard Debry
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
c     http://spacetown.free.fr/lib/atmodata


c     Function: compute_coagulation_free_transition
c
c     Computes coagulation kernels for monodispersed
c     aerosols in the free transition regime.
c     2005/3/23: cleaning (Bruno Sportisse).
c
c     Parameters:
c     dp - aerosol diameter (�m).
c     vmp - mean particle velocity (m/s).
c     stick - sticking probability 0< <1 ().
c     cdifp - diffusion coefficient (m^2/s).
c     deltap - particle Knudsen number (�m).
c
c     Returns:
c     kercg - coagulation kernel (m^3/s).
      subroutine compute_coagulation_free_transition(dp, cdifp, deltap,
     $     vmp, stick, kercg)

      double precision pi
      parameter(pi = 3.141592653589d0)

      double precision dp, cdifp, deltap
      double precision vmp, stick, kercg

      double precision beta, dpp

      dpp = dp * 1.d-06         ! convert �m to m

      beta = 1.d0 /
     &     ( dp / (dp + deltap) ! adim
     &     + 8.d0 * stick       ! adim
     &     * cdifp              ! m2.s - 1
     &     / vmp                ! m.s - 1
     &     / dpp )              ! m
      
      kercg = 8.d0 * pi         ! adim
     &     * cdifp              ! m2.s - 1
     &     * dpp                ! m
     &     * beta               ! adim
      
      end
      
      
c     Function: compute_coagulation_free_molecular
c
c     Computes coagulation kernels for monodispersed
c     aerosols in the free molecular regime.
c     2005/3/23: cleaning (Bruno Sportisse).
c
c     Parameters:
c     dp - aerosol diameter (�m).
c     vmp - mean particle velocity (m/s).
c     stick - sticking probability 0< <1 ().
c
c     Returns:
c     kercg - coagulation kernel (m^3/s).
      subroutine compute_coagulation_free_molecular(dp, vmp, stick,
     $     kercg)

      double precision pi
      parameter(pi = 3.141592653589d0)

      double precision dp, vmp, stick, kercg

      double precision dpp

      dpp = dp * 1.d-06         ! convert �m to m
      
      kercg = pi                ! adim
     &     * dpp * dpp          ! m2
     &     * vmp                ! m.s - 1
     &     * stick              ! adim

      end


c     Function: compute_coagulation_continuous
c
c     Computes coagulation kernels for monodispersed
c     aerosols in the continuous regime.
c     2005/3/23: cleaning (Bruno Sportisse).
c
c     Parameters:
c     dp - aerosol diameter (�m).
c     cdifp - diffusion coefficient (m^2/s).
c
c     Returns:
c     kercg - coagulation kernel (m^3/s).
      subroutine compute_coagulation_continuous(dp, cdifp, kercg)

      double precision pi
      parameter(pi = 3.141592653589d0)

      double precision dp, cdifp, kercg

      double precision dpp

      dpp = dp * 1.d-06         ! convert �m to m
      
      kercg = 8.d0 * pi         ! adim
     &     * cdifp              ! m2.s - 1
     &     * dpp                ! m

      end
