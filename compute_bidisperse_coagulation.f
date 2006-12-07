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


c     Function: compute_bidisperse_coagulation_kernel
c
c     Computes coagulation kernels for bidispersed aerosols.
c     2005/3/23: cleaning (Bruno Sportisse).
c
c     Parameters:
c     temp - temperature (Kelvin).
c     airfmp - air free mean path (�m).
c     d1 - diameter of the first coagulating particle (�m).
c     d2 - diameter of the second coagulating particle (�m).
c     m1 - mass of the first coagulating particle (�g).
c     m2 - mass of the second coagulating particle (�g).
c
c     Returns:
c     kercg - coagulation kernel (m^3/s).
      subroutine compute_bidisperse_coagulation_kernel(temp, airfmp,
     $     d1, d2, m1, m2, kercg)

      double precision temp, airfmp
      double precision d1, m1, d2, m2, kercg

      double precision pi, frac3, dmin, mmin, knmin
      double precision stick, kb, muair, knmax

      parameter(pi = 3.141592653589d0,
     $     frac3 = 0.333333333333d0,
     $     kb = 1.381d - 23,    ! boltzmann constant j.k - 1
     $     muair = 1.725d - 05, ! air dynamic viscosity kg.m - 1.s - 1
     $     stick = 1.d0,        ! sticking probability 0< <1
     $     dmin = 1.d-03,       ! min aero diam �m
     $     mmin = 1.d-17,       ! min aero mass �g
     $     knmin = 0.02d0,      ! min threshold of knudsen number
     $     knmax = 10.d0)

c     Mmin corresponds roughly to an aerosol of diameter Dmin and a
c     specific mass equal to 1.D-06 �g.�m-3.

      double precision kn1, kn2, knp
      double precision cdif1, cdif2, cdifp
      double precision delta1, delta2, deltap
      double precision vm1, vm2, vmp
      double precision cor1, cor2, l1, l2
      double precision dp, lambdap
      double precision d1c, d2c, m1c, m2c

c     eventually correct inputs
      d1c = dmax1(d1, dmin)
      d2c = dmax1(d2, dmin)
      m1c = dmax1(m1, mmin)
      m2c = dmax1(m2, mmin)
      
c     knudsen number for particles in air
      kn1 = 2.d0 * airfmp / d1c ! adim
      kn2 = 2.d0 * airfmp / d2c ! adim
      
c     diffusion coefficient
      cdif1 = Kb                ! J.K-1
     $     * temp               ! k
     $     / 3.d0 / pi          ! adim
     $     / muair              ! air viscosity kg.m - 1.s - 1
     $     / d1c                ! �m
     $     * 1.d06              ! convert �m - > m
      
      cor1 = ( 5.D0+( kn1*( 4.D0
     $     + 6.d0 * kn1 * ( 1.d0 + 3.d0 * kn1))) )
     $     / ( 5.d0 + kn1 * ((8.d0 + pi) * kn1 - 1.d0) )

      cdif1=cdif1*cor1          ! m2.s-1

      cdif2 = Kb                ! J.K-1
     $     * temp               ! k
     $     / 3.d0 / pi          ! adim
     $     / muair              ! air viscosity kg.m - 1.s - 1
     $     / d2c                ! �m
     $     * 1.d06              ! convert �m - > m
      

      cor2 = ( 5.D0+( kn2*( 4.D0
     $     + 6.d0 * kn2 * ( 1.d0 + 3.d0 * kn2))) )
     $     / ( 5.d0 + kn2 * ((8.d0 + pi) * kn2 - 1.d0) )
      
      cdif2 = cdif2 * cor2      ! m2.s-1

c     mean quadratic aerosol velocity  ! m.s-1
      vm1 = dsqrt( 8.d0 / pi    ! adim
     $     * kb                 ! j.k - 1
     $     * temp               ! k
     $     / m1c                ! �g
     $     * 1.d09 )            ! convert �g - > kg
      
      vm2 = dsqrt( 8.d0 / pi    ! adim
     $     * kb                 ! j.k - 1
     $     * temp               ! k
     $     / m2c                ! �g
     $     * 1.d09 )            ! convert �g - > kg

c     average values
      cdifp=(cdif1+cdif2)*5.D-01 ! m2.s-1
      vmp = dsqrt(vm1 * vm1 + vm2 * vm2) ! m.s - 1
      dp = (d1c + d2c) * 5.d-01 ! �m
      
      lambdap = 8.d0 / pi
     $     * cdifp              ! m2.s - 1
     $     / vmp                ! m.s - 1
     $     * 1.d06              ! convert m to �m
      
      knp = 2.d0
     $     * lambdap            ! �m
     $     / dp                 ! �m

      if (knp.le.knmin) then
         
         call compute_coagulation_continuous(dp, cdifp, kercg)
         
      elseif(knp.ge.knmax) then
         
         call compute_coagulation_free_molecular(dp, vmp, stick, kercg)

      else
                                ! free mean path of aerosols
         l1 = 8.d0 / pi         ! adim
     $        * cdif1           ! m2.s - 1
     $        / vm1             ! m.s - 1
     $        * 1.d06           ! convert m - > �m

         l2 = 8.d0 / pi         ! adim
     $        * cdif2           ! m2.s - 1
     $        / vm2             ! m.s - 1
     $        * 1.d06           ! convert m - > �m

                                ! aerosol knudsen number
         kn1 = l1 / d1c         ! adim
         kn2 = l2 / d2c         ! adim

         delta1 = d1c           ! �m
     $        * (((1.d0 + kn1) ** 3.d0
     $        - (1.d0 + kn1 * kn1) ** 1.5d0 )
     $        * frac3 / kn1 - 1.d0)
         
         delta2 = d2c           ! �m
     $        * (((1.d0 + kn2) ** 3.d0
     $        - (1.d0 + kn2 * kn2) ** 1.5d0 )
     $        * frac3 / kn2 - 1.d0)


         if ( delta1 * delta1 + delta2 * delta2.le.0.d0) then
            write(6, * )'coagulation kernel: sqrt(<0)'
            stop
         endif

         deltap = dsqrt( delta1 * delta1
     $        + delta2 * delta2 )
         
         call compute_coagulation_free_transition(dp, cdifp, deltap,
     $        vmp, stick, kercg)

      endif

      end