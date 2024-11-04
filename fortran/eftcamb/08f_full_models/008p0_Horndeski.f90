!----------------------------------------------------------------------------------------
!
! This file is part of EFTCAMB.
!
! Copyright (C) 2013-2023 by the EFTCAMB authors
!
! The EFTCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file eftcamb/LICENSE at
! the top level of the EFTCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 007p0_Horndeski.f90
!! This file contains the definition of the full Horndeski model.


!----------------------------------------------------------------------------------------
!> This module contains the full Horndeski model.

!> @author Gen Ye

module EFTCAMB_Horndeski

   use precision
   use IniObjects
   use MpiUtils
   use FileUtils
   use constants, only : c, const_pi
   use equispaced_linear_interpolation_1D
   use EFTCAMB_cache
   use EFT_def
   use EFTCAMB_mixed_algorithms
   use EFTCAMB_rootfind
   use EFTCAMB_abstract_model_full
   use EFTCAMB_abstract_parametrizations_1D
   use EFTCAMB_power_law_parametrizations_1D
   use EFTCAMB_constant_parametrization_1D
   use EFTCAMB_taylor_series_1D
   use EFTCAMB_pade_series_1D
   use EFTCAMB_fourier_parametrizations_1D
   use EFTCAMB_exponential_parametrizations_2_1D
   use EFTCAMB_double_exponential_parametrizations_1D
   use EFTCAMB_cosine_parametrizations_1D
   use EFTCAMB_axion_parametrizations_1D
   use EFTCAMB_abstract_model_designer
   use MassiveNu
!    use Horndeski_coefficients
   use iso_c_binding

   implicit none

   interface

      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that computes the energy constraint E3 H^3 + E2 H^2 + E1 H + E0 = 0
      subroutine horndeski_hubble_coefficients(a,dphi,rhotot,G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx,E0,E1,E2,E3) bind(C)
         use iso_c_binding, only: c_double
         implicit none

         real(c_double), value :: a      !<- the scale factor
         real(c_double), value :: dphi   !<- conformal time derivative of the scalar field
         real(c_double), value :: rhotot !<- \f$ \rho_{matter} + p_{matter} \f$
         real(c_double), value :: G2     !<- Horndeski Lagrangian function \f$ G2 \f$
         real(c_double), value :: G2p    !<- Horndeski Lagrangian function \f$ G2_\phi \f$
         real(c_double), value :: G2x    !<- Horndeski Lagrangian function \f$ G2_X \f$
         real(c_double), value :: G2pp
         real(c_double), value :: G2px
         real(c_double), value :: G2xx
         real(c_double), value :: G2ppp
         real(c_double), value :: G2ppx
         real(c_double), value :: G2pxx
         real(c_double), value :: G2xxx
         real(c_double), value :: G3p
         real(c_double), value :: G3x
         real(c_double), value :: G3pp
         real(c_double), value :: G3px
         real(c_double), value :: G3xx
         real(c_double), value :: G3ppp
         real(c_double), value :: G3ppx
         real(c_double), value :: G3pxx
         real(c_double), value :: G3xxx
         real(c_double), value :: G4
         real(c_double), value :: G4p
         real(c_double), value :: G4x
         real(c_double), value :: G4pp
         real(c_double), value :: G4px
         real(c_double), value :: G4xx
         real(c_double), value :: G4ppp
         real(c_double), value :: G4ppx
         real(c_double), value :: G4pxx
         real(c_double), value :: G4xxx
         real(c_double), value :: G4pppp
         real(c_double), value :: G4pppx
         real(c_double), value :: G4ppxx
         real(c_double), value :: G4pxxx
         real(c_double), value :: G4xxxx
         real(c_double), value :: G5
         real(c_double), value :: G5p
         real(c_double), value :: G5x
         real(c_double), value :: G5pp
         real(c_double), value :: G5px
         real(c_double), value :: G5xx
         real(c_double), value :: G5ppp
         real(c_double), value :: G5ppx
         real(c_double), value :: G5pxx
         real(c_double), value :: G5xxx
         real(c_double), value :: G5pppp
         real(c_double), value :: G5pppx
         real(c_double), value :: G5ppxx
         real(c_double), value :: G5pxxx
         real(c_double), value :: G5xxxx

         real(c_double) :: E0
         real(c_double) :: E1
         real(c_double) :: E2
         real(c_double) :: E3

      end subroutine horndeski_hubble_coefficients

      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that computes the EFT functions alpha_i and a few derivatives
      subroutine horndeski_eft_alphas(a,H,dH,ddH,dphi,ddphi,dddphi,G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx,alphaM,dalphaM,alphaB,dalphaB,alphaK,dalphaK,alphaT,dalphaT,M,dM,ddM) bind(C)
         use iso_c_binding, only: c_double
         implicit none

         real(c_double), value :: a      !<- the scale factor
         real(c_double), value :: H      !<- the physcial Hubble parameter
         real(c_double), value :: dH     !<- conformal time derivative of the physical Hubble
         real(c_double), value :: ddH    !<- second conformal time derivative of the physical Hubble
         real(c_double), value :: dphi   !<- conformal time derivative of the scalar field
         real(c_double), value :: ddphi  !<- double conformal time derivative of the scalar field
         real(c_double), value :: dddphi !<- third conformal time derivative of the scalar field
         real(c_double), value :: G2     !<- Horndeski Lagrangian function \f$ G2 \f$
         real(c_double), value :: G2p    !<- Horndeski Lagrangian function \f$ G2_\phi \f$
         real(c_double), value :: G2x    !<- Horndeski Lagrangian function \f$ G2_X \f$
         real(c_double), value :: G2pp
         real(c_double), value :: G2px
         real(c_double), value :: G2xx
         real(c_double), value :: G2ppp
         real(c_double), value :: G2ppx
         real(c_double), value :: G2pxx
         real(c_double), value :: G2xxx
         real(c_double), value :: G3p
         real(c_double), value :: G3x
         real(c_double), value :: G3pp
         real(c_double), value :: G3px
         real(c_double), value :: G3xx
         real(c_double), value :: G3ppp
         real(c_double), value :: G3ppx
         real(c_double), value :: G3pxx
         real(c_double), value :: G3xxx
         real(c_double), value :: G4
         real(c_double), value :: G4p
         real(c_double), value :: G4x
         real(c_double), value :: G4pp
         real(c_double), value :: G4px
         real(c_double), value :: G4xx
         real(c_double), value :: G4ppp
         real(c_double), value :: G4ppx
         real(c_double), value :: G4pxx
         real(c_double), value :: G4xxx
         real(c_double), value :: G4pppp
         real(c_double), value :: G4pppx
         real(c_double), value :: G4ppxx
         real(c_double), value :: G4pxxx
         real(c_double), value :: G4xxxx
         real(c_double), value :: G5
         real(c_double), value :: G5p
         real(c_double), value :: G5x
         real(c_double), value :: G5pp
         real(c_double), value :: G5px
         real(c_double), value :: G5xx
         real(c_double), value :: G5ppp
         real(c_double), value :: G5ppx
         real(c_double), value :: G5pxx
         real(c_double), value :: G5xxx
         real(c_double), value :: G5pppp
         real(c_double), value :: G5pppx
         real(c_double), value :: G5ppxx
         real(c_double), value :: G5pxxx
         real(c_double), value :: G5xxxx

         real(c_double) :: alphaM          !<- The EFT function alpha_M
         real(c_double) :: dalphaM         !<- conformal time of the EFT function alpha_M
         real(c_double) :: alphaB          !<- The EFT function alpha_B
         real(c_double) :: dalphaB         !<- conformal time of the EFT function alpha_B
         real(c_double) :: alphaK          !<- The EFT function alpha_K
         real(c_double) :: dalphaK         !<- conformal time of the EFT function alpha_K
         real(c_double) :: alphaT          !<- The EFT function alpha_T
         real(c_double) :: dalphaT         !<- conformal time of the EFT function alpha_T
         real(c_double) :: M               !<- The effective planck mass squared
         real(c_double) :: dM              !<- conformal time of the effective planck mass squared
         real(c_double) :: ddM             !<- second conformal time of the effective planck mass squared

      end subroutine horndeski_eft_alphas


      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that computes the background EoM coefficients
      subroutine horndeski_backeq_coefficients(a,H,dphi,rptot,G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx,cp_1,ch_1,cs_1,cp_2,ch_2,cs_2) bind(C)
         use iso_c_binding, only: c_double
         implicit none

         real(c_double), value :: a      !<- the scale factor
         real(c_double), value :: H      !<- the physcial Hubble parameter
         real(c_double), value :: dphi   !<- conformal time derivative of the scalar field
         real(c_double), value :: rptot  !<- \f$ \rho_{matter} + p_{matter} \f$
         real(c_double), value :: G2     !<- Horndeski Lagrangian function \f$ G2 \f$
         real(c_double), value :: G2p    !<- Horndeski Lagrangian function \f$ G2_\phi \f$
         real(c_double), value :: G2x    !<- Horndeski Lagrangian function \f$ G2_X \f$
         real(c_double), value :: G2pp
         real(c_double), value :: G2px
         real(c_double), value :: G2xx
         real(c_double), value :: G2ppp
         real(c_double), value :: G2ppx
         real(c_double), value :: G2pxx
         real(c_double), value :: G2xxx
         real(c_double), value :: G3p
         real(c_double), value :: G3x
         real(c_double), value :: G3pp
         real(c_double), value :: G3px
         real(c_double), value :: G3xx
         real(c_double), value :: G3ppp
         real(c_double), value :: G3ppx
         real(c_double), value :: G3pxx
         real(c_double), value :: G3xxx
         real(c_double), value :: G4
         real(c_double), value :: G4p
         real(c_double), value :: G4x
         real(c_double), value :: G4pp
         real(c_double), value :: G4px
         real(c_double), value :: G4xx
         real(c_double), value :: G4ppp
         real(c_double), value :: G4ppx
         real(c_double), value :: G4pxx
         real(c_double), value :: G4xxx
         real(c_double), value :: G4pppp
         real(c_double), value :: G4pppx
         real(c_double), value :: G4ppxx
         real(c_double), value :: G4pxxx
         real(c_double), value :: G4xxxx
         real(c_double), value :: G5
         real(c_double), value :: G5p
         real(c_double), value :: G5x
         real(c_double), value :: G5pp
         real(c_double), value :: G5px
         real(c_double), value :: G5xx
         real(c_double), value :: G5ppp
         real(c_double), value :: G5ppx
         real(c_double), value :: G5pxx
         real(c_double), value :: G5xxx
         real(c_double), value :: G5pppp
         real(c_double), value :: G5pppx
         real(c_double), value :: G5ppxx
         real(c_double), value :: G5pxxx
         real(c_double), value :: G5xxxx

         real(c_double) :: cp_1   !<- EoM cp_i * phi'' + ch_i * H' + cs_i = 0, i=1,2
         real(c_double) :: ch_1
         real(c_double) :: cs_1
         real(c_double) :: cp_2
         real(c_double) :: ch_2
         real(c_double) :: cs_2

      end subroutine horndeski_backeq_coefficients


      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that computes the conformal time derivative of background EoM coefficients
      subroutine horndeski_backeq_coefficients_p(a,H,dH,dphi,ddphi,rptot,drptot,G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx,cp_1p,ch_1p,cs_1p,cp_2p,ch_2p,cs_2p) bind(C)
         use iso_c_binding, only: c_double
         implicit none

         real(c_double), value :: a      !<- the scale factor
         real(c_double), value :: H      !<- the physcial Hubble parameter
         real(c_double), value :: dH     !<- conformal time derivative of the physical Hubble
         real(c_double), value :: dphi   !<- conformal time derivative of the scalar field
         real(c_double), value :: ddphi  !<- double conformal time derivative of the scalar field
         real(c_double), value :: rptot  !<- \f$ \rho_{matter} + p_{matter} \f$
         real(c_double), value :: drptot !<- conformal time derivative of \f$ \rho_{matter} + p_{matter} \f$
         real(c_double), value :: G2     !<- Horndeski Lagrangian function \f$ G2 \f$
         real(c_double), value :: G2p    !<- Horndeski Lagrangian function \f$ G2_\phi \f$
         real(c_double), value :: G2x    !<- Horndeski Lagrangian function \f$ G2_X \f$
         real(c_double), value :: G2pp
         real(c_double), value :: G2px
         real(c_double), value :: G2xx
         real(c_double), value :: G2ppp
         real(c_double), value :: G2ppx
         real(c_double), value :: G2pxx
         real(c_double), value :: G2xxx
         real(c_double), value :: G3p
         real(c_double), value :: G3x
         real(c_double), value :: G3pp
         real(c_double), value :: G3px
         real(c_double), value :: G3xx
         real(c_double), value :: G3ppp
         real(c_double), value :: G3ppx
         real(c_double), value :: G3pxx
         real(c_double), value :: G3xxx
         real(c_double), value :: G4
         real(c_double), value :: G4p
         real(c_double), value :: G4x
         real(c_double), value :: G4pp
         real(c_double), value :: G4px
         real(c_double), value :: G4xx
         real(c_double), value :: G4ppp
         real(c_double), value :: G4ppx
         real(c_double), value :: G4pxx
         real(c_double), value :: G4xxx
         real(c_double), value :: G4pppp
         real(c_double), value :: G4pppx
         real(c_double), value :: G4ppxx
         real(c_double), value :: G4pxxx
         real(c_double), value :: G4xxxx
         real(c_double), value :: G5
         real(c_double), value :: G5p
         real(c_double), value :: G5x
         real(c_double), value :: G5pp
         real(c_double), value :: G5px
         real(c_double), value :: G5xx
         real(c_double), value :: G5ppp
         real(c_double), value :: G5ppx
         real(c_double), value :: G5pxx
         real(c_double), value :: G5xxx
         real(c_double), value :: G5pppp
         real(c_double), value :: G5pppx
         real(c_double), value :: G5ppxx
         real(c_double), value :: G5pxxx
         real(c_double), value :: G5xxxx

         real(c_double) :: cp_1p          !<- EoM cp_i * phi'' + ch_i * H' + cs_i = 0, i=1,2
         real(c_double) :: ch_1p
         real(c_double) :: cs_1p
         real(c_double) :: cp_2p
         real(c_double) :: ch_2p
         real(c_double) :: cs_2p

      end subroutine horndeski_backeq_coefficients_p


      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that computes the coefficients of perturbation equations to linear order in phidot
      subroutine horndeski_perturbeq_coefficients_lin(a,H,dphi,rptot,drptot,G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx,spi1,spi2,spi3,spih,spie,spim,s00,s00k,s00p,s0i,s0ip,sii,siik,siip,siipp,sij) bind(C)
         use iso_c_binding, only: c_double
         implicit none

         real(c_double), value :: a      !<- the scale factor
         real(c_double), value :: H      !<- the physcial Hubble parameter
         real(c_double), value :: dphi   !<- conformal time derivative of the scalar field
         real(c_double), value :: rptot  !<- \f$ \rho_{matter} + p_{matter} \f$
         real(c_double), value :: drptot !<- conformal time derivative of \f$ \rho_{matter} + p_{matter} \f$
         real(c_double), value :: G2     !<- Horndeski Lagrangian function \f$ G2 \f$
         real(c_double), value :: G2p    !<- Horndeski Lagrangian function \f$ G2_\phi \f$
         real(c_double), value :: G2x    !<- Horndeski Lagrangian function \f$ G2_X \f$
         real(c_double), value :: G2pp
         real(c_double), value :: G2px
         real(c_double), value :: G2xx
         real(c_double), value :: G2ppp
         real(c_double), value :: G2ppx
         real(c_double), value :: G2pxx
         real(c_double), value :: G2xxx
         real(c_double), value :: G3p
         real(c_double), value :: G3x
         real(c_double), value :: G3pp
         real(c_double), value :: G3px
         real(c_double), value :: G3xx
         real(c_double), value :: G3ppp
         real(c_double), value :: G3ppx
         real(c_double), value :: G3pxx
         real(c_double), value :: G3xxx
         real(c_double), value :: G4
         real(c_double), value :: G4p
         real(c_double), value :: G4x
         real(c_double), value :: G4pp
         real(c_double), value :: G4px
         real(c_double), value :: G4xx
         real(c_double), value :: G4ppp
         real(c_double), value :: G4ppx
         real(c_double), value :: G4pxx
         real(c_double), value :: G4xxx
         real(c_double), value :: G4pppp
         real(c_double), value :: G4pppx
         real(c_double), value :: G4ppxx
         real(c_double), value :: G4pxxx
         real(c_double), value :: G4xxxx
         real(c_double), value :: G5
         real(c_double), value :: G5p
         real(c_double), value :: G5x
         real(c_double), value :: G5pp
         real(c_double), value :: G5px
         real(c_double), value :: G5xx
         real(c_double), value :: G5ppp
         real(c_double), value :: G5ppx
         real(c_double), value :: G5pxx
         real(c_double), value :: G5xxx
         real(c_double), value :: G5pppp
         real(c_double), value :: G5pppx
         real(c_double), value :: G5ppxx
         real(c_double), value :: G5pxxx
         real(c_double), value :: G5xxxx

         real(c_double) :: spi1
         real(c_double) :: spi2
         real(c_double) :: spi3
         real(c_double) :: spih
         real(c_double) :: spie
         real(c_double) :: spim
         real(c_double) :: s00
         real(c_double) :: s00k
         real(c_double) :: s00p
         real(c_double) :: s0i
         real(c_double) :: s0ip
         real(c_double) :: sii
         real(c_double) :: siik
         real(c_double) :: siip
         real(c_double) :: siipp
         real(c_double) :: sij

      end subroutine horndeski_perturbeq_coefficients_lin


      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that computes the coefficients of perturbation equations to linear order in phidot
      subroutine horndeski_perturbeq_coefficients(a,H,dphi,rptot,drptot,G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx,spi1,spi2,spi3,spih,spie,spim,s00,s00k,s00p,s0i,s0ip,sii,siik,siip,siipp,sij) bind(C)
         use iso_c_binding, only: c_double
         implicit none

         real(c_double), value :: a      !<- the scale factor
         real(c_double), value :: H      !<- the physcial Hubble parameter
         real(c_double), value :: dphi   !<- conformal time derivative of the scalar field
         real(c_double), value :: rptot  !<- \f$ \rho_{matter} + p_{matter} \f$
         real(c_double), value :: drptot !<- conformal time derivative of \f$ \rho_{matter} + p_{matter} \f$
         real(c_double), value :: G2     !<- Horndeski Lagrangian function \f$ G2 \f$
         real(c_double), value :: G2p    !<- Horndeski Lagrangian function \f$ G2_\phi \f$
         real(c_double), value :: G2x    !<- Horndeski Lagrangian function \f$ G2_X \f$
         real(c_double), value :: G2pp
         real(c_double), value :: G2px
         real(c_double), value :: G2xx
         real(c_double), value :: G2ppp
         real(c_double), value :: G2ppx
         real(c_double), value :: G2pxx
         real(c_double), value :: G2xxx
         real(c_double), value :: G3p
         real(c_double), value :: G3x
         real(c_double), value :: G3pp
         real(c_double), value :: G3px
         real(c_double), value :: G3xx
         real(c_double), value :: G3ppp
         real(c_double), value :: G3ppx
         real(c_double), value :: G3pxx
         real(c_double), value :: G3xxx
         real(c_double), value :: G4
         real(c_double), value :: G4p
         real(c_double), value :: G4x
         real(c_double), value :: G4pp
         real(c_double), value :: G4px
         real(c_double), value :: G4xx
         real(c_double), value :: G4ppp
         real(c_double), value :: G4ppx
         real(c_double), value :: G4pxx
         real(c_double), value :: G4xxx
         real(c_double), value :: G4pppp
         real(c_double), value :: G4pppx
         real(c_double), value :: G4ppxx
         real(c_double), value :: G4pxxx
         real(c_double), value :: G4xxxx
         real(c_double), value :: G5
         real(c_double), value :: G5p
         real(c_double), value :: G5x
         real(c_double), value :: G5pp
         real(c_double), value :: G5px
         real(c_double), value :: G5xx
         real(c_double), value :: G5ppp
         real(c_double), value :: G5ppx
         real(c_double), value :: G5pxx
         real(c_double), value :: G5xxx
         real(c_double), value :: G5pppp
         real(c_double), value :: G5pppx
         real(c_double), value :: G5ppxx
         real(c_double), value :: G5pxxx
         real(c_double), value :: G5xxxx

         real(c_double) :: spi1
         real(c_double) :: spi2
         real(c_double) :: spi3
         real(c_double) :: spih
         real(c_double) :: spie
         real(c_double) :: spim
         real(c_double) :: s00
         real(c_double) :: s00k
         real(c_double) :: s00p
         real(c_double) :: s0i
         real(c_double) :: s0ip
         real(c_double) :: sii
         real(c_double) :: siik
         real(c_double) :: siip
         real(c_double) :: siipp
         real(c_double) :: sij

      end subroutine horndeski_perturbeq_coefficients

   end interface

   private

   public EFTCAMB_Hdsk

   !----------------------------------------------------------------------------------------
   !> This is the type that contains the definition of the full Horndeski model.
   type, extends ( EFTCAMB_full_model ) :: EFTCAMB_Hdsk

      ! model and parameterization:
      integer                             :: model           !< Model selection flag of covariant Horndeski models
      real(dl), allocatable, dimension(:) :: model_params    !< model parameters

      real(dl) :: maxfrac_hdsk = 0._dl   !< max of Horndeski field energy fraction
      real(dl) :: z_maxfrac_hdsk = 0._dl !< redshift where Horndeski field energy fraction peaks

      ! background IC:
      real(dl) :: Lambda_CC = 0._dl    !< The cosmological constant, non-zero in some EDE models
      real(dl) :: H_ini                !< initial value of Hubble
      real(dl) :: phi_ini              !< initial value of the scalar field
      real(dl) :: phidot_ini           !< initial conformal time derivative of the scalar field
      logical  :: model_specifc_ic     !< flag whether using the model specific initial condition, might not be defined for every model
      logical  :: evolve_hubble        !< flag whether ingrating H' instead of applying the energy constraint
      logical  :: evolve_delta_phi     !< flag whether computing the delta phi equation coefficients
      logical  :: evolve_metric_h      !< flag whether integrating metric perturbation h instead of using constraint

      ! shooting flag
      logical  :: is_shooting          !< flag whether only solve background without computing EFT functions

      ! In some functions we abuse notation and store the conformal time derivative (instead of x) in yp, ypp, etc., they will be marked with (*)

      ! the background solutions
      type( equispaced_linear_interpolate_function_1D ) :: tau            !< conformal time
      type( equispaced_linear_interpolate_function_1D ) :: H              !< physical Hubble
      type( equispaced_linear_interpolate_function_1D ) :: dH             !< conformal time derivative of Hubble
      type( equispaced_linear_interpolate_function_1D ) :: ddH            !< second conformal time derivative of Hubble
      type( equispaced_linear_interpolate_function_1D ) :: rhom           !< energy density of matter (everything except the Horndeski field)
      type( equispaced_linear_interpolate_function_1D ) :: presm          !< pressure of matter
      type( equispaced_linear_interpolate_function_1D ) :: rhomdot        !< conformal time derivative of energy density of matter
      type( equispaced_linear_interpolate_function_1D ) :: presmdot       !< conformal time derivative of pressure of matter
      type( equispaced_linear_interpolate_function_1D ) :: phi            !< scalar field value
      type( equispaced_linear_interpolate_function_1D ) :: dphi           !< conformal time derivative of phi
      type( equispaced_linear_interpolate_function_1D ) :: ddphi          !< second conformal time derivative of phi
      type( equispaced_linear_interpolate_function_1D ) :: dddphi         !< third conformal time derivative of phi
      type( equispaced_linear_interpolate_function_1D ) :: tau_hdsk       !< time scale of the horndeski field

      ! the EFT functions, TODO: implement analytic derivatives
      type( equispaced_linear_interpolate_function_1D ) :: EFTlambda      !< (*) The EFT function Lambda.
      type( equispaced_linear_interpolate_function_1D ) :: EFTc           !< (*) The EFT function c.
      type( equispaced_linear_interpolate_function_1D ) :: EFTomega       !< (*) The EFT function Omega.
      type( equispaced_linear_interpolate_function_1D ) :: EFTgamma1      !< (*) The EFT function Gamma1.
      type( equispaced_linear_interpolate_function_1D ) :: EFTgamma2      !< (*) The EFT function Gamma2.
      type( equispaced_linear_interpolate_function_1D ) :: EFTgamma3      !< (*) The EFT function Gamma3.
      ! the alternative EFT functions (computed analytically)
      type( equispaced_linear_interpolate_function_1D ) :: EFTM
      type( equispaced_linear_interpolate_function_1D ) :: EFTMdot
      type( equispaced_linear_interpolate_function_1D ) :: EFTMdotdot
      type( equispaced_linear_interpolate_function_1D ) :: EFTalphaM
      type( equispaced_linear_interpolate_function_1D ) :: EFTalphaMdot
      type( equispaced_linear_interpolate_function_1D ) :: EFTalphaB
      type( equispaced_linear_interpolate_function_1D ) :: EFTalphaBdot
      type( equispaced_linear_interpolate_function_1D ) :: EFTalphaK
      type( equispaced_linear_interpolate_function_1D ) :: EFTalphaKdot
      type( equispaced_linear_interpolate_function_1D ) :: EFTalphaT
      type( equispaced_linear_interpolate_function_1D ) :: EFTalphaTdot

      ! the interpolated delta_phi perturbation equation coefficients
      type( equispaced_linear_interpolate_function_1D ) :: spi1      !< coefficients of delta_phi equation: deltaphi'' + 2(1+spi1)aH deltaphi' + ((1 + spi2)k^2 + spi3*a^2H^2) deltaphi = spih h' + spie k^2 eta/aH + spim deltap_m
      type( equispaced_linear_interpolate_function_1D ) :: spi2
      type( equispaced_linear_interpolate_function_1D ) :: spi3
      type( equispaced_linear_interpolate_function_1D ) :: spih
      type( equispaced_linear_interpolate_function_1D ) :: spie
      type( equispaced_linear_interpolate_function_1D ) :: spim
      type( equispaced_linear_interpolate_function_1D ) :: s00       !< coefficients of delta_phi sources in the 00 equation: (s00k k^2 + s00) deltaphi + s00p deltaphi'
      type( equispaced_linear_interpolate_function_1D ) :: s00k
      type( equispaced_linear_interpolate_function_1D ) :: s00p
      type( equispaced_linear_interpolate_function_1D ) :: s0i       !< coefficients of delta_phi sources in the 0i equation: s0i deltaphi + s0ip deltaphi'
      type( equispaced_linear_interpolate_function_1D ) :: s0ip
      type( equispaced_linear_interpolate_function_1D ) :: sii       !< coefficients of delta_phi sources in the trace equation: (sii + siik k^2) delta_phi + siip delta_phi' + siipp delta_phi''
      type( equispaced_linear_interpolate_function_1D ) :: siik
      type( equispaced_linear_interpolate_function_1D ) :: siip
      type( equispaced_linear_interpolate_function_1D ) :: siipp
      type( equispaced_linear_interpolate_function_1D ) :: sij       !< coefficients of delta_phi sources in the non-diagnal equation: sij delta_phi


      ! some precision parameters:
      integer  :: interpolation_num_points = 6000                 !< Number of points sampled by the background solver code.
      real(dl) :: x_initial                = log(1d-14)           !< log(a start)
      real(dl) :: x_final                  = 0._dl                !< log(a final)

      real(dl) :: a_pertcutoff_before      = -1                   !< Turnoff scalar field perturbation before this time, minus means do not turn off
      real(dl) :: a_pertcutoff_after       = -1                   !< Turnoff scalar field perturbation after this time, minus means do not turn off
      real(dl) :: min_rel_timescale = 1d-5              !< Lowerbound of horndeski field timescale relative to Hubble time scale
      real(dl) :: hubble_friction = 5._dl               !< friction term to dynamically turn off scalar field perturbation

   contains

      ! initialization of the model:
      procedure :: read_model_selection            => EFTCAMBHdskReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
      procedure :: allocate_model_selection        => EFTCAMBHdskAllocateModelSelection      !< subroutine that allocates the model selection.
      procedure :: init_model_parameters           => EFTCAMBHdskInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
      procedure :: init_model_parameters_from_file => EFTCAMBHdskInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.

      ! utility functions:
      procedure :: compute_param_number  => EFTCAMBHdskComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
      procedure :: feedback              => EFTCAMBHdskFeedback                   !< subroutine that prints on the screen feedback information about the model.
      procedure :: parameter_names       => EFTCAMBHdskParameterNames             !< subroutine that returns the i-th parameter name of the model.
      procedure :: parameter_names_latex => EFTCAMBHdskParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
      procedure :: parameter_values      => EFTCAMBHdskParameterValues            !< subroutine that returns the i-th parameter value.

      ! CAMB related procedures:
      procedure :: compute_background_EFT_functions  => EFTCAMBHdskBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
      procedure :: compute_secondorder_EFT_functions => EFTCAMBHdskSecondOrderEFTFunctions  !< subroutine that computes the value of the second order EFT functions at a given time.
      procedure :: compute_dtauda                    => EFTCAMBHdskComputeDtauda            !< function that computes dtauda = 1/sqrt(a^2H^2).
      procedure :: compute_adotoa                    => EFTCAMBHdskComputeAdotoa            !< subroutine that computes adotoa = H.
      procedure :: compute_H_derivs                  => EFTCAMBHdskComputeHubbleDer         !< subroutine that computes the two derivatives wrt conformal time of H.

      ! background solver:
      procedure :: initialize_background           => EFTCAMBHdskInitBackground               !< subroutine that initializes the background of Hdsk.
      procedure :: solve_background_equations      => EFTCAMBHdskSolveBackgroundEquations     !< subroutine that solves the Hdsk background equations.
      procedure :: find_initial_conditions         => EFTCAMBHdskFindInitialConditions        !< subroutine that searches for the scalar field IC based on the specific model

      ! horndeski background
      procedure :: horndeski_lagrangian_coefficients => EFTCAMBHdskLagrangianCoefficients     !< subroutine that defines the covariant model

   end type

   ! ---------------------------------------------------------------------------------------------

   ! define debug files
   type(TTextFile) :: file_debug_1, file_debug_2, file_debug_3, file_debug_4

   ! ---------------------------------------------------------------------------------------------

contains

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that reads the parameters of the model from file.
   subroutine EFTCAMBHdskReadModelSelectionFromFile( self, Ini, eft_error )

      implicit none

      class(EFTCAMB_Hdsk) :: self      !< the base class
      type(TIniFile)    :: Ini       !< Input ini file
      integer           :: eft_error !< error code: 0 all fine, 1 initialization failed

      ! read model selection flags:
      self%model = Ini%Read_Int( 'Horndeski_model', 1 )
      ! read number of model parameters:
      self%parameter_number = Ini%Read_Int( 'Horndeski_parameter_number', 1 ) + 2
      ! read IC flag:
      self%model_specifc_ic = Ini%Read_Logical( 'Horndeski_model_specific_ic', .False. )
      ! read hubble evolution flag
      self%evolve_hubble = Ini%Read_Logical( 'Horndeski_evolve_hubble', .True. )
      ! read shooting flag
      self%is_shooting = Ini%Read_Logical( 'Horndeski_shooting', .False. )
      ! read perturbation equation flag:
      self%evolve_delta_phi = Ini%Read_Logical( 'EFT_evolve_delta_phi', .False. )
      self%evolve_metric_h  = Ini%Read_Logical( 'EFT_evolve_metric_h', .False. )
      ! read perturbation turn off flags
      self%a_pertcutoff_after = Ini%Read_Double( 'a_pertcutoff_after', -1._dl )
      self%a_pertcutoff_before = Ini%Read_Double( 'a_pertcutoff_before', -1._dl )
      self%min_rel_timescale = Ini%Read_Double( 'Horndeski_min_rel_timescale', 1d-5)
      self%hubble_friction = Ini%Read_Double( 'Horndeski_hubble_friction', 5._dl)
      ! read precision parameters
      self%interpolation_num_points = Ini%Read_Int('Horndeski_background_interpolation_num_points', 6000)
      self%x_initial = log(Ini%Read_Double('Horndeski_background_a_ini', 1d-14))

   end subroutine EFTCAMBHdskReadModelSelectionFromFile

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that allocates the model selection.
   subroutine EFTCAMBHdskAllocateModelSelection( self, Ini, eft_error )

      implicit none

      class(EFTCAMB_Hdsk) :: self      !< the base class
      type(TIniFile)    :: Ini       !< Input ini file
      integer           :: eft_error !< error code: 0 all fine, 1 initialization failed

      integer :: temp_feedback

      ! get feedback flag:
      temp_feedback = Ini%Read_Int('feedback_level', 0)

      select case (self%model)
       case(1)
         self%name = "Power-law quintessence"
       case(2)
         self%name = "Axion-like quintessence"
       case(3)
         self%name = "Horndeski EDE XboxPhi"
       case(4)
         self%name = "Horndeski EDE X2"
       case(5)
         self%name = "Horndeski EDE Phi2X"
       case(6)
         self%name = "Scaling cubic Galileon"
       case default
         self%name = "Covariant Horndeski Theory #"//integer_to_string(self%model)
      end select

   end subroutine EFTCAMBHdskAllocateModelSelection

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that initializes the model parameters based on the values found in an input array.
   subroutine EFTCAMBHdskInitModelParameters( self, array )

      implicit none

      class(EFTCAMB_Hdsk)                             :: self   !< the base class
      real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters of the model.

      integer                                         :: i

      ! read the field IC
      self%phi_ini = array(1)
      self%phidot_ini = array(2)

      ! read the model parameters:
      if ( allocated(self%model_params) ) deallocate(self%model_params)
      allocate( self%model_params(self%parameter_number-2) )
      do i = 1, self%parameter_number - 2
         self%model_params(i) = array(i+2)
      end do

   end subroutine EFTCAMBHdskInitModelParameters

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that reads the parameters of the model from file.
   subroutine EFTCAMBHdskInitModelParametersFromFile( self, Ini, eft_error )

      implicit none

      class(EFTCAMB_Hdsk)  :: self      !< the base class
      type(TIniFile)     :: Ini       !< Input ini file
      integer            :: eft_error !< error code: 0 all fine, 1 initialization failed

      integer :: i

      ! read the model parameters:
      if ( allocated(self%model_params) ) deallocate(self%model_params)
      allocate( self%model_params(self%parameter_number-2) )
      do i = 1, self%parameter_number-2
         self%model_params(i) = Ini%Read_Double( 'Horndeski_param'//integer_to_string(i), 0._dl )
      end do

      ! read the field IC:
      self%phi_ini = Ini%Read_Double( "Honrdeski_phi_ini", 0._dl )
      if (.not. self%model_specifc_ic) then
         self%phidot_ini = Ini%Read_Double( "Honrdeski_phidot_ini", 0._dl )
      end if

   end subroutine EFTCAMBHdskInitModelParametersFromFile

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that computes the number of parameters of the model.
   subroutine EFTCAMBHdskComputeParametersNumber( self )

      implicit none

      class(EFTCAMB_Hdsk)  :: self   !< the base class

   end subroutine EFTCAMBHdskComputeParametersNumber

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that prints on the screen feedback information about the model.
   subroutine EFTCAMBHdskFeedback( self, print_params )

      implicit none

      class(EFTCAMB_Hdsk)   :: self         !< the base class
      logical, optional     :: print_params   !< optional flag that decised whether to print numerical values of the parameters.

      integer :: i

      write(*,*)
      write(*,'(a,a)')    '   Model               =  ', self%name
      write(*,'(a,I3)')   '   Number of params    ='  , self%parameter_number

      write(*,*)
      write(*,'(a24,F12.6)') '   phi_ini          =', self%phi_ini
      if ( .not. self%model_specifc_ic ) then
         write(*,'(a24,F12.6)') '   phidot_ini          =', self%phidot_ini
      end if

      write(*,*)
      write(*,'(a26,F12.6)') 'Max DE energy fraction =  ', self%maxfrac_hdsk
      write(*,'(a26,F12.6)') 'Max DE energy fraction z =', self%z_maxfrac_hdsk

      write(*,*)
      do i = 1, self%parameter_number-2
         write(*,*) '   Horndeski parameter #', i, ' = ', self%model_params(i)
      end do

   end subroutine EFTCAMBHdskFeedback

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that returns the i-th parameter name of the model
   subroutine EFTCAMBHdskParameterNames( self, i, name )

      implicit none

      class(EFTCAMB_Hdsk)         :: self   !< the base class
      integer     , intent(in)    :: i      !< the index of the parameter
      character(*), intent(out)   :: name   !< the output name of the i-th parameter

      ! check validity of input:
      if ( i<=0 .or. i>self%parameter_number ) then
         write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
         write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
         call MpiStop('EFTCAMB error')
      else if ( i==1 ) then
         name = TRIM('phi_ini')
         return
      else if ( i==2 ) then
         name = TRIM('phidot_ini')
         return
      else
         name = TRIM('Horndeski_param')//integer_to_string(i-2)
         return
      end if

   end subroutine EFTCAMBHdskParameterNames

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that returns the i-th parameter name of the model
   subroutine EFTCAMBHdskParameterNamesLatex( self, i, latexname )

      implicit none

      class(EFTCAMB_Hdsk)         :: self       !< the base class
      integer     , intent(in)    :: i          !< The index of the parameter
      character(*), intent(out)   :: latexname  !< the output latex name of the i-th parameter

      ! check validity of input:
      if ( i<=0 .or. i>self%parameter_number ) then
         write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
         write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
         call MpiStop('EFTCAMB error')
      else if ( i==1 ) then
         latexname = TRIM('\phi_\mathrm{ini}')
         return
      else if ( i==2 ) then
         latexname = TRIM('\dot{\phi}_\mathrm{ini}')
         return
      else
         latexname = TRIM('Horndeski param')//integer_to_string(i-2)
         return
      end if

   end subroutine EFTCAMBHdskParameterNamesLatex

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that returns the i-th parameter name of the model
   subroutine EFTCAMBHdskParameterValues( self, i, value )

      implicit none

      class(EFTCAMB_Hdsk)         :: self   !< the base class
      integer , intent(in)        :: i      !< The index of the parameter
      real(dl), intent(out)       :: value  !< the output value of the i-th parameter

      ! check validity of input:
      if ( i<=0 .or. i>self%parameter_number ) then
         write(*,'(a,I3)') 'EFTCAMB error: no parameter corresponding to number ', i
         write(*,'(a,I3)') 'Total number of parameters is ', self%parameter_number
         call MpiStop('EFTCAMB error')
      else if ( i==1 ) then
         value = self%phi_ini
         return
      else if ( i==2 ) then
         value = self%phidot_ini
         return
      else
         value = self%model_params(i-2)
         return
      end if

   end subroutine EFTCAMBHdskParameterValues

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that computes the value of the background EFT functions at a given time.
   subroutine EFTCAMBHdskBackgroundEFTFunctions( self, a, eft_par_cache, eft_cache )

      implicit none

      class(EFTCAMB_Hdsk)                           :: self          !< the base class
      real(dl), intent(in)                          :: a             !< the input scale factor
      type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
      type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

      real(dl) :: x, mu
      integer  :: ind

      ! protect against calling zero and convert to log:
      if ( a < exp( self%H%x(1) ) .or. self%is_shooting ) then
         eft_cache%EFTc         = 0._dl
         eft_cache%EFTLambda    = 0._dl
         eft_cache%EFTcdot      = 0._dl
         eft_cache%EFTLambdadot = 0._dl
         eft_cache%EFTOmegaV    = 0._dl
         eft_cache%EFTOmegaP    = 0._dl
         eft_cache%EFTOmegaPP   = 0._dl
         eft_cache%EFTOmegaPPP  = 0._dl
         eft_cache%phi_scf      = 0._dl
         eft_cache%dphi         = 0._dl
         eft_cache%ddphi        = 0._dl
         eft_cache%dddphi       = 0._dl
         eft_cache%dtau_hdsk    = 1._dl / eft_par_cache%h0_Mpc
      else
         x = log(a)
         ! All the EFT functions are sampled on the same time grid.
         call self%H%precompute(x, ind, mu )
         ! EFT functions, remeber we simply store the a (or conformal time for Lambda and c) derivative of them in yp, ypp. They are NOT wrt x!
         eft_cache%EFTc         = self%EFTc%value( x, index=ind, coeff=mu )
         eft_cache%EFTLambda    = self%EFTlambda%value( x, index=ind, coeff=mu )
         eft_cache%EFTcdot      = self%EFTc%first_derivative( x, index=ind, coeff=mu )
         eft_cache%EFTLambdadot = self%EFTlambda%first_derivative( x, index=ind, coeff=mu )
         eft_cache%EFTOmegaV    = self%EFTomega%value( x, index=ind, coeff=mu )
         eft_cache%EFTOmegaP    = self%EFTomega%first_derivative( x, index=ind, coeff=mu )
         eft_cache%EFTOmegaPP   = self%EFTomega%second_derivative( x, index=ind, coeff=mu )
         eft_cache%EFTOmegaPPP  = self%EFTomega%third_derivative( x, index=ind, coeff=mu )
         eft_cache%phi_scf      = self%phi%value( x, index=ind, coeff=mu )
         eft_cache%dphi         = self%dphi%value( x, index=ind, coeff=mu )
         eft_cache%ddphi        = self%ddphi%value( x, index=ind, coeff=mu )
         eft_cache%dddphi       = self%dddphi%value( x, index=ind, coeff=mu )
         eft_cache%dtau_hdsk    = self%tau_hdsk%value( x, index=ind, coeff=mu )
      end if

      ! if ((abs(a-1._dl) < 0.001) .or. (abs(a*10._dl-1._dl) < 0.01) .or. (abs(a*1.d2-1._dl) < 0.01) .or. (abs(a*1.d3-1._dl) < 0.01) .or. (abs(a*1.d4-1._dl) < 0.01) .or. (abs(a*1.d5-1._dl) < 0.01)) then
      !    call self%compute_adotoa(a, eft_par_cache, eft_cache)
      !    call self%compute_H_derivs(a, eft_par_cache, eft_cache)
      !    write(*,*) a, eft_cache%adotoa, eft_cache%Hdot, eft_cache%Hdotdot, eft_cache%Hdotdotdot, self%compute_dtauda(a, eft_par_cache, eft_cache), eft_cache%EFTc, eft_cache%EFTLambda, eft_cache%EFTcdot, eft_cache%EFTLambdadot 
      ! end if

   end subroutine EFTCAMBHdskBackgroundEFTFunctions

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that computes the value of the background EFT functions at a given time.
   subroutine EFTCAMBHdskSecondOrderEFTFunctions( self, a, eft_par_cache, eft_cache )

      implicit none

      class(EFTCAMB_Hdsk)                           :: self          !< the base class
      real(dl), intent(in)                          :: a             !< the input scale factor
      type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
      type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

      real(dl) :: x, mu
      integer  :: ind

      x = log(a)
      
      ! protect against calling zero and convert to log:
      if ( x < self%H%x(1) .or. self%is_shooting ) then
         
         ! alpha EFT functions
         eft_cache%alphaB       = 0._dl
         eft_cache%alphaBdot    = 0._dl
         eft_cache%alphaT       = 0._dl
         eft_cache%alphaTdot    = 0._dl
         eft_cache%alphaK       = 0._dl
         eft_cache%alphaKdot    = 0._dl
         eft_cache%alphaM       = 0._dl
         eft_cache%alphaMdot    = 0._dl
         eft_cache%Meff2        = 1._dl

         ! gamma EFT functions
         eft_cache%EFTGamma1V  = 0._dl
         eft_cache%EFTGamma1P  = 0._dl
         eft_cache%EFTGamma2V  = 0._dl
         eft_cache%EFTGamma2P  = 0._dl
         eft_cache%EFTGamma3V  = 0._dl
         eft_cache%EFTGamma3P  = 0._dl
         eft_cache%EFTGamma4V  = 0._dl
         eft_cache%EFTGamma4P  = 0._dl
         eft_cache%EFTGamma4PP = 0._dl
         eft_cache%EFTGamma5V  = 0._dl
         eft_cache%EFTGamma5P  = 0._dl
         eft_cache%EFTGamma6V  = 0._dl
         eft_cache%EFTGamma6P  = 0._dl

      ! deltaphi functions
         eft_cache%spi1  = self%hubble_friction
         eft_cache%spi2  = 0._dl
         eft_cache%spi3  = 0._dl
         eft_cache%spih  = 0._dl
         eft_cache%spie  = 0._dl
         eft_cache%spim  = 0._dl
         eft_cache%s00   = 0._dl
         eft_cache%s00k  = 0._dl
         eft_cache%s00p  = 0._dl
         eft_cache%s0i   = 0._dl
         eft_cache%s0ip  = 0._dl
         eft_cache%sii   = 0._dl
         eft_cache%siik  = 0._dl
         eft_cache%siip  = 0._dl
         eft_cache%siipp = 0._dl
         eft_cache%sij   = 0._dl

      else
         
         ! All the EFT functions are sampled on the same time grid.
         call self%H%precompute(x, ind, mu )

         ! alpha EFT functions
         eft_cache%alphaB       = self%EFTalphaB%value( x, index=ind, coeff=mu )
         eft_cache%alphaBdot    = self%EFTalphaBdot%value( x, index=ind, coeff=mu )
         eft_cache%alphaT       = self%EFTalphaT%value( x, index=ind, coeff=mu )
         eft_cache%alphaTdot    = self%EFTalphaTdot%value( x, index=ind, coeff=mu )
         eft_cache%alphaK       = self%EFTalphaK%value( x, index=ind, coeff=mu )
         eft_cache%alphaKdot    = self%EFTalphaKdot%value( x, index=ind, coeff=mu )
         eft_cache%alphaM       = self%EFTalphaM%value( x, index=ind, coeff=mu )
         eft_cache%alphaMdot    = self%EFTalphaMdot%value( x, index=ind, coeff=mu )
         eft_cache%Meff2        = self%EFTM%value( x, index=ind, coeff=mu )

         ! gamma EFT functions, remeber we simply store the conformal time derivative of them in yp, ypp. They are NOT wrt x!
         eft_cache%EFTGamma1V  = self%EFTgamma1%value( x, index=ind, coeff=mu )
         eft_cache%EFTGamma1P  = self%EFTgamma1%first_derivative( x, index=ind, coeff=mu )
         eft_cache%EFTGamma2V  = self%EFTgamma2%value( x, index=ind, coeff=mu )
         eft_cache%EFTGamma2P  = self%EFTgamma2%first_derivative( x, index=ind, coeff=mu )
         eft_cache%EFTGamma3V  = self%EFTgamma3%value( x, index=ind, coeff=mu )
         eft_cache%EFTGamma3P  = self%EFTgamma3%first_derivative( x, index=ind, coeff=mu )
         eft_cache%EFTGamma4V  = -eft_cache%EFTGamma3V
         eft_cache%EFTGamma4P  = -eft_cache%EFTGamma3P
         eft_cache%EFTGamma4PP = -self%EFTgamma3%second_derivative( x, index=ind, coeff=mu )
         eft_cache%EFTGamma5V  = +0.5_dl*eft_cache%EFTGamma3V
         eft_cache%EFTGamma5P  = +0.5_dl*eft_cache%EFTGamma3P
         eft_cache%EFTGamma6V  = 0._dl
         eft_cache%EFTGamma6P  = 0._dl

         ! deltaphi functions
         eft_cache%spi1  = self%spi1%value( x, index=ind, coeff=mu )
         eft_cache%spi2  = self%spi2%value( x, index=ind, coeff=mu )
         eft_cache%spi3  = self%spi3%value( x, index=ind, coeff=mu )
         eft_cache%spih  = self%spih%value( x, index=ind, coeff=mu )
         eft_cache%spie  = self%spie%value( x, index=ind, coeff=mu )
         eft_cache%spim  = self%spim%value( x, index=ind, coeff=mu )
         eft_cache%s00   = self%s00%value( x, index=ind, coeff=mu )
         eft_cache%s00k  = self%s00k%value( x, index=ind, coeff=mu )
         eft_cache%s00p  = self%s00p%value( x, index=ind, coeff=mu )
         eft_cache%s0i   = self%s0i%value( x, index=ind, coeff=mu )
         eft_cache%s0ip  = self%s0ip%value( x, index=ind, coeff=mu )
         eft_cache%sii   = self%sii%value( x, index=ind, coeff=mu )
         eft_cache%siik  = self%siik%value( x, index=ind, coeff=mu )
         eft_cache%siip  = self%siip%value( x, index=ind, coeff=mu )
         eft_cache%siipp = self%siipp%value( x, index=ind, coeff=mu )
         eft_cache%sij   = self%sij%value( x, index=ind, coeff=mu )
         ! convert x derivative to conformal time derivative
         eft_cache%sijdot= a*self%H%value( x, index=ind, coeff=mu )*self%sij%first_derivative( x, index=ind, coeff=mu )

      end if

      ! if ((abs(a-1._dl) < 0.001) .or. (abs(a*10._dl-1._dl) < 0.01) .or. (abs(a*1.d2-1._dl) < 0.01) .or. (abs(a*1.d3-1._dl) < 0.01) .or. (abs(a*1.d4-1._dl) < 0.01) .or. (abs(a*1.d5-1._dl) < 0.01)) then
      !    write(*,*) a, eft_cache%EFTGamma1V, eft_cache%EFTGamma1P, eft_cache%EFTGamma2V, eft_cache%EFTGamma2P
      ! end if

   end subroutine EFTCAMBHdskSecondOrderEFTFunctions

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that initializes the background of the quintessence model.
   subroutine EFTCAMBHdskInitBackground( self, params_cache, feedback_level, success, outroot )

      implicit none

      class(EFTCAMB_Hdsk)                           :: self           !< the base class
      type(TEFTCAMB_parameter_cache), intent(inout) :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters
      integer                       , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
      logical                       , intent(out)   :: success        !< wether the background initialization succeded or not
      character(LEN=*), optional    , intent(in)    :: outroot        !< the output root for the debug files

      integer  :: i, j, error_flag
      real(dl) :: phi, dphi, last_dphi, last_dphi_tau, tau, tau_hdsk, x, a, H, adotoa, Hdot, Hdotdot, cHdot, cHdotdot, rhov, presv, rhovdot, presvdot, rhom, presm, rhomdot, presmdot, tmp
      logical, dimension(:), allocatable :: delta_phi_lin_flag
      real(dl) :: RPH_PM_V, RPH_AT_V, RPH_PM_P, RPH_AT_P, RPH_PM_PP, RPH_AT_PP, RPH_PM_PPP, RPH_AT_PPP, RPH_AK_V, RPH_AB_V, RPH_AK_P, RPH_AB_P
      real(dl) :: daM, ddaM, dddaM, daT, ddaT, dddaT
      real(dl) :: spi1,spi2,spi3,spih,spie,spim,s00,s00k,s00p,s0i,s0ip,sii,siik,siip,siipp,sij
      real(dl) :: G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx
      real(dl) :: last_alphaB

      ! some feedback:
      if ( feedback_level>1 ) then
         write(*,'(a)') "***************************************************************"
         write(*,'(a)') ' EFTCAMB Horndeski background solver'
         write(*,'(a)')
      end if

      if ( DebugEFTCAMB .or. feedback_level>2 ) then
         call params_cache%print()
      end if

      ! initialize interpolating functions:
      call self%tau%initialize    ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%H%initialize      ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%dH%initialize      ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%ddH%initialize      ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%rhom%initialize      ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%presm%initialize      ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%rhomdot%initialize      ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%presmdot%initialize      ( self%interpolation_num_points, self%x_initial, self%x_final )

      call self%phi%initialize      ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%dphi%initialize      ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%ddphi%initialize      ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%dddphi%initialize      ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%tau_hdsk%initialize      ( self%interpolation_num_points, self%x_initial, self%x_final )

      call self%EFTc%initialize      ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%EFTlambda%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%EFTomega%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%EFTgamma1%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%EFTgamma2%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%EFTgamma3%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )

      call self%EFTM%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%EFTMdot%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%EFTMdotdot%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%EFTalphaB%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%EFTalphaBdot%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%EFTalphaK%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%EFTalphaKdot%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%EFTalphaM%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%EFTalphaMdot%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%EFTalphaT%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%EFTalphaTdot%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )

      call self%spi1%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%spi2%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%spi3%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%spih%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%spie%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%spim%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%s00%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%s00k%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%s00p%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%s0i%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%s0ip%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%sii%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%siik%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%siip%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%siipp%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )
      call self%sij%initialize ( self%interpolation_num_points, self%x_initial, self%x_final )

      ! find model specific ICs:
      call self%find_initial_conditions( params_cache, feedback_level, success )
      if ( success ) then
         if ( feedback_level>1 ) write(*,'(a,E13.4)') '   initial condition phidot_ini  = ', self%phidot_ini
      else if ( .not. success ) then
         return
      end if

      ! solve the background equations and store the solution:
      call self%solve_background_equations( params_cache, only_solve=self%is_shooting, success=success )
      if ( .not. success ) then
         if ( feedback_level>1 ) write(*,*) 'Background solve failed.'
         return
      end if

      params_cache%maxfrac_hdsk = self%maxfrac_hdsk
      params_cache%z_maxfrac_hdsk = self%z_maxfrac_hdsk

      ! stop here if no further output are required
      if ( self%is_shooting ) return
      
      ! first loop to initialize delta_phi_lin_flag and check alpha_B
      last_alphaB = self%EFTalphaB%y(1)
      j = 0
      if ( allocated(delta_phi_lin_flag) ) deallocate( delta_phi_lin_flag )
      allocate( delta_phi_lin_flag(self%H%num_points) )
      do i=1, self%H%num_points
         delta_phi_lin_flag(i) = .False.
         if ( (last_alphaB-2._dl)*(self%EFTalphaB%y(i)-2._dl) <= 0._dl ) j = j + 1
         last_alphaB = self%EFTalphaB%y(i)
      end do
      if ( j>0 .and. (.not. self%evolve_metric_h) ) then
         write(*,*) "alphaB crossing 2 detected, at which h'' constraint is numerically unstable. Consider setting EFT_evolve_metric_h = T."
      end if

      ! second loop to find dphi sign switching time and potential oscillation time scale
      last_dphi_tau = 0._dl
      last_dphi = 0._dl
      tau_hdsk = self%tau%y(self%tau%num_points) ! initialize to conformal time today
      do i=1, self%H%num_points
         a = Exp(self%H%x(i))
         tau = self%tau%y(i)
         H = self%H%y(i)
         dphi = self%dphi%y(i)
         ! only update tau_hdsk when scf perturbation is on and away from initial point
         if ( (a<self%a_pertcutoff_after .or. self%a_pertcutoff_after<0) .and. a>self%a_pertcutoff_before .and. tau>0.1_dl ) then
            if (dphi * last_dphi <= 0._dl) then
               ! we rely on a_pertcutoff_before to avoid relaxation period to the true background solution
               tau_hdsk = tau - last_dphi_tau
               last_dphi_tau = tau
               ! if a sign change detected, we force use of linear coefficients a few steps around it
               delta_phi_lin_flag(i) = .True.
               do j = 1, 5
                  delta_phi_lin_flag(MAX(i - j, 1)) = .True.
                  delta_phi_lin_flag(MIN(i + j, self%H%num_points)) = .True.
               end do
               if (feedback_level > 2) then
                  write(*,*) "phidot changing sign detected at z=", 1._dl / a - 1, "tau=", tau
                  if ( .not. self%evolve_delta_phi ) then
                     write(*,*) "pi equations are singular when phidot crosses 0. Consider setting EFT_evolve_delta_phi = T for this model."
                  end if
               end if
            end if
            ! when scf perturbation off, set to some safe big value
         else
            tau_hdsk = self%tau%y(self%tau%num_points)
         end if
         last_dphi = dphi
         self%tau_hdsk%y(i) = MAX(abs(tau_hdsk), abs(1._dl/a/H*self%min_rel_timescale))
      end do

      ! compute numeric derivatives (wrt x=loga) that might be needed later
      call self%EFTMdotdot%initialize_derivatives()
      call self%EFTalphaTdot%initialize_derivatives()
      call self%ddH%initialize_derivatives()

      ! third loop to compute EFT gammas and the delta_phi coefficients
      ! TODO: also check positivity here
      do i=1, self%H%num_points
         a = Exp(self%H%x(i))
         tau = self%tau%y(i)
         H = self%H%y(i)         ! physical Hubble
         Hdot = self%dH%y(i)     ! conformal time derivative of physical Hubble
         Hdotdot = self%ddH%y(i) ! second conformal time derivative of physical Hubble
         phi = self%phi%y(i)
         dphi = self%dphi%y(i)
         adotoa = a*H                   ! conformal Hubble
         cHdot =  adotoa**2 + a*Hdot    ! conformal time derivative of conformal Hubble
         cHdotdot = 2._dl*adotoa*cHdot + adotoa*a*Hdot + a*Hdotdot ! second conformal time derivative of conformal Hubble
         rhom = self%rhom%y(i)
         presm = self%presm%y(i)
         rhomdot = self%rhomdot%y(i)
         presmdot = self%presmdot%y(i)
         rhov = 3._dl*H**2 - rhom                         ! rho_DE = 3H^2 - rho_m
         presv = -(2._dl*Hdot/a + rhom + presm + rhov)    ! p_DE = -(2H'/a + rho_m + p_m + rho_DE)
         rhovdot = -3._dl*a*H*(rhov + presv)               ! rho_DE' = -3aH(rho_DE + p_DE)
         presvdot = -(2._dl*Hdotdot/a - 2._dl*H*Hdot + rhomdot + presmdot + rhovdot) ! p_DE' = -(2H''/a - 2HH' + rho_m' + p_m' + rho_DE')

         ! 1) obtain the required alpha's
         ! store some conformal time derivative for convenience
         daM = self%EFTMdot%y(i)
         ddaM = self%EFTMdotdot%y(i)
         dddaM = adotoa*self%EFTMdotdot%yp(i)
         daT = self%EFTalphaTdot%y(i)
         ddaT = adotoa*self%EFTalphaTdot%yp(i)
         dddaT = adotoa**2*self%EFTalphaTdot%ypp(i) + cHdot*self%EFTalphaTdot%yp(i)

         RPH_PM_V               = self%EFTM%y(i) - 1._dl
         RPH_PM_P               = daM/(a*adotoa)
         RPH_PM_PP              = ddaM/(a*adotoa)**2 - cHdot*daM/(a**2*adotoa**3) - daM/(a**2*adotoa)
         RPH_PM_PPP             = dddaM/(a*adotoa)**3 - 3._dl*cHdot*ddaM/(a**3*adotoa**4) - 3._dl*ddaM/(a**3*adotoa**2) - cHdotdot*daM/(a**3*adotoa**4) + 3._dl*cHdot**2*daM/(a**3*adotoa**5) + 3._dl*cHdot*daM/(a*adotoa)**3 + 2._dl*daM/(a**3*adotoa)
         RPH_AT_V               = self%EFTalphaT%y(i)
         RPH_AT_P               = daT/(a*adotoa)
         RPH_AT_PP              = ddaT/(a*adotoa)**2 - cHdot*daT/(a**2*adotoa**3) - daT/(a**2*adotoa)
         RPH_AT_PPP             = dddaT/(a*adotoa)**3 - 3._dl*cHdot*ddaT/(a**3*adotoa**4) - 3._dl*ddaT/(a**3*adotoa**2) - cHdotdot*daT/(a**3*adotoa**4) + 3._dl*cHdot**2*daT/(a**3*adotoa**5) + 3._dl*cHdot*daT/(a*adotoa)**3 + 2._dl*daT/(a**3*adotoa)
         RPH_AK_V               = self%EFTalphaK%y(i)
         RPH_AB_V               = self%EFTalphaB%y(i)
         RPH_AK_P               = self%EFTalphaKdot%y(i)/(a*adotoa)
         RPH_AB_P               = self%EFTalphaBdot%y(i)/(a*adotoa)
         
         ! 1.5) A convention difference between the eftcamb alpha_B and 1404.3713? A strange issue or not?
         RPH_AB_V = -0.5_dl * RPH_AB_V
         RPH_AB_P = -0.5_dl * RPH_AB_P

         ! 2) compute the EFT functions:
         ! NOTE: we abuse notation and store directly the a (or conformal time for Lambda and c) derivative to yp, ypp, yppp. NOT wrt x!
         self%EFTomega%y(i)     = RPH_PM_V +RPH_AT_V*(1._dl +RPH_PM_V)
         self%EFTomega%yp(i)    = (RPH_PM_V+1.0_dl)*RPH_AT_P +(1._dl +RPH_AT_V)*RPH_PM_P
         self%EFTomega%ypp(i)   = (RPH_PM_V+1.0_dl)*RPH_AT_PP +2._dl*RPH_PM_P*RPH_AT_P +(1._dl +RPH_AT_V)*RPH_PM_PP
         self%EFTomega%yppp(i)  = (RPH_PM_V+1.0_dl)*RPH_AT_PPP +3._dl*RPH_PM_PP*RPH_AT_P +3._dl*RPH_PM_P*RPH_AT_PP &
         & +(1._dl +RPH_AT_V)*RPH_PM_PPP
         self%EFTc%y(i)         = ( adotoa**2 - cHdot )*( self%EFTomega%y(i) + 0.5_dl*a*self%EFTomega%yp(i) ) &
         & -0.5_dl*( a*adotoa )**2*self%EFTomega%ypp(i)&
         & +0.5_dl*a*a*(rhov + presv)
         self%EFTlambda%y(i)    = -self%EFTomega%y(i)*( 2._dl*cHdot +adotoa**2 ) &
         & -a*self%EFTomega%yp(i)*( 2._dl*adotoa**2 + cHdot ) &
         & -( a*adotoa )**2*self%EFTomega%ypp(i) &
         & +a*a*presv
         self%EFTc%yp(i)        = +0.5_dl*a*a*(rhovdot + presvdot) &
         & -self%EFTomega%y(i)*( cHdotdot -4._dl*adotoa*cHdot +2._dl*adotoa**3 ) &
         & +0.5_dl*a*self%EFTomega%yp(i)*( -cHdotdot +adotoa*cHdot +adotoa**3) &
         & +0.5_dl*a**2*adotoa*self%EFTomega%ypp(i)*( adotoa**2 -3._dl*cHdot ) &
         & -0.5_dl*(a*adotoa)**3*self%EFTomega%yppp(i)
         self%EFTlambda%yp(i)   = -2._dl*self%EFTomega%y(i)*( cHdotdot -adotoa*cHdot -adotoa**3 ) &
         & -a*self%EFTomega%yp(i)*( +cHdotdot +5._dl*adotoa*cHdot -adotoa**3  ) &
         & -a**2*self%EFTomega%ypp(i)*adotoa*( +2._dl*adotoa**2 +3._dl*cHdot )&
         & -(a*adotoa)**3*self%EFTomega%yppp(i) &
         & +a*a*presvdot
         self%EFTgamma1%y(i)  = 0.25_dl*( RPH_AK_V*(1._dl +RPH_PM_V)*adotoa**2 &
         & -2._dl*self%EFTc%y(i) )/(params_cache%h0_Mpc**2*a**2)
         self%EFTgamma1%yp(i)  = - 0.5_dl*( RPH_AK_V*(1._dl +RPH_PM_V)*adotoa**2 &
         & -2._dl*self%EFTc%y(i) )/(params_cache%h0_Mpc**2*a**3) &
         & +0.25_dl*( RPH_AK_P*(1._dl +RPH_PM_V)*adotoa**2 &
         & +RPH_AK_V*RPH_PM_P*adotoa**2 &
         & +2._dl*RPH_AK_V*(1._dl +RPH_PM_V)*cHdot/a &
         & -4._dl*self%EFTc%y(i)/a -2._dl*self%EFTc%yp(i)/a/adotoa )/(params_cache%h0_Mpc**2*a**2)
         self%EFTgamma2%y(i)  = ( +2._dl*RPH_AB_V*(1._dl +RPH_PM_V) &
         & -a*self%EFTomega%yp(i) )*adotoa/(params_cache%h0_Mpc*a)
         self%EFTgamma2%yp(i)  = -( +2._dl*RPH_AB_V*(1._dl +RPH_PM_V) &
         & -a*self%EFTomega%yp(i) )*adotoa/(params_cache%h0_Mpc*a**2) &
         & -( -2._dl*(1._dl +RPH_PM_V)*( RPH_AB_P*adotoa**2 &
         & + RPH_AB_V*cHdot/a) &
         & - 2._dl*RPH_AB_V*adotoa**2*RPH_PM_P &
         & + self%EFTomega%yp(i)*( adotoa**2 +cHdot ) &
         & + a*adotoa**2*self%EFTomega%ypp(i) )/(params_cache%h0_Mpc*a*adotoa)
         self%EFTgamma3%y(i)  = -RPH_AT_V*(1._dl +RPH_PM_V)
         self%EFTgamma3%yp(i)  = -RPH_PM_P*RPH_AT_V -(1._dl +RPH_PM_V)*RPH_AT_P
         self%EFTgamma3%ypp(i) = -(1._dl +RPH_PM_V)*RPH_AT_PP - RPH_PM_PP*RPH_AT_V - 2._dl*RPH_PM_P*RPH_AT_P

         ! write(*,*) a, self%EFTlambda%y(i), self%EFTlambda%yp(i), self%EFTc%y(i), self%EFTc%yp(i)
         ! write(*,*) a, self%EFTgamma1%y(i), self%EFTgamma1%yp(i), self%EFTgamma2%y(i), self%EFTgamma2%yp(i), self%EFTgamma3%y(i), self%EFTgamma3%yp(i), self%EFTgamma3%ypp(i), 4._dl*self%model_params(1)*dphi/a/params_cache%h0_Mpc

         ! 3) compute the delta_phi perturbation coefficients
         if ( self%evolve_delta_phi .and. (a<self%a_pertcutoff_after .or. self%a_pertcutoff_after<0) .and. a>self%a_pertcutoff_before ) then
            ! we could also check for positivity here
            call self%horndeski_lagrangian_coefficients( a, phi, 0.5_dl*(dphi/a)**2, G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx)
            if (delta_phi_lin_flag(i)) then
               call horndeski_perturbeq_coefficients_lin(a, H, dphi, rhom+presm, rhomdot+presmdot, G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx, spi1, spi2, spi3, spih, spie, spim, s00, s00k, s00p, s0i, s0ip, sii, siik, siip, siipp, sij)
            else
               call horndeski_perturbeq_coefficients(a, H, dphi, rhom+presm, rhomdot+presmdot, G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx, spi1, spi2, spi3, spih, spie, spim, s00, s00k, s00p, s0i, s0ip, sii, siik, siip, siipp, sij)
            end if
         else
            spi1 = self%hubble_friction
            spi2 = 0._dl
            spi3 = 0._dl
            spih = 0._dl
            spie = 0._dl
            spim = 0._dl
            s00 = 0._dl
            s00k = 0._dl
            s00p = 0._dl
            s0i = 0._dl
            s0ip = 0._dl
            sii = 0._dl
            siik = 0._dl
            siip = 0._dl
            siipp = 0._dl
            sij = 0._dl
         end if
         self%spi1%y(i) = spi1
         self%spi2%y(i) = spi2
         self%spi3%y(i) = spi3
         self%spih%y(i) = spih
         self%spie%y(i) = spie
         self%spim%y(i) = spim
         self%s00%y(i) = s00
         self%s00k%y(i) = s00k
         self%s00p%y(i) = s00p
         self%s0i%y(i) = s0i
         self%s0ip%y(i) = s0ip
         self%sii%y(i) = sii
         self%siik%y(i) = siik
         self%siip%y(i) = siip
         self%siipp%y(i) = siipp
         self%sij%y(i) = sij

      end do

      ! compute numeric derivative of sij, needed for newtonian sources
      if ( self%evolve_delta_phi ) then
         call self%sij%initialize_derivatives()
      end if

      ! debug output
      ! file_debug_1%unit = 32
      ! call file_debug_1%CreateFile( TRIM(outroot)//'background_hdsk.dat' )
      ! write(file_debug_1%unit, *) '# z   tau   H   Hdot    rhom    presm    rhomdot    presmdot   phi   phidot   phidotdot   alphaK   alphaB   spi1   spi2   spi3   spih   spie   spim   s00   s00k   s00p   s0i   s0ip   sii   siik   siip   siipp   sij'
      ! do i = 1, self%H%num_points
      !    write(file_debug_1%unit, *) exp(-self%H%x(i))-1._dl, self%tau%y(i), self%H%y(i), self%dH%y(i), self%rhom%y(i), self%presm%y(i), self%rhomdot%y(i), self%presmdot%y(i), self%phi%y(i), self%dphi%y(i), self%ddphi%y(i), self%EFTalphaK%y(i), self%EFTalphaB%y(i), self%spi1%y(i), self%spi2%y(i), self%spi3%y(i), self%spih%y(i), self%spie%y(i), self%spim%y(i), self%s00%y(i), self%s00k%y(i), self%s00p%y(i), self%s0i%y(i), self%s0ip%y(i), self%sii%y(i), self%siik%y(i), self%siip%y(i), self%siipp%y(i), self%sij%y(i)
      ! end do
      ! call file_debug_1%close()

      deallocate( delta_phi_lin_flag )

   end subroutine EFTCAMBHdskInitBackground

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that solves the Horndeski background equations.
   subroutine EFTCAMBHdskSolveBackgroundEquations( self, params_cache, only_solve, success )

      implicit none

      class(EFTCAMB_Hdsk)                          :: self          !< the base class.
      type(TEFTCAMB_parameter_cache), intent(in)   :: params_cache  !< a EFTCAMB parameter cache containing cosmological parameters.
      logical , optional                           :: only_solve    !< logical flag that tells the code wether to compute only the solution or also the EFT functions.
      logical , intent(out)                        :: success       !< whether the calculation ended correctly or not

      integer :: num_eq = 3   !<  Number of equations

      ! y = tau, phi, phi', (H)
      real(dl), allocatable :: y(:), ydot(:)

      ! odepack quantities:
      integer  :: itol, itask, istate, iopt, LRN, LRS, LRW, LIS, LIN, LIW, JacobianMode, i
      real(dl) :: rtol, atol, t1, t2
      real(dl), allocatable :: rwork(:)
      integer,  allocatable :: iwork(:)

      ! background quantities that are shared among the derivs and output routines:
      real(dl) :: a, a2, grhob_t, grhoc_t, grhor_t, grhog_t, grhonu_tot, gpinu_tot, grhonudot_tot, gpinudot_tot
      real(dl) :: grhonu, gpinu, grhonudot, gpinudot, grhormass_t, grho_matter, gpres_matter
      real(dl) :: phi, phi_prime, phi_prime_prime, phi_prime_prime_prime, H, H_prime, H_prime_prime, adotoa
      real(dl) :: cp1, ch1, cs1, cp2, ch2, cs2 !< Scalar field and Hubble equations: cp_i * phi'' + ch_i * H' + cs_i = 0, i=1,2
      real(dl) :: G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx !< Lagrangian coefficients
      real(dl) :: rho_tot, pres_tot, rhodot_tot, presdot_tot
      real(dl) :: energy_conservation
      integer  :: nu_i

      integer :: error_flag

      ! digest the input:
      if ( .not. present(only_solve) ) only_solve = .False.

      if ( self%evolve_hubble ) then
         num_eq = 4
      end if

      success = .True.

      self%maxfrac_hdsk = 0._dl
      self%z_maxfrac_hdsk = 0._dl

      ! initial time
      t1  = self%H%x(1)
      a   = exp(t1)

      ! Set initial conditions:
      allocate(y(num_eq))
      allocate(ydot(num_eq))
      y(1) = 1.0_dl/a/self%H_ini
      y(2) = self%phi_ini
      y(3) = self%phidot_ini
      if ( self%evolve_hubble ) y(4) = self%H_ini

      ! Initialize DLSODA:
      ! set-up the relative and absolute tollerances:
      itol = 1
      rtol = 1.d-10
      atol = 1.d-14
      ! initialize task to do:
      itask  = 1
      istate = 1
      iopt   = 1
      ! initialize the work space:
      LRN = 20 + 16*num_eq
      LRS = 22 + 9*num_eq + num_eq**2
      LRW = max(LRN,LRS)
      LIS = 20 + num_eq
      LIN = 20
      LIW = max(LIS,LIN)
      ! allocate the arrays:
      allocate(rwork(LRW))
      allocate(iwork(LIW))
      ! optional lsoda input:
      RWORK(5) = 0._dl  ! the step size to be attempted on the first step. The default value is determined by the solver.
      RWORK(6) = 0._dl  ! the maximum absolute step size allowed. The default value is infinite.
      RWORK(7) = 0._dl  ! the minimum absolute step size allowed. The default value is 0.
      IWORK(5) = 0      ! flag to generate extra printing at method switches. IXPR = 0 means no extra printing (the default). IXPR = 1 means print data on each switch.
      IWORK(6) = 1000   ! maximum number of (internally defined) steps allowed during one call to the solver. The default value is 500.
      IWORK(7) = 0      ! maximum number of messages printed (per problem) warning that T + H = T on a step (H = step size). This must be positive to result in a non-default value.  The default value is 10.
      IWORK(8) = 0      ! the maximum order to be allowed for the nonstiff (Adams) method.  the default value is 12.
      IWORK(9) = 0      ! the maximum order to be allowed for the stiff (BDF) method.  The default value is 5.
      ! additional lsoda stuff:
      CALL XSETF(0) ! suppress odepack printing
      ! Jacobian mode: 1=fullJacobian, 2=not provided
      JacobianMode = 2

      ! compute output EFT functions:
      call output( num_eq, 1, t1, y )

      ! solve the equations:
      do i=1, self%H%num_points-1

         ! set the time step:
         t1 = self%H%x(i)
         t2 = MIN(self%H%x(i+1), self%H%x_final)

         ! perform one integration step
         call DLSODA ( derivs, num_eq, y, t1, t2, itol, rtol, atol, itask, istate, iopt, RWORK, LRW, IWORK, LIW, jacobian, JacobianMode)
         ! check istate for LSODA good completion:
         if ( istate < 0 ) then
            if ( istate == -1 ) then
               istate = 1
            else
               success = .False.
               return
            end if
         end if

         ! compute output EFT functions if needed:
         call output( num_eq, i+1, t2, y )

      end do

      return

   contains

      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that computes y' given y for the background
      subroutine derivs( num_eq, x, y, ydot )

         implicit none

         integer , intent(in)                     :: num_eq !< number of equations in the ODE system
         real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed
         real(dl), intent(in) , dimension(num_eq) :: y      !< input status of the system
         real(dl), intent(out), dimension(num_eq) :: ydot   !< value of the derivative at time x

         !< local variables
         real(dl) :: E0, E1, E2, E3 !< Hubble constraint E0 + E1*H + E2*H^2 + E3*H^3 = 0

         ! 0) store variables:
         phi       = y(2)
         phi_prime = y(3)

         ! 1) convert x to a:
         a  = Exp(x)
         a2 = a*a

         ! 2) compute matter densities and pressure:
         grhob_t = params_cache%grhob/a         ! 8\pi G_N \rho_b a^2: baryons background density
         grhoc_t = params_cache%grhoc/a         ! 8\pi G_N \rho_{cdm} a^2: cold dark matter background density
         grhor_t = params_cache%grhornomass/a2  ! 8\pi G_N \rho_{\nu} a^2: massless neutrinos background density
         grhog_t = params_cache%grhog/a2        ! 8\pi G_N \rho_{\gamma} a^2: radiation background density

         grhonu_tot = 0._dl
         gpinu_tot  = 0._dl
         if ( params_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, params_cache%Nu_mass_eigenstates
               grhonu      = 0._dl
               gpinu       = 0._dl
               grhormass_t = params_cache%grhormass(nu_i)/a2
               call ThermalNuBack%rho_P( a*params_cache%nu_masses(nu_i), grhonu, gpinu)
               grhonu_tot = grhonu_tot + grhormass_t*grhonu ! 8\pi G_N \rho_mnu a^2: massive neutrino background density
               gpinu_tot  = gpinu_tot  + grhormass_t*gpinu  ! 8\pi G_N P_mnu a^2   : massive neutrino background pressure
            end do
         end if

         grho_matter  = grhob_t + grhoc_t + grhog_t + grhor_t + grhonu_tot
         gpres_matter = gpinu_tot + ( grhog_t + grhor_t )/3._dl

         rho_tot = grho_matter / a2
         pres_tot = gpres_matter / a2

         ! 3) compute Horndeski lagrangian functions
         call self%horndeski_lagrangian_coefficients( a, phi, 0.5_dl*(phi_prime/a)**2, G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx)

         ! 4) compute Hubble: E0 + E1*H + E2*H^2 + E3*H^3 = 0
         call horndeski_hubble_coefficients(a, phi_prime, rho_tot, G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx, E0, E1, E2, E3)
         if ( .not. self%evolve_hubble ) then
            if ( IsNaN(E0) .or. IsNaN(E1) .or. IsNaN(E2) .or. IsNaN(E3) .or. (E2==0._dl) .or. (E1*E1-4._dl*E2*E0<0._dl) ) then
               success = .False.
               return
            end if
            if ( abs(E3) < 1d-30 ) then
               ! use exact root when E3 = 0
               H = 0.5_dl * (-E1 + sqrt(E1 * E1 - 4._dl * E2 * E0)) / E2
            else
               write(*,*) "E3 = ", E3, "is non-zero."
               write(*,*) "You have G5x != 0 or G5xx != 0. Horndeski_evolve_hubble must be enabled in this case."
               call MpiStop('EFTCAMB error')
            end if
         else
            H = y(4)
         end if
         ! sanity check
         if (H <= 0._dl ) then
            success = .False.
            return
         end if
         
         adotoa = a*H

         ! 5) compute H' and phi'': cp_i * phi'' + ch_i * H' + cs_i = 0, i=1,2
         call horndeski_backeq_coefficients(a, H, phi_prime, rho_tot+pres_tot, G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx, cp1, ch1, cs1, cp2, ch2, cs2)

         phi_prime_prime = (ch2*cs1 - ch1*cs2)/(-(ch2*cp1) + ch1*cp2)
         H_prime = (cp2*cs1 - cp1*cs2)/(ch2*cp1 - ch1*cp2)

         ydot(1) = 1._dl / adotoa
         ydot(2) = phi_prime / adotoa
         ydot(3) = phi_prime_prime / adotoa
         if ( self%evolve_hubble ) then
            ! add friction term to enforce Hubble constraint, i.e. dynamically damp deviation from energy_conservation = 0
            energy_conservation = E0 + E1*H + E2*H*H + E3*H*H*H
            if ( E1 + 2._dl*E2*H + 3._dl*E3*H*H >= 0._dl ) then ! energy_conservation go up when H gets bigger
               ydot(4) = (H_prime - 5._dl * a * energy_conservation) / adotoa
            else
               ydot(4) = (H_prime + 5._dl * a * energy_conservation) / adotoa
            end if
         end if

      end subroutine

      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that computes the Jacobian of the system.
      subroutine jacobian( num_eq, x, y, ml, mu, pd, nrowpd )

         implicit none

         integer                            :: num_eq !< number of components of the Jacobian
         integer                            :: ml     !< ignored
         integer                            :: mu     !< ignored
         integer                            :: nrowpd !< ignored
         real(dl)                           :: x      !< time at which the Jacobian is computed
         real(dl), dimension(num_eq)        :: y      !< input status of the system
         real(dl), dimension(nrowpd,num_eq) :: pd     !< output Jacobian

      end subroutine jacobian

      ! ---------------------------------------------------------------------------------------------
      !> Subroutine that takes the solution of the background equations and computes additional quantities
      !! and the values of the EFT functions.
      subroutine output( num_eq, ind, x, y )

         implicit none

         integer , intent(in)                     :: num_eq !< number of equations in the ODE system.
         integer , intent(in)                     :: ind    !< index of the EFT functions interpolation tables to fill.
         real(dl), intent(in)                     :: x      !< time at which the derivatives of the system are computed.
         real(dl), intent(in) , dimension(num_eq) :: y      !< input status of the system.

         ! local variables
         logical :: is_open
         real(dl) :: cp1p, ch1p, cs1p, cp2p, ch2p, cs2p
         real(dl) :: alphaM, dalphaM, alphaB, dalphaB, alphaK, dalphaK, alphaT, dalphaT, M, dM, ddM

         ! 1) call derivs to make sure everything is initialized at the correct time step.
         ! this also computes a, a2, phi, phi_prime, phi_prime_prime, H, H_prime, adotoa, rho_tot, pres_tot, G_i's and c[hps]_i's
         call derivs( num_eq, x, y, ydot )
         ! write(*,*) 'a, phidot, x, y1, y2, H2, v1, v2', a, phi_prime, phi_prime/sqrt(6._dl)/H/a, sqrt(self%model_params(5)*exp(-phi*self%model_params(2)))/sqrt(3._dl)/H, sqrt(self%model_params(6)*exp(-phi*self%model_params(3)))/sqrt(3._dl)/H, H, self%model_params(5)*exp(-phi*self%model_params(2)), self%model_params(6)*exp(-phi*self%model_params(3))

         ! 2) update the horndeski maxfrac and redshift
         if (self%maxfrac_hdsk < 1._dl - (rho_tot+self%Lambda_CC)/(3._dl*H**2)) then
            self%maxfrac_hdsk = 1._dl - (rho_tot+self%Lambda_CC)/(3._dl*H**2)
            self%z_maxfrac_hdsk = 1._dl/a - 1._dl
         end if

         ! 3) log basic background quantities
         self%tau%y(ind) = y(1)
         self%H%y(ind) = H
         self%dH%y(ind) = H_prime
         self%rhom%y(ind) = rho_tot
         self%presm%y(ind) = pres_tot
         self%phi%y(ind) = phi
         self%dphi%y(ind) = phi_prime

         ! stop here if no further output are required
         if ( only_solve ) return

         ! 4) compute additional background quantities
         grhonu_tot = 0._dl
         gpinu_tot  = 0._dl
         grhonudot_tot = 0._dl
         gpinudot_tot  = 0._dl
         if ( params_cache%Num_Nu_Massive /= 0 ) then
            do nu_i = 1, params_cache%Nu_mass_eigenstates
               grhonu      = 0._dl
               gpinu       = 0._dl
               grhonudot   = 0._dl
               gpinudot    = 0._dl
               grhormass_t = params_cache%grhormass(nu_i)/a2
               call ThermalNuBack%rho_P( a*params_cache%nu_masses(nu_i), grhonu, gpinu)
               grhonu_tot = grhonu_tot + grhormass_t*grhonu ! 8\pi G_N \rho_mnu a^2: massive neutrino background density
               gpinu_tot  = gpinu_tot  + grhormass_t*gpinu  ! 8\pi G_N P_mnu a^2   : massive neutrino background pressure
               grhonudot = ThermalNuBack%drho(a*params_cache%nu_masses(nu_i), adotoa)
               grhonudot_tot = grhonudot_tot + grhormass_t*( grhonudot -4._dl*adotoa*grhonu)
               gpinudot = ThermalNuBack%pidot(a*params_cache%nu_masses(nu_i),adotoa, gpinu )
               gpinudot_tot  = gpinudot_tot  + grhormass_t*( gpinudot - 4._dl*adotoa*gpinu)
            end do
         end if

         rhodot_tot = -3._dl * adotoa * (rho_tot + pres_tot)
         presdot_tot = (- 4._dl/3._dl * adotoa * (grhog_t + grhor_t) + gpinudot_tot)/a2

         call horndeski_backeq_coefficients_p( a, H, H_prime, phi_prime, phi_prime_prime, rho_tot+pres_tot, rhodot_tot+presdot_tot, G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx, cp1p, ch1p, cs1p, cp2p, ch2p, cs2p )

         phi_prime_prime_prime = ((ch2p*cp1 + ch2*cp1p - ch1p*cp2 - ch1*cp2p)*(ch2*cs1 - ch1*cs2) - (ch2*cp1 - ch1*cp2)*(ch2p*cs1 + ch2*cs1p - ch1p*cs2 - ch1*cs2p))/(ch2*cp1 - ch1*cp2)**2
         H_prime_prime = ((ch2p*cp1 + ch2*cp1p - ch1p*cp2 - ch1*cp2p)*(-(cp2*cs1) + cp1*cs2) - (ch2*cp1 - ch1*cp2)*(-(cp2p*cs1) - cp2*cs1p + cp1p*cs2 + cp1*cs2p))/(ch2*cp1 - ch1*cp2)**2

         
         self%ddH%y(ind) = H_prime_prime
         self%rhomdot%y(ind) = rhodot_tot
         self%presmdot%y(ind) = presdot_tot
         self%ddphi%y(ind) = phi_prime_prime
         self%dddphi%y(ind) = phi_prime_prime_prime

         ! 5) compute the EFT function alphas
         call horndeski_eft_alphas(a, H, H_prime, H_prime_prime, phi_prime, phi_prime_prime, phi_prime_prime_prime, G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx, alphaM, dalphaM, alphaB, dalphaB, alphaK, dalphaK, alphaT, dalphaT, M, dM, ddM)
         self%EFTalphaM%y(ind) = alphaM
         self%EFTalphaMdot%y(ind) = dalphaM
         self%EFTalphaK%y(ind) = alphaK
         self%EFTalphaKdot%y(ind) = dalphaK
         self%EFTalphaB%y(ind) = alphaB
         self%EFTalphaBdot%y(ind) = dalphaB
         self%EFTalphaT%y(ind) = alphaT
         self%EFTalphaTdot%y(ind) = dalphaT
         self%EFTM%y(ind) = M
         self%EFTMdot%y(ind) = dM
         self%EFTMdotdot%y(ind) = ddM

         ! 6) computation of EFT gammas and the delta_phi coefficients are done seperately after the background integration

      end subroutine

   end subroutine EFTCAMBHdskSolveBackgroundEquations


   ! ---------------------------------------------------------------------------------------------
   !> Function that computes dtauda = 1/a**2/H.
   function EFTCAMBHdskComputeDtauda( self, a, eft_par_cache, eft_cache )
      implicit none
      class(EFTCAMB_Hdsk)                           :: self          !< the base class.
      real(dl), intent(in)                          :: a             !< the input scale factor.
      type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
      type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.
      real(dl)                                      :: EFTCAMBHdskComputeDtauda !< the output dtauda

      real(dl) :: x, mu, a1
      integer  :: ind

      x = Log(a)

      if (x > self%H%x(1)) then
         call self%H%precompute(x, ind, mu )
         EFTCAMBHdskComputeDtauda = 1._dl / a / a / self%H%value( x, index=ind, coeff=mu )   
      else 
         a1 = exp(self%H%x(1))
         EFTCAMBHdskComputeDtauda = 1._dl / a1 / a1 / self%H%y(1)
      end if
   
   end function EFTCAMBHdskComputeDtauda


   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that computes adotoa = a*H.
   subroutine EFTCAMBHdskComputeAdotoa( self, a, eft_par_cache, eft_cache )
      implicit none
      class(EFTCAMB_Hdsk)                           :: self          !< the base class.
      real(dl), intent(in)                          :: a             !< the input scale factor.
      type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
      type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

      real(dl) :: x, mu, a1
      integer  :: ind

      x = Log(a)

      if (x > self%H%x(1)) then
         call self%H%precompute(x, ind, mu )
         ! remeber what we have are the physical H and its conformal derivatives
         eft_cache%adotoa = a*self%H%value( x, index=ind, coeff=mu )
      ! extrapolate
      else
         a1 = exp(self%H%x(1))
         eft_cache%adotoa = (a1/a)*(a1*self%H%y(1))
      end if

   end subroutine EFTCAMBHdskComputeAdotoa


   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that computes the two derivatives wrt conformal time of H.
   subroutine EFTCAMBHdskComputeHubbleDer( self, a, eft_par_cache, eft_cache )
      implicit none
      class(EFTCAMB_Hdsk)                           :: self          !< the base class.
      real(dl), intent(in)                          :: a             !< the input scale factor.
      type(TEFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
      type(TEFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

      real(dl) :: x, mu, H, dH, ddH, dddH, a1
      integer  :: ind

      x = Log(a)

      if (x > self%H%x(1)) then
         call self%H%precompute( x, ind, mu )
         ! remeber what we have are the physical H and its conformal derivatives
         H = self%H%value( x, ind, mu )
         dH = self%dH%value( x, index=ind, coeff=mu )
         ddH = self%ddH%value( x, index=ind, coeff=mu )
         ! we should not need Hdotdotdot for the positivity bound of Horndeski
         dddH = a*H*self%ddH%first_derivative( x, index=ind, coeff=mu ) ! d/dx -> d/dtau
         eft_cache%Hdot =  a*dH + (a*H)**2
         eft_cache%Hdotdot = a*ddH + 3._dl*a**2*H*dH + 2._dl*(a*H)**3
         eft_cache%Hdotdotdot = a*dddH + 4._dl*a**2*H*ddH + 3._dl*a**2*dH**2 + 12._dl*a**3*H**2*dH + 6._dl*(a*H)**4
      else
         a1 = exp(self%H%x(1))
         H = self%H%y(1)
         dH = self%dH%y(1)
         ddH = self%ddH%y(1)
         ! we should not need Hdotdotdot for the positivity bound of Horndeski
         dddH = a1*H*self%ddH%yp(1) ! d/dx -> d/dtau
         eft_cache%Hdot =  (a1/a)**2*(a1*dH + (a1*H)**2)
         eft_cache%Hdotdot = (a1/a)**3*(a1*ddH + 3._dl*a1**2*H*dH + 2._dl*(a1*H)**3)
         eft_cache%Hdotdotdot = (a1/a)**4*(a1*dddH + 4._dl*a1**2*H*ddH + 3._dl*a1**2*dH**2 + 12._dl*a1**3*H**2*dH + 6._dl*(a1*H)**4)
      end if

   end subroutine EFTCAMBHdskComputeHubbleDer


   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that finds the model specific initial conditions.
   subroutine EFTCAMBHdskFindInitialConditions( self, params_cache, feedback_level, success )

      implicit none

      class(EFTCAMB_Hdsk)                          :: self           !< the base class
      type(TEFTCAMB_parameter_cache), intent(in)   :: params_cache   !< a EFTCAMB parameter cache containing cosmological parameters.
      integer                      , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
      logical , intent(out)                        :: success        !< whether the calculation ended correctly or not

      integer :: nu_i

      real(dl) :: a_ini, H_ini
      real(dl) :: rhob_t, rhoc_t, rhog_t, rhor_t, grhonu_tot, grhonu, gpinu, grhormass_t
      real(dl) :: G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx

      success = .True.

      a_ini = exp( self%x_initial )

      ! Estimate inital Hubble, assuming radiation dominance
      rhob_t = params_cache%grhob/a_ini**3        
      rhoc_t = params_cache%grhoc/a_ini**3        
      rhor_t = params_cache%grhornomass/a_ini**4  
      rhog_t = params_cache%grhog/a_ini**4         
      grhonu_tot = 0._dl
      if ( params_cache%Num_Nu_Massive /= 0 ) then
         do nu_i = 1, params_cache%Nu_mass_eigenstates
            grhonu      = 0._dl
            gpinu       = 0._dl
            grhormass_t = params_cache%grhormass(nu_i)/a_ini**2
            call ThermalNuBack%rho_P( a_ini*params_cache%nu_masses(nu_i), grhonu, gpinu)
            grhonu_tot = grhonu_tot + grhormass_t*grhonu
         end do
      end if
      H_ini = sqrt((rhob_t + rhoc_t + rhog_t + rhor_t + grhonu_tot / a_ini**2)/3._dl)
      
      self%H_ini = H_ini

      ! assume field frozen at first, i.e. phidot = 0
      call self%horndeski_lagrangian_coefficients( a_ini, self%phi_ini, 0._dl, G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx)

      ! use the input phidot_ini when model specific IC is not requested
      if ( .not. self%model_specifc_ic ) return

      ! set the IC depending on model
      select case (self%model)
       case(1)
         ! tracking solution H dphi/dt = -V_phi/5, same as canonical scalar field
         self%phidot_ini = a_ini * G2p / 5._dl / H_ini
       case(2)
         ! tracking solution H dphi/dt = -V_phi/5, same as canonical scalar field
         self%phidot_ini = a_ini * G2p / 5._dl / H_ini
       case(3)
         ! tracking solution H dphi/dt =  (-1 + sqrt(1 - 12/5*(xi*V_phi)))/(6 xi)
         if (abs(self%model_params(1)*G2p) > 0.01 ) then
            self%phidot_ini = a_ini * (-1._dl + sqrt(1._dl + 2.4_dl * self%model_params(1) * G2p)) / (6._dl * self%model_params(1) * H_ini)
         else
            self%phidot_ini = a_ini * G2p / 5._dl / H_ini
         end if
         ! should be accurate enough for now, TODO: shooting
         self%Lambda_CC = params_cache%grhov
       case(4)
         ! tracking solution H dphi/dt = -V_phi/5, same as canonical scalar field
         self%phidot_ini = a_ini * G2p / 5._dl / H_ini
       case default
         write(*, '(a,I3)') "Model specific initial condition not implemented for Model #", self%model
         write(*, *) "Default to phidot_ini = 0"
         write(*, *) "If it is OK for you, set Horndeski_model_specific_ic to False and pass in Honrdeski_phidot_ini explicitly."
         self%phidot_ini = 0._dl
      end select

   end subroutine EFTCAMBHdskFindInitialConditions

   ! ---------------------------------------------------------------------------------------------
   !> Subroutine that computes all needed honrdeski functions given phi and X
   subroutine EFTCAMBHdskLagrangianCoefficients(self, a, phi, X, G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx)
      
      implicit none

      class(EFTCAMB_Hdsk)          :: self         !< the base class

      real(dl), intent(in)         :: a            !<- the scale factor
      real(dl), intent(in)         :: phi          !<- the scalar field value in units of Mp
      real(dl), intent(in)         :: X            !<- the scalar field kinetic term \f$ -1/2 \phi_\mu \phi^\nu \f$          
      real(dl), intent(out)        :: G2,G2p,G2x,G2pp,G2px,G2xx,G2ppp,G2ppx,G2pxx,G2xxx,G3p,G3x,G3pp,G3px,G3xx,G3ppp,G3ppx,G3pxx,G3xxx,G4,G4p,G4x,G4pp,G4px,G4xx,G4ppp,G4ppx,G4pxx,G4xxx,G4pppp,G4pppx,G4ppxx,G4pxxx,G4xxxx,G5,G5p,G5x,G5pp,G5px,G5xx,G5ppp,G5ppx,G5pxx,G5xxx,G5pppp,G5pppx,G5ppxx,G5pxxx,G5xxxx !<- Horndeski Lagrangian functions

      ! model specific convenient variables
      real(dl) :: v0, v1, v2, n, m2, f, xi, lambda, theta, beta1, beta2

      G2 = 0._dl
      G2p = 0._dl
      G2x = 1._dl ! canonical scalar field
      G2pp = 0._dl
      G2px = 0._dl
      G2xx = 0._dl
      G2ppp = 0._dl
      G2ppx = 0._dl
      G2pxx = 0._dl
      G2xxx = 0._dl
      G3p = 0._dl
      G3x = 0._dl
      G3pp = 0._dl
      G3px = 0._dl
      G3xx = 0._dl
      G3ppp = 0._dl
      G3ppx = 0._dl
      G3pxx = 0._dl
      G3xxx = 0._dl
      G4 = 0.5_dl ! standard GR value
      G4p = 0._dl
      G4x = 0._dl
      G4pp = 0._dl
      G4px = 0._dl
      G4xx = 0._dl
      G4ppp = 0._dl
      G4ppx = 0._dl
      G4pxx = 0._dl
      G4xxx = 0._dl
      G4pppp = 0._dl
      G4pppx = 0._dl
      G4ppxx = 0._dl
      G4pxxx = 0._dl
      G4xxxx = 0._dl
      G5 = 0._dl
      G5p = 0._dl
      G5x = 0._dl
      G5pp = 0._dl
      G5px = 0._dl
      G5xx = 0._dl
      G5ppp = 0._dl
      G5ppx = 0._dl
      G5pxx = 0._dl
      G5xxx = 0._dl
      G5pppp = 0._dl
      G5pppx = 0._dl
      G5ppxx = 0._dl
      G5pxxx = 0._dl
      G5xxxx = 0._dl

      ! power law potential
      if (self%model == 1) then

         v0 = self%model_params(1)
         n = self%model_params(2)

         G2 = X - v0 * phi**n
         G2x = 1._dl
         if (n == 1) then
            G2p = -v0
            G2pp = 0._dl
            G2ppp = 0._dl
         else if (n == 2) then
            G2p = -2._dl*v0*phi
            G2pp = -2._dl*v0
            G2ppp = 0._dl
         else if (n == 3) then
            G2p = -3._dl*v0*phi*phi
            G2pp = -6._dl*v0*phi
            G2ppp = -6._dl*v0
         else
            G2p = -n * v0 * phi**(n - 1._dl)
            G2pp = -n * (n - 1._dl) * v0 * phi**(n - 2._dl)
            G2ppp = -n * (n - 1._dl) * (n - 2._dl) * v0 * phi**(n - 3._dl)
         end if

      ! axion-like potential
      else if (self%model == 2) then
         
         m2 = self%model_params(1)
         f = self%model_params(2)
         n = self%model_params(3)

         theta = phi/f

         G2 = X - m2 * f * f * (1._dl - cos(theta))**n
         G2x = 1._dl
         if (n==1) then
            G2p = -f*m2*sin(theta)
            G2pp = -m2*cos(theta)
            G2ppp = m2*sin(theta)/f
         else if (n==2) then
            G2p = -2._dl*f*m2*(1._dl - cos(theta))*sin(theta)
            G2pp = -2._dl*m2*(1._dl - cos(theta))*cos(theta) - 2._dl*m2*sin(theta)**2
            G2ppp = (2._dl*m2*(1._dl - cos(theta))*sin(theta))/f - (6._dl*m2*cos(theta)*sin(theta))/f
         else if (n==3) then
            G2p = -3._dl*f*m2*(1._dl - cos(theta))**2*sin(theta)
            G2pp = -3._dl*m2*(1._dl - cos(theta))**2*cos(theta) - 6._dl*m2*(1._dl - cos(theta))*sin(theta)**2
            G2ppp = (3._dl*m2*(1._dl - cos(theta))**2*sin(theta))/f - (18._dl*m2*(1._dl - cos(theta))*cos(theta)*sin(theta))/f - (6._dl*m2*sin(theta)**3)/f
         else
            G2p = -(f*m2*n*(1._dl - cos(theta))**(-1._dl + n)*sin(theta))
            G2pp = -(m2*n*(1._dl - cos(theta))**(-1._dl + n)*cos(theta)) - m2*(-1._dl + n)*n*(1._dl - cos(theta))**(-2._dl + n)*sin(theta)**2
            G2ppp = (m2*n*(1._dl - cos(theta))**(-1._dl + n)*sin(theta))/f - (3._dl*m2*(-1._dl + n)*n*(1._dl - cos(theta))**(-2._dl + n)*cos(theta)*sin(theta))/f - (m2*(-2._dl + n)*(-1._dl + n)*n*(1._dl - cos(theta))**(-3._dl + n)*sin(theta)**3)/f
         end if

      ! XboxPhi operator
      else if (self%model == 3) then

         xi = self%model_params(1)
         v0 = self%model_params(2)

         G2 = X - v0 * phi * phi * phi * phi - self%Lambda_CC
         G2x = 1._dl
         G2p = -4._dl * v0 * phi * phi * phi
         G2pp = -12._dl * v0 * phi * phi
         G2ppp = -24._dl * v0 * phi

         G3x = xi
      
      ! X^2 operator
      else if (self%model == 4) then
      
         xi = self%model_params(1)
         v0 = self%model_params(2)

         ! P(X) = X + xi X^2 - v0*phi^4
         G2 = X + xi * X * X - v0 * phi * phi * phi * phi
         G2x = 1._dl + 2._dl * xi * X
         G2xx = 2._dl * xi
         G2p = -4._dl * v0 * phi * phi * phi
         G2pp = -12._dl * v0 * phi * phi
         G2ppp = -24._dl * v0 * phi

      ! phi^2 X operator
      else if (self%model == 5) then
      
         xi = self%model_params(1)
         v0 = self%model_params(2)

         ! P(X) = (1 + xi*phi^2) X - V
         G2 = X + xi * phi * phi * X - v0 * phi * phi * phi * phi
         G2x = 1._dl + xi * phi * phi
         G2p = 2._dl * xi * phi * X - 4._dl * v0 * phi * phi * phi
         G2pp = 2._dl * xi * X - 12._dl * v0 * phi * phi
         G2px = 2._dl * xi
         G2ppp = -24._dl * v0 * phi

      ! scaling cubic galileon, arXiv:2112.06892
      else if (self%model == 6) then
      
         v0 = self%model_params(1)
         beta1 = self%model_params(2)
         beta2 = self%model_params(3)
         lambda = self%model_params(4)
         v1 = self%model_params(5)*exp(-beta1*phi)
         v2 = self%model_params(6)*exp(-beta2*phi)

         G2 = X - v1 - v2
         G2x = 1._dl
         G2p = beta1*v1 + beta2*v2
         G2pp = -beta1**2*v1 - beta2**2*v2
         G2ppp = beta1**3*v1 + beta2**3*v2

         G3p = lambda*v0
         G3x = v0/X
         G3xx = -v0/X/X
         G3xxx = 2._dl*v0/X/X/X

      else

         write(*,'(a,I3)') 'Undefined Horndeski model', self%model
         call MpiStop('EFTCAMB error')

      end if

   end subroutine EFTCAMBHdskLagrangianCoefficients

end module EFTCAMB_Horndeski
