!! This module contains preprocessor definitions.  It is used to remove
!  the need for the preprocessor in normal files.
!

module defines_mod


!! from units.h
!
#if EXTERNAL_POT>0
      logical, parameter :: external_pot = .true.
#else
      logical, parameter :: external_pot = .false.
#endif
#if FLUID_RESTART>0
      logical, parameter :: restart_fluid = .true.
#else
      logical, parameter :: restart_fluid = .false.
#endif
#if WIGGLE_RESTART>0
      logical, parameter :: restart_wiggle = .true.
#else
      logical, parameter :: restart_wiggle = .false.
#endif
#if ROTATING>0
      real*8, parameter :: omega_frame=5d-2  
#endif


!! from hydroparam.h
!
! change the following for passive array size modification.

#if PASSIVE>0
      integer, parameter :: PAS=PASSIVE
#else
      integer, parameter :: PAS = 0    ! give a default value
#endif


!! from globals.h
!
#if PASSIVE>0
      real*8 :: passflux
      common/passivearray/passflux(JMAX+2,KMAX+2,LMAX,PAS)
#endif


end module defines_mod
