c------------------------------------------------------------------------
c initialize AfkQed common block variables as beget would have:
c**********************************************************************
c   IR - 1-pi0, 2-eta, 3-eta'
c   IMODE - decay mode, 0 - no decay
c   for IR=1
c     1 - pi0 -> 2gamma
c   for IR=2
c     1 - eta -> 2gamma
c     2 - eta -> 3pi0
c     3 - eta -> pi+ pi- pi0
c     4 - eta -> pi+ pi- gamma
c   for IR=3
c     1 - eta-prime -> 2gamma
c     2 - eta-prime -> pi+ pi- eta, eta -> 2gamma
c     3 - eta-prime -> pi+ pi- eta, eta -> 3pi0
c     4 - eta-prime -> pi+ pi- eta, eta -> pi+ pi- pi0
c     5 - eta-prime -> pi+ pi- eta, eta -> pi+ pi- gamma
c     6 - eta-prime -> pi0 pi0 eta, eta -> 2gamma
c     7    - eta-prime -> pi0 pi0 eta, eta -> 3pi0
c     8    - eta-prime -> pi0 pi0 eta, eta -> pi+ pi- pi0
c     9    - eta-prime -> pi0 pi0 eta, eta -> pi+ pi- gamma
c     10   - eta-prime -> rho0 gamma
c     11   - eta-prime -> pi+ pi- pi0
c   for IR=6,7
c     1 - f1 -> a0(980)pi -> eta pi+pi-
c     2 - f1 -> rho gamma -> pi+pi-gamma
c     3 - f1 -> eta f0(500) -> eta pi+pi-
c**********************************************************************  
      program test

      implicit none
      Double Precision EB0,FMAX0,rmax0,rmin0,t2max0
      integer IR0,IMODE0,KVDM0,RAD0,i
c---------------------------------------------------------------------      
      EB0=10.58/2     !10.58/2.
c      EB0=189./2    ! L3 at LEP
c      IR0=6          !  f1 m=1
      IR0=7          
      IMODE0=0       !  number of decay mode
      KVDM0=1        !  KVDM factor
c      FMax0=0.2      ! for IR0  =7  for notag MC
c      FMax0=1. ! for IR0  =6  for notag MC 
c      FMax0=0.02      ! for IR0  =6  for t2max0=-1.5d0
      FMax0=1.6      ! for IR0  =7  for t2max0=-1.5d0
c      RAD0=0         ! for notag MC   
      RAD0=1
      rmax0=0.5
      rmin0=1.d-4
c     t2max0=-0.0001
      t2max0=-2.
      FMax0=0.05
      RAD0=1
      call GGRESI(EB0,IR0,IMODE0,KVDM0,RAD0,FMAX0,
     &  rmax0,rmin0,t2max0)

      do i=1,10000
        call GGRESPS(EB0)
c        if(i/100*100.eq.i)print *,i
      enddo  

      CALL GGREND
      Return
      End

      include 'ggrps.f'
      include 'ggfsr.f'
      include 'ggrspc.f'

