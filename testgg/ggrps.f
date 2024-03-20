      SUBROUTINE GGRESI(EB0,IR0,IMODE0,KVDM0,RAD0,
     &  FMAX0,rmax0,rmin0,t2max0)
c**********************************************************************
c                                          pc  -+    
c   MC generator of event  ee -> ee + R ( J  =0  )  
c                                                   
c   Input parameters:                               
c                                                   
c   EB - Energy of beam, GeV                        
c   IR - 1-pi0, 2-eta, 3-eta', 4-etac, 5-etab, 6-f1(TS+ST), 7-f1(TT)
c   IMODE - decay mode, 0 - no decay                                
c   for IR=1
c     1 - pi0 -> 2gamma
c     2 - pi0 -> e+ e- gamma
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
c     7 - eta-prime -> pi0 pi0 eta, eta -> 3pi0
c     8 - eta-prime -> pi0 pi0 eta, eta -> pi+ pi- pi0
c     9 - eta-prime -> pi0 pi0 eta, eta -> pi+ pi- gamma
c    10 - eta-prime -> rho0 gamma
c    11 - eta-prime -> pi+ pi- pi0
c   for IR=6,7
c     1 - f1 -> a0(980)pi -> eta pi+pi-
c     2 - f1 -> rho gamma -> pi+pi-gamma 
c     3 - f1 -> eta f0(500) -> eta pi+pi-
c      
c   KVDM - key  for Vector Dom. Model: 1 - included 
c   FMAX - max value of cross section
c                                                                     
c   Result of simulation: in common-block /GGREV/                    
c                                                                    
c   PPART(1-3,i) - Momenta (GeV/c) of produced particles             
c   PPART(4,i)   - their total energy (GeV)                          
c**********************************************************************
      IMPLICIT NONE
C#include "AfkQed/ggrps.inc"
      include 'ggrps.inc'
      real*8 EB0,FMAX0,rmax0,rmin0,t2max0
      integer IR0,IMODE0,KVDM0,RAD0
      REAL*8 RMASS(7),RGGG(7),RWIDTH(7)
      data RMASS/
     &  0.1349766D0,0.54751D0,0.95778D0,2.9804D0,9.4,1.2819,1.2819/
      data RGGG/0.007742D0,0.510D0,4.30D0,7.14D0,0.4,2.69,2.69/
      Data RWIDTH/0.0d0,0.0d0,0.0d0,34.3d-3,30.0d-3,22.7d-3,22.7d-3/
      integer i
c      external afkrndm_ranmar
      integer LUX,INT,k1,k2
      save

      OPEN(UNIT=13,FILE='vpol10.dat',STATUS='Old',ERR=101)
      goto 110
 101  Close(13)
      OPEN(UNIT=13,FILE='RELEASE/AfkQed/vpol10.dat',
     &  STATUS='Old',ERR=102)
      goto 110
 102  Close(13)
      OPEN(UNIT=13,FILE='PARENT/AfkQed/vpol10.dat',
     &  STATUS='Old',ERR=103)
      goto 110
 103  print *,'GGRESI: Could not open vpol10.dat'
      stop
      
 110  do i=1,7330
        read(13,*) sets(i),setpol(i)
      enddo
      Close(13)

      FMAX=FMAX0
      IR=IR0
      IMODE=IMODE0
      KVDM=KVDM0
      EB=EB0
      RAD=RAD0
      RMAS=RMASS(IR)
      RWID=RWIDTH(IR)
      RG=RGGG(IR)
      rmax=rmax0
      rmin=rmin0
      Q2max=t2max0

      PI=DACOS(-1.D0)
      alpha=1./137.03604
      EM=0.51099892D-3
      mpi0=0.1349766D0
      mpi=0.13957018d0
      meta=0.54751D0
      Metap=0.95778
      Mks=0.497648
      Mkc=0.493677
      mrho=0.77526d0
      Mjpsi=3.096916
      Mups=9.4603
      ma0=0.9825
      wa0=0.0756
      mf0=0.510
      wf0=0.350
      wrho=0.1491
      g12=0.3
      
      SUM=0.
      ES=0.
      FM=0.
      NGT=0
      NOBR=0
      Nact=0
      
      PRINT 29, EB0,RMAS,RWID,RG,KVDM,RAD,rmax,Q2max,FMAX,IMODE
29    FORMAT(/,80(1H*),//,
     *' Simualtion of e+e- -> e+e- +R(0-+). 09.01.2008',/,
     *' ------------------------------------------------------',//,
     * ' Ebeam = ',F6.3,' GeV',/,
     * ' Mass of resonance  = ',F6.3,' GeV/c**2',/,
     * ' Total width of resonance  = ',F10.7,' GeV/c**2',/,
     * ' Resonance gg width = ',F6.3,' keV',/,
     * ' KVDM = ',I2,/,
     * ' RAD = ',I2,/,
     * ' max. Eph/Ebeam = ',F6.3,/,
     * ' Q2max = ',F10.6,' GeV^2',/,  
     * ' FMAX = ', 1pE10.2,/,
     * ' Resonance decay mode = ',I3,':',/)

c
      IF(IR.EQ.1) THEN
        IF(IMODE.EQ.1)then
          print *, ' pi0 -> 2gamma'
        else
          IMODE=0
        endif  
      END IF
c
      IF(IR.EQ.2) THEN
        IF(IMODE.EQ.1) then
          print *, ' eta -> 2gamma'
        else IF(IMODE.EQ.2) then
          print *, ' eta -> 3pi0'
        else IF(IMODE.EQ.3) then
          print *, ' eta -> pi+ pi- pi0'
        else IF(IMODE.EQ.4) then
          print *, ' eta -> pi+ pi- gamma'
        else
          IMODE=0
        endif
      END IF
c
      IF(IR.EQ.3) THEN
        IF(IMODE.EQ.1) then
          print *,'eta-prime -> 2gamma'
        else IF(IMODE.EQ.2) then
          print *,'eta-prime -> pi+pi-eta (eta -> 2gamma)'
        else IF(IMODE.EQ.3) then
          print *,'eta-prime -> pi+pi-eta (eta -> 3pi0)'
        else IF(IMODE.EQ.4) then
          print *,'eta-prime -> pi+pi-eta (eta -> pi+pi-pi0)'
        else IF(IMODE.EQ.5) then
          print *,'eta-prime -> pi+pi-eta (eta -> pi+pi-gamma)'
        else IF(IMODE.EQ.6) then
          print *,'eta-prime -> pi0pi0eta (eta -> 2gamma)'
        else IF(IMODE.EQ.7) then
          print *,'eta-prime -> pi0pi0eta (eta -> 3pi0)'
        else IF(IMODE.EQ.8) then
          print *,'eta-prime -> pi0pi0eta (eta -> pi+pi-pi0)'
        else IF(IMODE.EQ.9) then
          print *,'eta-prime -> pi0pi0eta (eta -> pi+pi-gamma)'
        else IF(IMODE.EQ.10) then
          print *,'eta-prime -> rho0gamma (rho0 ->pi+pi-)'
        else
          IMODE=0
        END IF
      endif
      IF(IR.EQ.4) THEN
        if(IMODE.EQ.1) then
          print *,' eta_c -> KS K+ pi-'
        else
          IMODE=0
        endif
      endif
      IF(IR.EQ.6) THEN
        if(IMODE.EQ.1) then
          print *,' f1 -> a0 pi -> eta pi+pi- (M=1)'
        else if(IMODE.EQ.2) then
          print *,' f1 -> rho gamma (M=1)'
        else if(IMODE.EQ.3) then
          print *,' f1 -> f0(500)eta -> eta pi+pi- (M=1)'
        else
          IMODE=0
        endif
      END IF
      IF(IR.EQ.7) THEN
        if(IMODE.EQ.1) then
          print *,' f1 -> a0 pi -> eta pi+pi- (M=0)'
        else if(IMODE.EQ.2) then
          print *,' f1 -> rho gamma (M=0)'
        else if(IMODE.EQ.3) then
          print *,' f1 -> f0(500)eta -> eta pi+pi- (M=0)'
        else
          IMODE=0
        endif
      endif
      print 21
21    FORMAT(/,80(1H*))
      
c      call afkrndminit( afkrndm_ranmar )
      lux=3
      int=6789
      k1=0
      k2=0
      call RLUXGO(LUX,INT,k1,k2)

      end

      SUBROUTINE GGREND
      IMPLICIT none
c#include "AfkQed/ggrps.inc"
      include 'ggrps.inc'
      REAL*8 ERSEC,TSEC
c
      ERSEC= ES-SUM**2/NOBR
      IF(ERSEC.LT.0.)ERSEC=0.
      ERSEC= SQRT(ERSEC)/NOBR
      TSEC=SUM/NOBR
      print 3,NOBR,NACT,NGT,FM,TSEC,ERSEC
3     FORMAT(/,' Results of calculation:',/,
     &         ' -----------------------',//,
     & ' NOBR=',I20,2X,//,
     & ' NACT=',I20,2X,//,
     & ' NGT=',I6,2X,//,
     & ' Fm=',D10.2,//,
     &  ' Cross sec. = (',D10.4,' +-',D8.2,') nb',/)
c
      open(1, file='cross.txt', status='old', action='write')
      write(1, *) TSEC
      close(1)
c      
      end
c---------------------------------------------------------------------
      SUBROUTINE GGRESPS(EB0)
      IMPLICIT none
c#include "AfkQed/ggrps.inc"
      include 'ggrps.inc'
      REAL*8 EB0
      REAL*8 R0
      parameter (R0=2.817940325D-13)
      REAL*8 E,S,T1,T2,S1,S2,WEI,W,q2maxp
      REAL*8 P1P2,P1K2,P2K1,AK1K2,B1,B2,B3,B,A
      REAL*8 FAC,FAC1,FAC2,FAC3,SS,FI
      REAL*8 RND,bw
      REAL*8 GGRFVDM,ggpol
      real*8 nu,X,u1,u2,cophi,GGRD4,dd
      real*8 ll,rr,beta,beta0,C(3),GG(3),ptot,gamma
      Real*8  max,min, ptl(4),pa,pa0
      data GG/0.,0.,1./
      integer i

3     NOBR=NOBR+1        ! counter of calls
      if (RAD.eq.1) then
        Call GRNDM(rnd,1)
        beta0=alpha/pi*(dlog(-Q2max*(1.-rmax)/(em**2))-1.)
        rr=rmax*exp(dlog(rnd)/beta0)
      else
        beta0=0.
        rr=0.
      endif
      if(rr.gt.rmin)then
        S=4.d0*Eb0**2*(1.d0-rr)
        Q2maxp=q2max*(1.d0-rr)
      else
        S=4.d0*Eb0**2
        Q2maxp=q2max
      endif  
      E=DSQRT(s)/2.d0  
c generate resonance mass using Breit-Wigner distribution
      if(ir.le.3)then
        rm=rmas
      else if(ir.le.5)then
        min=atan(((rmas-15.*rwid)**2-rmas**2)/rmas/rwid)/pi
        max=atan(((rmas+15.*rwid)**2-rmas**2)/rmas/rwid)/pi
        CALL GRNDM(RND,1)
        rnd=(max-min)*rnd+min
        rm=sqrt(rmas**2+rmas*rwid*tan(pi*rnd))
      else
        IF(IMODE.eq.1) then
          CALL GGRf1M1(3,0)
        else if(IMODE.eq.2)then
          CALL GGRf1M2(3,0)
        else if(IMODE.eq.3)then
          CALL GGRf1M3(3,0)
        else  
   4      CALL GRNDM(RND,1)
          min=ma0+mpi
          max=rmas+15.*rwid
          rm=(max-min)*rnd+min
          if(rm.ge.ma0+mpi)then
            pa=sqrt((1.-(ma0+mpi)**2/rm**2)*
     &        (1.-(ma0-mpi)**2/rm**2))
            pa0=sqrt((1.-(ma0+mpi)**2/rmas**2)*
     &        (1.-(ma0-mpi)**2/rmas**2))
            bw=(rmas*rwid)**2/((rm**2-rmas**2)**2+(rmas*rwid)**2)*
     &        (pa/pa0)**3
          else
            bw=0.
          endif
          CALL GRNDM(RND,1)
          if(rnd*1.01.gt.bw)goto 4
        endif
      endif  

      CALL GGRINV(S,RM,T1,T2,S1,S2,WEI,Q2maxp)
      ss=0.
      if(wei.ne.0.)then
        W=RM
        P1P2=S/2.-EM**2
        P1K2=1./2.*(S1-T2-EM**2)
        P2K1=1./2.*(S2-T1-EM**2)
        AK1K2=1./2.*(RM**2-T1-T2)
c**** B from paper of H.Terazava, Rev.Mod.Phys. 43(1973)615
        if(ir.le.5)then
c        B1=(4.*P1P2-2.*P1K2-2.*P2K1+AK1K2)**2+
c     &    AK1K2**2-T1*T2-16.*EM**4
c        B2=P1P2*AK1K2-P1K2*P2K1
c        B3=T1*(2*P1K2-AK1K2)**2+T2*(2*P2K1-AK1K2)**2+
c     &    4.*EM**2*AK1K2**2
c        B=1./4.*T1*T2*B1-4.*B2**2+EM**2*B3
c using Budnev's formalizm 
C sigma_TT=8pi^2(Gamma_{2gamma}/RM) A FAC2 delta(rm**2-rmas**2)
c tau_TT=-2sigma_TT        
c        A=2\sqrt{X}/RM^2
c b=(2rho_1^{++}2rho_2^{++} t1t2 A +
c    2|rho_1^{+-}rho_2^{+-}| t1t2 cos(2phi) (-2A)) \sqrt(X) RM**2/16        
          nu=(w**2-t1-t2)/2.
          X=nu**2-t1*t2
          u1=s2-em**2-t1
          u2=s1-em**2-t2
          b=((u2-nu)**2/X*t1+t1+4.*em**2)*
     &      ((u1-nu)**2/X*t2+t2+4.*em**2)*X/8.
          cophi=(-2.*s+u1+u2-nu+4*em**2+
     &      nu*(u2-nu)*(u1-nu)/x)/
     &      sqrt(t1*t2*((u2-nu)**2/X-1.+4.*em**2/t1)*
     &      ((u1-nu)**2/X-1.+4.*em**2/t2))
          cophi=cophi**2*2.-1.
          B=B-
     &      abs(((u2-nu)**2/X*t1-t1+4.*em**2))*
     &      abs(((u1-nu)**2/X*t2-t2+4.*em**2))*
     &      cophi*X/8.
        else if(ir.eq.6)then
c f1 polarization =1 sigma_TS and sigma_ST
c sigma_TS=8pi^2(Gamma_{2gamma}/RM) A FAC2 delta(rm**2-rmas**2)
c A_TS=3*(t1-nu)**2/X*(-t2)/sqrt(X)        
c sigma_ST=8pi^2(Gamma_{2gamma}/RM) A FAC2 delta(rm**2-rmas**2)
c A_ST=3*(t2-nu)**2/X*(-t1)/sqrt(X)
c tau_TS=8pi^2(Gamma_{2gamma}/RM) A FAC2 delta(rm**2-rmas**2)
c A_INF=3./2.*(t1-nu)*(t2-nu)/X*sqrt(t1*t2/X)        
c b=(2rho_1^{++}rho_2^{00} t1t2 A_TS + (rho_1^{00}2rho_2^{++} t1t2 A_ST +
c    -8|rho_1^{+0}rho_2^{+0}| t1t2 cos(phi) A_INF) \sqrt(X) RM**2/16        
          nu=(w**2-t1-t2)/2.
          X=nu**2-t1*t2
          u1=s2-em**2-t1
          u2=s1-em**2-t2
          b=((u2-nu)**2/X*t1+t1+4.*em**2)*
     &      ((u1-nu)**2/X*t2-t2)*
     &      3./16.*(t1-nu)**2/X*(-t2)*RM**2+
     &      ((u2-nu)**2/X*t1-t1)*
     &      ((u1-nu)**2/X*t2+t2+4.*em**2)*
     &      3./16.*(t2-nu)**2/X*(-t1)*RM**2
          cophi=(-2.*s+u1+u2-nu+4*em**2+
     &      nu*(u2-nu)*(u1-nu)/x)/
     &      sqrt(t1*t2*((u2-nu)**2/X-1.+4.*em**2/t1)*
     &      ((u1-nu)**2/X-1.+4.*em**2/t2))
          B=B-4.*abs((u2-nu)*(u1-nu)/X)*
     &      sqrt(
     &      abs(((u2-nu)**2/X*t1-t1+4.*em**2))*
     &      abs(((u1-nu)**2/X*t2-t2+4.*em**2))*t1*t2)*
     &      cophi*3./32.*(t1-nu)*(t2-nu)/X*sqrt(t1*t2)*RM**2
        else if(ir.eq.7)then        
c f1 polarization =0 sigma_TT
C sigma_TT=8pi^2(Gamma_{2gamma}/RM) A FAC2 delta(rm**2-rmas**2)
c tau_TT=-2sigma_TT        
c        A=3./2.*(t2-t1)**2/X*nu**2/sqrt(X)/RM**2
c b=(2rho_1^{++}2rho_2^{++} t1t2 A +
c    2|rho_1^{+-}rho_2^{+-}| t1t2 cos(2phi) (-2A)) \sqrt(X) RM**2/16        
          nu=(w**2-t1-t2)/2.
          X=nu**2-t1*t2
          u1=s2-em**2-t1
          u2=s1-em**2-t2
          b=((u2-nu)**2/X*t1+t1+4.*em**2)*
     &      ((u1-nu)**2/X*t2+t2+4.*em**2)*
     &      3./32.*(t2-t1)**2/X*nu**2
          cophi=(-2.*s+u1+u2-nu+4*em**2+
     &      nu*(u2-nu)*(u1-nu)/x)/
     &      sqrt(t1*t2*((u2-nu)**2/X-1.+4.*em**2/t1)*
     &      ((u1-nu)**2/X-1.+4.*em**2/t2))
          cophi=cophi**2*2.-1.
          B=B-
     &      abs(((u2-nu)**2/X*t1-t1+4.*em**2))*
     &      abs(((u1-nu)**2/X*t2-t2+4.*em**2))*
     &      cophi*3./32.*(t2-t1)**2/X*nu**2
        else
          print *,'not correct IR=',IR
          stop
        endif  

        
        FAC=1./(4.*PI)/E**4 *B/T1/T2
        FAC1=(R0*EM)**2*1.E33
        If(KVDM.ne.0)then
          FAC2=GGRFVDM(T1,T2)
        else
          FAC2=1.
        endif
        if(RAD.eq.1) then
          ll=dlog(-t2/em**2)-1.
          beta=alpha/pi*ll
          FAC3=(1.+
     &      alpha/pi*(0.75*ll-0.25))*
     &      beta/beta0*dexp(beta0*dlog(rmax))*
     &      dexp((beta-beta0)*dlog(rr))/
     &      (1.-ggpol(-t2))**2/(1.-ggpol(-t1))**2
        else
          FAC3=1.
        endif
        SS=FAC1*FAC*FAC2*FAC3*(RG*1.E-6)/RM**3*WEI*(RM/RMAS)**3 ! cr.sec. in nb
c        SS=FAC1*FAC*FAC2*FAC3*(RG*1.E-6)/RM**3*WEI
      endif
      SUM=SUM+SS
      ES=ES+SS**2
      IF(SS.GT.FM) then
        FM=SS      !  define maximum weight
      endif  
      IF(SS.GT.FMAX) NGT=NGT+1
      CALL GRNDM(RND,1)
      IF(FMAX*RND.GT.SS) GOTO 3

      Nact=Nact+1
c
c   Generation ee-> eeR
c
c   Calculation of lab.mom-s of scattered e+, e- and produced resonance
c      
      CALL GRNDM(RND,1)
      FI=2.*PI*RND
      CALL GGRLMOM(S,S1,S2,T1,T2,RM,FI)      
      NPART=3

      if(rr.gt.rmin) then
c   Lorenz transformation of mom-s of scattered e+, e- and produced
c   resonance into lab. system
        gamma=(2.d0-rr)/(2.d0*dsqrt(1.d0-rr))
c 
        ptot=DSQRT(PPART(1,1)**2+PPART(2,1)**2+PPART(3,1)**2)
        C(1)=PPART(1,1)/ptot
        C(2)=PPART(2,1)/ptot
        C(3)=PPART(3,1)/ptot
        CALL GGRLOR(GG,gamma,C,PPART(4,1),em)
        ptot=dsqrt(PPART(4,1)**2-em**2)
        PPART(1,1)=C(1)*ptot
        PPART(2,1)=C(2)*ptot
        PPART(3,1)=C(3)*ptot
c        
        ptot=DSQRT(PPART(1,2)**2+PPART(2,2)**2+PPART(3,2)**2)
        C(1)=PPART(1,2)/ptot
        C(2)=PPART(2,2)/ptot
        C(3)=PPART(3,2)/ptot
        CALL GGRLOR(GG,gamma,C,PPART(4,2),em)
        ptot=dsqrt(PPART(4,2)**2-em**2)
        PPART(1,2)=C(1)*ptot
        PPART(2,2)=C(2)*ptot
        PPART(3,2)=C(3)*ptot
c      
        ptot=DSQRT(PPART(1,3)**2+PPART(2,3)**2+PPART(3,3)**2)
        C(1)=PPART(1,3)/ptot
        C(2)=PPART(2,3)/ptot
        C(3)=PPART(3,3)/ptot
        CALL GGRLOR(GG,gamma,C,PPART(4,3),rm)
        ptot=dsqrt(PPART(4,3)**2-rm**2)       
        PPART(1,3)=C(1)*ptot
        PPART(2,3)=C(2)*ptot
        PPART(3,3)=C(3)*ptot
c photon momentum
        NPART=4
        type(4)=22
        mother(4)=0
        PPART(4,4)=eb0*rr
        PPART(1,4)=0.d0
        PPART(2,4)=0.d0
        PPART(3,4)=-PPART(4,4)
      endif            
c FSR generation
      if(RAD.eq.1)call ggfsr(t2)
c      
      IF(IR.EQ.1) THEN
        TYPE(3)=111
        if(imode.ne.0)CALL GGRDEC1(3)      ! pi0 decay
      else  IF(IR.EQ.2) THEN
        TYPE(3)=221
        if(imode.ne.0)CALL GGRETD(3,IMODE) ! eta decay
      else IF(IR.EQ.3) then
        TYPE(3)=331
        if(imode.ne.0)CALL GGRET1D(3,IMODE) ! eta' decay
      else IF(IR.EQ.4) then
        TYPE(3)=441
        IF(IMODE.NE.0) CALL GGRETCD(3,IMODE)
      else IF(IR.EQ.5) then
        TYPE(3)=551
      else IF(IR.EQ.6.or.IR.EQ.7) then
        TYPE(3)=20223
        IF(IMODE.eq.1) then
          CALL GGRf1M1(3,1)
        else if(IMODE.eq.2)then
          CALL GGRf1M2(3,1)
        else if(IMODE.eq.3)then
          CALL GGRf1M3(3,1)
        endif
      endif  
c      ptl(1)=0.
c      ptl(2)=0.
c      ptl(3)=0.
c      ptl(4)=0.
c      do i=1,npart
c        print *,i,type(i),mother(i),PPART(4,i),PPART(1,i),
c     &    PPART(2,i),PPART(3,i),
c     &    sqrt(PPART(4,i)**2-PPART(1,i)**2-PPART(2,i)**2-PPART(3,i)**2)
c     &    (PPART(4,i)**2-PPART(1,i)**2-PPART(2,i)**2-PPART(3,i)**2)
c        if(type(i).ne.20223)then
c          ptl(1)=ptl(1)+PPART(1,i)
c          ptl(2)=ptl(2)+PPART(2,i)
c          ptl(3)=ptl(3)+PPART(3,i)
c          ptl(4)=ptl(4)+PPART(4,i)
c        endif
c      enddo
c      print *,'****** ',ptl
c      print 98,t1,t2,rm,PPART(1,3)**2+PPART(2,3)**2
 98   format(4f15.5)
      RETURN
      END

c
      SUBROUTINE GGRINV(S,RM,T1,T2,S1,S2,WEI,Q2maxp)
      Implicit None
      real*8 S,RM,T1,T2,S1,S2,WEI,Q2maxp
      REAL*8 alpha,PI,EM,mpi0,mpi,meta,metap,mks,mkc,
     &  mrho,mjpsi,mups
      COMMON /GGRCON/alpha,PI,EM,mpi0,mpi,meta,metap,mks,mkc,
     &  mrho,mjpsi,mups  
      REAL*8 RND
      REAL*8 lambda,x,y,z
      real*8 m,w,m2,w2,beta
      real*8 wm2,t2min,t2max,dt2
      real*8 y2,q,a1,b1,del1,t1min,t1max,dt1
      real*8 y1,nupKW,ds1,Xs1
      real*8 ffa,ffb,ffdel
c      real*8 u1,u2,pht,cpht,nu,twoKW
c
      lambda(x,y,z) = (x - y - z)**2  - 4d0*y*z
c
      w     = rm
      m     = em
      m2    = m*m
      W2    = W*W
      beta  = dsqrt(1d0 - 4d0*m2/s)
c t2 generation    
      Wm2   = (W + 2d0*m)**2
      t2min=(4d0*m2-s+w2+2.d0*m*w-beta*dsqrt((s-w2)*(s-wm2)))/2d0
      t2max=m2*w2*wm2/s/t2min
C tagged electron Q2
      t2max=min(t2max,Q2maxp)
      if(t2min.ge.t2max)then
        wei=0.
        return
      endif
      dt2   = dlog(t2max/t2min)
      CALL GRNDM(RND,1)
      t2    = t2min*dexp(RND*dt2)
c t1 generation
      y2    = dsqrt(1d0 - 4d0*m2/t2)
      Q     = 4d0*(s+t2-4d0*m2)/(1d0+beta*y2) - t2 - W2
      a1    = 2d0*(Q+t2+2d0*m2+W2)
      b1    = Q**2-W2**2+2d0*W2*t2-t2**2-8d0*m2*t2-8d0*m2*W2
      Del1  = (Q+t2-W2+4d0*m*W)*(Q+t2-W2-4d0*m*W)*
     &  (Q**2-2d0*Q*t2+2d0*Q*W2+t2**2+W2**2-16d0*m2*t2-2d0*W2*t2)
      t1min = - (b1 + dsqrt(Del1))/2d0/a1
      t1max = 4d0*m2*(W2-t2)**2/(a1*t1min)
C untagged electron Q2 
      t1min = max(t1min,-0.6)
c      t1min = max(t1min,-1.0)
      if(t1min.ge.t1max)then
        wei=0.
        return
      endif  
      dt1   = dlog(t1max/t1min)
      CALL GRNDM(RND,1)
      t1  = t1min*dexp(RND*dt1)
c s1 generation
      y1    = dsqrt(1d0 - 4d0*m2/t1)     
      nupKW = (W2-t1-t2)/2d0 + dsqrt(lambda(W2,t1,t2))/2d0
      ds1   = dlog( s*(1d0+beta)**2/(nupKW*(1d0+y1)*(1d0+y2)) )
      CALL GRNDM(RND,1)
      Xs1   = nupKW*(1d0+y1)*dexp(ds1*RND)
      s1    = Xs1/2d0 + m2 + t2 + 2d0*m2*t2/Xs1
c s2 generation
      ffa = -2.D0*m**2*s1-2.D0*m**2*t2+m**4+s1**2-2.D0*t2*s1+t2**2
      ffb = -2.D0*m**6+4.D0*m**4*s1-4.D0*m**4*w**2+2.D0*S*t2*w**2-
     &  2.D0*m**2*s1**2+2.D0*m**2*S*w**2+2.D0*t1*S*s1-2.D0*m**2*t2*S+
     &  8.D0*t2*m**4-2.D0*m**2*t2**2-2.D0*t2*m**2*t1+4.D0*m**2*s1*w**2-
     &  2.D0*t2**2*S+2.D0*t1*m**4+2.D0*t1*t2*S+2.D0*S*t2*s1+
     &  2.D0*t1*t2*s1-2.D0*S*s1*w**2-2.D0*t1*s1**2-2.D0*S*m**2*t1
      ffdel = 16.D0*(m**2*t2**2+m**2*w**4-t1*s1*w**2-
     &  2.D0*m**2*t2*w**2+t1*t2*w**2-t1*m**2*w**2-2.D0*t1*m**2*s1-
     &  t1*t2*s1+t1**2*s1+t1*s1**2-t2*m**2*t1+t1*m**4)*(m**2*s1**2-
     &  2.D0*m**4*s1-S*t2*s1-3.D0*m**2*t2*S+t2*S**2+m**6+t2**2*S)
      If(ffdel.lt.0d0) then
        wei = 0d0
      else  
        CALL GRNDM(RND,1)
        s2 = (-ffb - dsqrt(ffdel)*dcos(RND*Pi))/(2d0*ffa)
        wei = dt1*dt2*ds1*4.d0*pi
      endif

      End
c
      REAL*8 FUNCTION GGRD4(S,S1,S2,T1,T2,W)
c
c   D4 - Granm determinant for reaction ee -> ee + R (with mass W)
c
      IMPLICIT none
c#include "AfkQed/ggrps.inc"
      include 'ggrps.inc'
      REAL*8 S,S1,S2,T1,T2,W
      REAL*8 EM2,U1,U2,W2,X,H
c
      EM2=EM**2
      U1=S2-T1-EM2
      U2=S1-T2-EM2
      W2=W**2
      X=0.5*(W2-T1-T2)
      H=-4.*T1*T2*(S-U1)*(S-U2)+(U1*U2+T1*T2-2.*S*X)**2
     &  +4.*EM2*(4.*S*(T1*T2-X**2)-T1*(U2+T2)**2-T2*(U1+T1)**2
     &  +2.*X*(U1+T1)*(U2+T2))
      GGRD4=H/16.
      RETURN
      END
c
      REAL*8 FUNCTION GGRGF(X,Y,Z,U,V,W)
c
c   Kinematical G-function
c
      IMPLICIT none
      REAL*8 X,Y,Z,U,V,W,GFUN
c
      GFUN=X**2*Y+X*Y**2+Z**2*U+Z*U**2+V**2*W+V*W**2+X*Z*W
     *     +X*U*V+Y*Z*V+Y*U*W-X*Y*(Z+U+V+W)
     *     -Z*U*(X+Y+V+W)-V*W*(X+Y+Z+U)
      GGRGF=GFUN
      RETURN
      END
c
      Real*8 function GGRFVDM(T1,T2)
c
c  VDM factor
c
c  T  - momentum transfer in (GeV/c)**2
c  AM - vector meson mass in GeV/c**2
c
      IMPLICIT none
c#include "AfkQed/ggrps.inc"      
      include 'ggrps.inc'
      real*8 T1,T2
      real*8 AM

      if(IR.le.3)then
        AM=Mrho
      else if(IR.eq.4)then
        AM=Mjpsi
      else if(IR.eq.5)then
        AM=Mups
      else
        AM=Mrho
      endif

c      GGRFVDM=1./(1.-T1/AM**2)**2/
c     &  (1.-T2/AM**2)**2
      GGRFVDM=1./(1.-(T1+t2)/AM**2)**2
c
      Return
      End
c
      SUBROUTINE GGRLMOM(S,S1,S2,T1,T2,W,FI)
c*********************************************************************
c                                                                    *
c  Calculation of lab. moments of produced e-, e+, resonance         *
c                                                                    *
c  Angle Teta - from direction of electrons (in +z axis)             *
c                                                                    *
c   Parameters of produced e-,e+ resonance in laboratory:            *
c   PPART(1-3,i) - Px,Py,Pz (GeV/c) of e- (i=1), e+(i=2), res.(i=3)  *
c   PPART(4,i)   - their total energy (GeV)                          *
c                                                                    *
c*********************************************************************
      IMPLICIT NONE
c#include "AfkQed/ggrps.inc"
      include 'ggrps.inc'
      REAL*8 S,S1,S2,T1,T2,W,FI
      REAL*8 E,EM2,PB,
     &  EPL,PPL,CTPL,STPL,
     &  EMI,PMI,CTMI,STMI,
     &  EW,PPW,CTW,STW,
     &  PTE,PTP,PTW,CFPL,CFMI,SIG,SFPL,SFMI,FIPL,FIMI  
      REAL*8 RND
      save
c
      EM2=EM**2
      E=DSQRT(S)/2.
      PB=DSQRT((E+EM)*(E-EM))
c
      EPL=(S-S1+EM2)/4./E             !  e+
      PPL=DSQRT((EPL-EM)*(EPL+EM))
      CTPL=-(E*EPL-EM2+0.5*T2)/PB/PPL
      IF(CTPL.le.-1.D0)then
        CTPL=-1.D0
        STPL=0.D0
      else if(CTPL.ge.1.D0)then
        CTPL=1.D0
        STPL=0.D0
      else
        STPL=sqrt((1.-CTPL)*(1.+CTPL))
      endif  
c
      EMI=(S-S2+EM2)/4./E             !  e-
      PMI=DSQRT((EMI-EM)*(EMI+EM))
      CTMI=(E*EMI-EM2+0.5*T1)/PB/PMI
      IF(CTMI.le.-1.D0)then
        CTMI=-1.D0
        STMI=0.D0
      else if(CTMI.ge.1.D0)then
        CTMI=1.D0
        STMI=0.D0
      else
        STMI=sqrt((1.-CTMI)*(1.+CTMI))
      endif  
c
      EW=(S1+S2-2*EM2)/4./E           !  resonance
      PPW=DSQRT((EW-W)*(EW+W))
      CTW=(E*EW-0.5*(S1+T1-T2-EM2))/PB/PPW
      IF(CTW.le.-1.D0)then
        CTW=-1.D0
        STW=0.D0
      else if(CTW.ge.1.D0)then
        CTW=1.D0
        STW=0.D0
      else
        STW=sqrt((1.-CTW)*(1.+CTW))
      endif        
c
      PTE=PMI*STMI
      PTP=PPL*STPL
      PTW=PPW*STW
c
      CFPL=(PTE**2-PTP**2-PTW**2)/2./PTP/PTW
      CFMI=(PTP**2-PTE**2-PTW**2)/2./PTE/PTW
      SIG=1.
c
      CALL GRNDM(RND,1)
c
      IF(RND.GT.0.5) SIG=-1.
      if(1.-CFPL**2.le.0.d0) then
        SFPL=0.d0
      else
        SFPL=SIG*DSQRT(1.-CFPL**2)
      endif
      SFMI=-PTP/PTE*SFPL
      FIPL=DATAN2(SFPL,CFPL)+FI
      FIMI=DATAN2(SFMI,CFMI)+FI      
c   Scattered electron
      PPART(3,1)=PMI*CTMI
      PPART(2,1)=PTE*DSIN(FIMI)
      PPART(1,1)=PTE*DCOS(FIMI)
      PPART(4,1)=EMI
      MPART(1)=EM
      TYPE(1)=11
      MOTHER(1)=0
c   Scattered positron
      PPART(3,2)=PPL*CTPL
      PPART(2,2)=PTP*DSIN(FIPL)
      PPART(1,2)=PTP*DCOS(FIPL)
      PPART(4,2)=EPL
      MPART(2)=EM
      TYPE(2)=-11
      MOTHER(2)=0
c   Produced resonance
      PPART(3,3)=PPW*CTW
      PPART(2,3)=PTW*DSIN(FI)
      PPART(1,3)=PTW*DCOS(FI)
      PPART(4,3)=EW
      MPART(3)=W
      TYPE(3)=0
      MOTHER(3)=0
      RETURN
      END
c
c=======================================================================
c
      Subroutine GGRETCD(IND,IMOD)
c*********************************************************************
c                                                                    *
c  Simulation of ETAC decay                                          *
c                                                                    *
c  Input: IND - position of etac in /GGREV/                          *
c                                                                    *
c  IMOD =1: etac -> KS K+ pi-                                        *
c                                                                    *
c*********************************************************************
      Implicit NONE
c#include "AfkQed/ggrps.inc"     
      include 'ggrps.inc'
      Integer IND,IMOD
      Real *8 CC,R,P
      Integer NP
      COMMON /GGRARI/  CC(4),R(20)
      COMMON /GGRPRT/  P(4,20),NP
      Integer i,ty(3),IND0
      Real*8 FMAT,FMAJ,RND,F,PMASS
      Real*8 GAM,PREST,GG(3),C(3),ENER,PTOT
c
      IND0=IND
      FMAT=1.                              ! used to avoid warning
c
      IF(IMOD.EQ.1)THEN    !   etac -> KS K+pi- phase space

        DO I=1,3
          CC(I)=0.
        END DO

        CC(4)=MPART(IND)    !  inv. mass of system of particles
        NP=3
        R(1)=mks
        R(2)=mkc
        R(3)=mpi
        ty(1)=310
        CALL GRNDM(RND,1)
        if(rnd.gt.0.5)then
          ty(2)=321
          ty(3)=-211
        else
          ty(2)=-321
          ty(3)=211
        endif
        FMAT=1.   !  uniform distribution
        FMAJ=15.
 1      CALL GGRSPC(F)
        CALL GRNDM(RND,1)
        IF(FMAJ*RND.GT.F*FMAT ) GOTO 1
        IF(F*FMAT.GT.FMAJ) PRINT 920, F*FMAT,FMAJ
 920    Format(' GGRETCD: F,Fmaj=',2(1pE10.2))

        GAM=PPART(4,ind)/MPART(ind) ! Lorentz-fact. of resonance in lab.
        PREST=
     &    DSQRT(PPART(1,ind)**2+PPART(2,ind)**2+PPART(3,ind)**2)
        GG(1)=PPART(1,ind)/PREST ! direction cosines in lab.
        GG(2)=PPART(2,ind)/PREST
        GG(3)=PPART(3,ind)/PREST

        DO I=1,3
          npart=npart+1
          MOTHER(npart)=IND0
          TYPE(npart)=ty(i)
          MPART(npart)=R(I)
          PMASS=R(I)
          C(1)=P(1,I)
          C(2)=P(2,I)
          C(3)=P(3,I)
          ENER=P(4,I)

          CALL GGRLOR(GG,GAM,C,ENER,PMASS)

          PTOT=DSQRT(ENER**2-R(I)**2)
          PPART(4,npart)=ENER
          PPART(1,npart)=PTOT*C(1)
          PPART(2,npart)=PTOT*C(2)
          PPART(3,npart)=PTOT*C(3)
        END DO
      else
        return
      endif
      End
c
c=======================================================================
c
      Subroutine GGRF1M1(IND,Ireg)
c*********************************************************************
c  Simulation of the f1(1285) -> a0(980)pi -> eta  pi+ pi-
c  f1 polarization: M=1 state for ir=6        
c                   M=0 state for ir=7        
c                                                                    *
c  Input: IND - position of f1 in /GGREV/                          *
c*********************************************************************
      Implicit NONE
c#include "AfkQed/ggrps.inc"     
      include 'ggrps.inc'
      Integer IND,ireg
      Real*8 GAM,PREST,GG(3),p(4),ptot
      Real*8 c(3),sint,cost,sinp,cosp
      Real*8 RND,f,c_t,s_f,c_f,s_t,fi,c1(3)
      real*8 min,max,am,pet0,peta,bw
      real*8 E
      real*8 GAMa0,ga0(3),p1(3),p2(3),p3(3),e1,e2,e3
      real*8 qq,vm0,vm1
      complex*16 bw1,bw2,v(3)
      common /f1par/p1,p2,p3,e1,e2,e3
c
      if(ireg.eq.0)then
c f1 mass
 1      min=meta+2.*mpi
        max=rmas+15.*rwid
        CALL GRNDM(RND,1)
        rm=(max-min)*rnd+min
c a0 angles
        CALL GRNDM(RND,1)
        C_T=-1.+2.*RND
        S_T=DSQRT(DABS(1.D0-C_T**2))
        CALL GRNDM(RND,1)
        FI=2.*PI*RND
        S_F=DSIN(FI)
        C_F=DCOS(FI)
        ga0(1)=S_T*C_F
        ga0(2)=S_T*S_F
        ga0(3)=C_T
c a0 mass
        min=meta+mpi
        max=rm-mpi
        CALL GRNDM(RND,1)
        am=(max-min)*rnd+min
c a0 energy and gamma factor
        E=(rm**2+am**2-mpi**2)/2./rm
        GAMa0=E/am
c recoil pi
        E1=(rm**2-am**2+mpi**2)/2./rm
        p1(1)=-ga0(1)
        p1(2)=-ga0(2)
        p1(3)=-ga0(3)
c pi from a0 decay
        E2=(am**2+mpi**2-meta**2)/2./am
        CALL GRNDM(RND,1)
        C_T=-1.+2.*RND
        S_T=DSQRT(DABS(1.D0-C_T**2))
        CALL GRNDM(RND,1)
        FI=2.*PI*RND
        S_F=DSIN(FI)
        C_F=DCOS(FI)
        p2(1)=S_T*C_F
        p2(2)=S_T*S_F
        p2(3)=C_T
        CALL GGRLOR(Ga0,GAMa0,p2,E2,mpi)
c  eta
        E3=(am**2-mpi**2+meta**2)/2./am
        p3(1)=-S_T*C_F
        p3(2)=-S_T*S_F
        p3(3)=-C_T
        CALL GGRLOR(Ga0,GAMa0,p3,E3,meta)
c decay rate
        qq=rm**2+mpi**2-2.*rm*e1
        bw1=ma0**2/cmplx(qq-ma0**2, ma0*wa0)*sqrt(e1**2-mpi**2)/rm
        qq=rm**2+mpi**2-2.*rm*e2
        bw2=ma0**2/cmplx(qq-ma0**2, ma0*wa0)*sqrt(e2**2-mpi**2)/rm
        V(1)=bw1*p1(1)+bw2*p2(1)
        V(2)=bw1*p1(2)+bw2*p2(2)
        V(3)=bw1*p1(3)+bw2*p2(3)
        vm0=abs(V(3))**2
        vm1=abs(bw1)**2+abs(bw2)**2+2.*real(bw1*conjg(bw2))*
     &    (p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3))-vm0
        f=(rm-2.*mpi-meta)/rmas*
     &    sqrt((am**2-(meta+mpi)**2)*(am**2-(meta-mpi)**2))/2./am/ma0*
     &    sqrt((rm**2-(am+mpi)**2)*(rm**2-(am-mpi)**2))/2./rm/rmas*
     &    (rmas*rwid)**2/((rm**2-rmas**2)**2+(rmas*rwid)**2)

        if(ir.eq.6)then
          f=vm1*f
        else
          f=vm0*f
        endif
        CALL GRNDM(RND,1)
        if(f.gt.0.6)print *, 'f=',f, ' > fmax=0.6'
        if(rnd*0.6.gt.f)goto 1
      else
c        print *,rm,sqrt(rm**2+mpi**2-2.*rm*e1),p1(3),
c     &    sqrt(rm**2+mpi**2-2.*rm*e2),p2(3)
        
c find the direction of gamma gamma colision axis in the f1 rest frame
        GAM=PPART(4,ind)/MPART(ind) ! Lorentz-fact. of resonance in lab.
        PREST=
     &    DSQRT(PPART(1,ind)**2+PPART(2,ind)**2+PPART(3,ind)**2)
        GG(1)=-PPART(1,ind)/PREST ! direction cosines in lab.
        GG(2)=-PPART(2,ind)/PREST
        GG(3)=-PPART(3,ind)/PREST

        P(1)=-ppart(1,1)
        P(2)=-ppart(2,1)
        p(3)=sqrt(eb**2-em**2)-ppart(3,1)
        P(4)=eb-ppart(4,1)

        call GGRLOR1(GG,GAM,P)
        
        ptot=sqrt(p(1)**2+p(2)**2+p(3)**2)
        C(1)=P(1)/ptot
        c(2)=P(2)/ptot
        c(3)=p(3)/ptot
        cost=c(3)
        sint=1.-cost**2
        if(sint.gt.0.)then
          sint=sqrt(sint)
        else
          sint=0.
        endif
        if(sint.gt.1.d-4)then
          cosp=C(1)/sint
          sinp=C(2)/sint
        else
          cosp=1.
          sinp=0.
        endif
        GG(1)=-gg(1)! direction cosines of f1 in lab.
        GG(2)=-gg(2)
        GG(3)=-gg(3)

c pi-
        c1(1)=cost*p1(1)+sint*p1(3)
        c1(2)=p1(2)
        c1(3)=-sint*p1(1)+cost*p1(3)
        p1(1)=cosp*c1(1)-sinp*c1(2)
        p1(2)=sinp*c1(1)+cosp*c1(2)
        p1(3)=c1(3)
        CALL GGRLOR(GG,GAM,p1,E1,mpi)
        npart=npart+1
        PPART(4,npart)=E1
        PTOT=DSQRT(E1**2-mpi**2)
        PPART(1,npart)=PTOT*p1(1)
        PPART(2,npart)=PTOT*p1(2)
        PPART(3,npart)=PTOT*p1(3)
        MPART(npart)=mpi
        TYPE(npart)=-211
        MOTHER(npart)=ind
c pi+
        c1(1)=cost*p2(1)+sint*p2(3)
        c1(2)=p2(2)
        c1(3)=-sint*p2(1)+cost*p2(3)
        p2(1)=cosp*c1(1)-sinp*c1(2)
        p2(2)=sinp*c1(1)+cosp*c1(2)
        p2(3)=c1(3)
        CALL GGRLOR(GG,GAM,p2,E2,mpi)
        npart=npart+1
        PPART(4,npart)=E2                     !  its  parameters in lab
        PTOT=DSQRT(E2**2-mpi**2)
        PPART(1,npart)=PTOT*p2(1)
        PPART(2,npart)=PTOT*p2(2)
        PPART(3,npart)=PTOT*p2(3)
        MPART(npart)=mpi
        TYPE(npart)=211
        MOTHER(npart)=ind
c eta
        c1(1)=cost*p3(1)+sint*p3(3)
        c1(2)=p3(2)
        c1(3)=-sint*p3(1)+cost*p3(3)
        p3(1)=cosp*c1(1)-sinp*c1(2)
        p3(2)=sinp*c1(1)+cosp*c1(2)
        p3(3)=c1(3)
        CALL GGRLOR(GG,GAM,p3,E3,meta)
        npart=npart+1
        PPART(4,npart)=E3                     !  its  parameters in lab
        PTOT=DSQRT(E3**2-meta**2)
        PPART(1,npart)=PTOT*p3(1)
        PPART(2,npart)=PTOT*p3(2)
        PPART(3,npart)=PTOT*p3(3)
        MPART(npart)=meta
        TYPE(npart)=221
        MOTHER(npart)=ind

      endif
      RETURN
      END

      Subroutine GGRF1M3(IND,Ireg)
c*********************************************************************
c  Simulation of the f1(1285) -> eta f0(500) -> eta  pi+ pi-
c  f1 polarization: M=1 state for ir=6        
c                   M=0 state for ir=7        
c                                                                    *
c  Input: IND - position of f1 in /GGREV/                          *
c*********************************************************************
      Implicit NONE
c#include "AfkQed/ggrps.inc"     
      include 'ggrps.inc'
      Integer IND,ireg
      Real*8 GAM,PREST,GG(3),p(4),ptot
      Real*8 c(3),sint,cost,sinp,cosp
      Real*8 RND,f,c_t,s_f,c_f,s_t,fi,c1(3)
      real*8 min,max,am,pet0,peta,bw
      real*8 E
      real*8 GAMa0,ga0(3),p1(3),p2(3),p3(3),e1,e2,e3
      real*8 qq,vm0,vm1
      complex*16 bw1,bw2,v(3)
      common /f1par/p1,p2,p3,e1,e2,e3
c
      if(ireg.eq.0)then
c f1 mass
 1      min=meta+2.*mpi
        max=rmas+15.*rwid
        CALL GRNDM(RND,1)
        rm=(max-min)*rnd+min
c f0 angles
        CALL GRNDM(RND,1)
        C_T=-1.+2.*RND
        S_T=DSQRT(DABS(1.D0-C_T**2))
        CALL GRNDM(RND,1)
        FI=2.*PI*RND
        S_F=DSIN(FI)
        C_F=DCOS(FI)
        ga0(1)=S_T*C_F
        ga0(2)=S_T*S_F
        ga0(3)=C_T
c f0 mass
        min=mpi+mpi
        max=rm-meta
        CALL GRNDM(RND,1)
        am=(max-min)*rnd+min
c f0 energy and gamma factor
        E=(rm**2+am**2-meta**2)/2./rm
        GAMa0=E/am
c recoil eta
        E3 =(rm**2-am**2+meta**2)/2./rm
        p3(1)=-ga0(1)
        p3(2)=-ga0(2)
        p3(3)=-ga0(3)
c pi- from f0 decay
        E1=am/2.
        CALL GRNDM(RND,1)
        C_T=-1.+2.*RND
        S_T=DSQRT(DABS(1.D0-C_T**2))
        CALL GRNDM(RND,1)
        FI=2.*PI*RND
        S_F=DSIN(FI)
        C_F=DCOS(FI)
        p1(1)=S_T*C_F
        p1(2)=S_T*S_F
        p1(3)=C_T
        CALL GGRLOR(Ga0,GAMa0,p1,E1,mpi)
c  pi+ from f0 decay
        E2=am/2.
        p2(1)=-S_T*C_F
        p2(2)=-S_T*S_F
        p2(3)=-C_T
        CALL GGRLOR(Ga0,GAMa0,p2,E2,mpi)
c decay rate
        qq=rm**2+meta**2-2.*rm*e3
        bw1=mf0**2/cmplx(qq-mf0**2, mf0*wf0)*sqrt(e3**2-meta**2)/rm
        V(1)=bw1*p3(1)
        V(2)=bw1*p3(2)
        V(3)=bw1*p3(3)
        vm0=abs(V(3))**2
        vm1=abs(bw1)**2-vm0
        f=(rm-2.*mpi-meta)/rmas*
     &    sqrt(am**2-4.*mpi**2)/2./mf0*
     &    sqrt((rm**2-(am+meta)**2)*(rm**2-(am-meta)**2))/2./rm/rmas*
     &    (rmas*rwid)**2/((rm**2-rmas**2)**2+(rmas*rwid)**2)

        if(ir.eq.6)then
          f=vm1*f
        else
          f=vm0*f
        endif
        CALL GRNDM(RND,1)
        if(f.gt.0.01)print *, 'f=',f, ' > fmax=0.01'
        if(rnd*0.01.gt.f)goto 1
      else
c        print *,rm,sqrt(rm**2+mpi**2-2.*rm*e1),p1(3),
c     &    sqrt(rm**2+mpi**2-2.*rm*e2),p2(3)
c        print *,rm,sqrt(rm**2+meta**2-2.*rm*e3),p3(3)
c find the direction of gamma gamma colision axis in the f1 rest frame
        GAM=PPART(4,ind)/MPART(ind) ! Lorentz-fact. of resonance in lab.
        PREST=
     &    DSQRT(PPART(1,ind)**2+PPART(2,ind)**2+PPART(3,ind)**2)
        GG(1)=-PPART(1,ind)/PREST ! direction cosines in lab.
        GG(2)=-PPART(2,ind)/PREST
        GG(3)=-PPART(3,ind)/PREST

        P(1)=-ppart(1,1)
        P(2)=-ppart(2,1)
        p(3)=sqrt(eb**2-em**2)-ppart(3,1)
        P(4)=eb-ppart(4,1)

        call GGRLOR1(GG,GAM,P)
        
        ptot=sqrt(p(1)**2+p(2)**2+p(3)**2)
        C(1)=P(1)/ptot
        c(2)=P(2)/ptot
        c(3)=p(3)/ptot
        cost=c(3)
        sint=1.-cost**2
        if(sint.gt.0.)then
          sint=sqrt(sint)
        else
          sint=0.
        endif
        if(sint.gt.1.d-4)then
          cosp=C(1)/sint
          sinp=C(2)/sint
        else
          cosp=1.
          sinp=0.
        endif
        GG(1)=-gg(1)! direction cosines of f1 in lab.
        GG(2)=-gg(2)
        GG(3)=-gg(3)

c pi-
        c1(1)=cost*p1(1)+sint*p1(3)
        c1(2)=p1(2)
        c1(3)=-sint*p1(1)+cost*p1(3)
        p1(1)=cosp*c1(1)-sinp*c1(2)
        p1(2)=sinp*c1(1)+cosp*c1(2)
        p1(3)=c1(3)
        CALL GGRLOR(GG,GAM,p1,E1,mpi)
        npart=npart+1
        PPART(4,npart)=E1
        PTOT=DSQRT(E1**2-mpi**2)
        PPART(1,npart)=PTOT*p1(1)
        PPART(2,npart)=PTOT*p1(2)
        PPART(3,npart)=PTOT*p1(3)
        MPART(npart)=mpi
        TYPE(npart)=-211
        MOTHER(npart)=ind
c pi+
        c1(1)=cost*p2(1)+sint*p2(3)
        c1(2)=p2(2)
        c1(3)=-sint*p2(1)+cost*p2(3)
        p2(1)=cosp*c1(1)-sinp*c1(2)
        p2(2)=sinp*c1(1)+cosp*c1(2)
        p2(3)=c1(3)
        CALL GGRLOR(GG,GAM,p2,E2,mpi)
        npart=npart+1
        PPART(4,npart)=E2                     !  its  parameters in lab
        PTOT=DSQRT(E2**2-mpi**2)
        PPART(1,npart)=PTOT*p2(1)
        PPART(2,npart)=PTOT*p2(2)
        PPART(3,npart)=PTOT*p2(3)
        MPART(npart)=mpi
        TYPE(npart)=211
        MOTHER(npart)=ind
c eta
        c1(1)=cost*p3(1)+sint*p3(3)
        c1(2)=p3(2)
        c1(3)=-sint*p3(1)+cost*p3(3)
        p3(1)=cosp*c1(1)-sinp*c1(2)
        p3(2)=sinp*c1(1)+cosp*c1(2)
        p3(3)=c1(3)
        CALL GGRLOR(GG,GAM,p3,E3,meta)
        npart=npart+1
        PPART(4,npart)=E3                     !  its  parameters in lab
        PTOT=DSQRT(E3**2-meta**2)
        PPART(1,npart)=PTOT*p3(1)
        PPART(2,npart)=PTOT*p3(2)
        PPART(3,npart)=PTOT*p3(3)
        MPART(npart)=meta
        TYPE(npart)=221
        MOTHER(npart)=ind

      endif
      RETURN
      END

c=======================================================================
c
      Subroutine GGRF1M2(IND,ireg)
c*********************************************************************
c                                                                    *
c  Simulation of f1(1485)-> rho gamma ->  pi+ pi- gamma  decay
c  f1 polarization: ir=6 - m=1, ir=7 - m=0 
c                                                                    *
c  Input: IND - position of f1 in /GGREV/                          *
c*********************************************************************
      Implicit NONE
c#include "AfkQed/ggrps.inc"     
      include 'ggrps.inc'
      Integer IND,ireg
      Real*8 GAM,PREST,GG(3),p(4),ptot
      Real*8 c(3),sint,cost,sinp,cosp
      Real*8 RND,f,c_t,s_f,c_f,s_t,fi,c1(3)
      real*8 min,max,am,pet0,peta,bw
      real*8 E,SIGN
      real*8 GAMrho,grho(3),p1(3),p2(3),p3(3),e1,e2,e3
      real*8 a,de,dp(3),cothg,mgam,jj,jn
      common /f1parn/p1,p2,p3,e1,e2,e3,am,cothg
c
c
      if(ireg.eq.0)then
c f1 mass
 1      min=2.*mpi
        max=rmas+15.*rwid
        CALL GRNDM(RND,1)
        rm=(max-min)*rnd+min
c rho angles
        CALL GRNDM(RND,1)
        C_T=-1.+2.*RND
        S_T=DSQRT(DABS(1.D0-C_T**2))
        CALL GRNDM(RND,1)
        FI=2.*PI*RND
        S_F=DSIN(FI)
        C_F=DCOS(FI)
        grho(1)=S_T*C_F
        grho(2)=S_T*S_F
        grho(3)=C_T
c rho mass
        min=mpi+mpi
        max=rm-0.001
        CALL GRNDM(RND,1)
        am=(max-min)*rnd+min
c rho energy and gamma factor
        E=(rm**2+am**2)/2./rm
        GAMrho=E/am      
c photon energy momentum
        E1=(rm**2-am**2)/2./rm
        p1(1)=-grho(1)
        p1(2)=-grho(2)
        p1(3)=-grho(3)
        cothg=p1(3)
c pi- energy, momentum
        E2=am/2.
        CALL GRNDM(RND,1)
        C_T=-1.+2.*RND
        S_T=DSQRT(DABS(1.D0-C_T**2))
        CALL GRNDM(RND,1)
        FI=2.*PI*RND
        S_F=DSIN(FI)
        C_F=DCOS(FI)
        p2(1)=S_T*C_F
        p2(2)=S_T*S_F
        p2(3)=C_T
        CALL GGRLOR(Grho,GAMrho,p2,E2,mpi)
c  pi+
        E3=am/2.
        p3(1)=-S_T*C_F
        p3(2)=-S_T*S_F
        p3(3)=-C_T
        CALL GGRLOR(Grho,GAMrho,p3,E3,mpi)
c decay rate
        a=am/rm
        if(am.ge.mpi+mpi)then
          pet0=sqrt(mrho**2-4.*mpi**2)
          peta=sqrt(am**2-4.*mpi**2)
          bw=(mrho*wrho)**2/((am**2-mrho**2)**2+(mrho*wrho)**2*
     &      mrho**2/am**2*(peta/pet0)**3)
        else
          bw=0.
          goto 1
        endif
        dp(1)=sqrt(e2**2-mpi**2)*p2(1)-sqrt(e3**2-mpi**2)*p3(1)
        dp(2)=sqrt(e2**2-mpi**2)*p2(2)-sqrt(e3**2-mpi**2)*p3(2)
        dp(3)=sqrt(e2**2-mpi**2)*p2(3)-sqrt(e3**2-mpi**2)*p3(3)
        jn=dp(1)*p1(1)+dp(2)*p1(2)+dp(3)*p1(3)
        jj=dp(1)**2+dp(2)**2+dp(3)**2
      
        if(ir.eq.6)then
          f=jn**2+(1.-cothg**2)/8.*
     &      ((1.+a**2)**2*(1.-(1.-a**2)*G12)**2*jj-
     &      (1.-a**2)**2*(1.+(1.+a**2)*G12)**2*jn**2)-
     &      (1.+a**2)/2.*(1.-(1.-a**2)*G12)*jn*(jn-dp(3)*cothg)
        else
          f=jn**2+cothg**2/4.*
     &      ((1.+a**2)**2*(1.-(1.-a**2)*G12)**2*jj-
     &      (1.-a**2)**2*(1.+(1.+a**2)*G12)**2*jn**2)-
     &      (1.+a**2)*(1.-(1.-a**2)*G12)*jn*dp(3)*cothg
        endif
        f=f*bw*peta**3*(1.-a**2)**3/mrho**3*(rm-2.*mpi)/rmas*
     &    (rmas*rwid)**2/((rm**2-rmas**2)**2+(rmas*rwid)**2)

        CALL GRNDM(RND,1)
        if(f.gt.0.11)print *, 'f=',f, ' > fmax=0.11'
        if(rnd*0.11.gt.f)goto 1
        
      else
c find the direction of gamma gamma colision axis in the f1 rest frame
c        print *,rm,am,cothg

        GAM=PPART(4,ind)/MPART(ind) ! Lorentz-fact. of resonance in lab.
        PREST=
     &    DSQRT(PPART(1,ind)**2+PPART(2,ind)**2+PPART(3,ind)**2)
        GG(1)=-PPART(1,ind)/PREST ! direction cosines in lab.
        GG(2)=-PPART(2,ind)/PREST
        GG(3)=-PPART(3,ind)/PREST

        P(1)=-ppart(1,1)
        P(2)=-ppart(2,1)
        p(3)=sqrt(eb**2-em**2)-ppart(3,1)
        P(4)=eb-ppart(4,1)

        call GGRLOR1(GG,GAM,P)

        ptot=sqrt(p(1)**2+p(2)**2+p(3)**2)
        C(1)=P(1)/ptot
        c(2)=P(2)/ptot
        c(3)=p(3)/ptot
        cost=c(3)
        sint=1.-cost**2
        if(sint.gt.0.)then
          sint=sqrt(sint)
        else
          sint=0.
        endif
        if(sint.gt.1.d-4)then
          cosp=C(1)/sint
          sinp=C(2)/sint
        else
          cosp=1.
          sinp=0.
        endif
        GG(1)=-gg(1)! direction cosines of f1 in lab.
        GG(2)=-gg(2)
        GG(3)=-gg(3)
c photon
        c1(1)=p1(1)
        c1(2)=cost*p1(2)+sint*p1(3)
        c1(3)=-sint*p1(2)+cost*p1(3)
        p1(1)=cosp*c1(2)+sinp*c1(1)
        p1(2)=sinp*c1(2)-cosp*c1(1)
        p1(3)=c1(3)
        mgam=0.
        CALL GGRLOR(GG,GAM,p1,E1,mgam)
        npart=npart+1
        PPART(4,npart)=E1                 
        PTOT=E1
        PPART(1,npart)=PTOT*p1(1)
        PPART(2,npart)=PTOT*p1(2)
        PPART(3,npart)=PTOT*p1(3)
        MPART(npart)=0.
        TYPE(npart)=22
        MOTHER(npart)=ind
c pi-
        c1(1)=p2(1)
        c1(2)=cost*p2(2)+sint*p2(3)
        c1(3)=-sint*p2(2)+cost*p2(3)
        p2(1)=cosp*c1(2)+sinp*c1(1)
        p2(2)=sinp*c1(2)-cosp*c1(1)
        p2(3)=c1(3)
        
        CALL GGRLOR(GG,GAM,p2,E2,mpi)
        npart=npart+1
        PPART(4,npart)=E2                     !  its  parameters in lab
        PTOT=DSQRT(E2**2-mpi**2)
        PPART(1,npart)=PTOT*p2(1)
        PPART(2,npart)=PTOT*p2(2)
        PPART(3,npart)=PTOT*p2(3)
        MPART(npart)=mpi
        TYPE(npart)=-211
        MOTHER(npart)=ind
c pi+
        c1(1)=p3(1)
        c1(2)=cost*p3(2)+sint*p3(3)
        c1(3)=-sint*p3(2)+cost*p3(3)
        p3(1)=cosp*c1(2)+sinp*c1(1)
        p3(2)=sinp*c1(2)-cosp*c1(1)
        p3(3)=c1(3)

        CALL GGRLOR(GG,GAM,p3,E3,mpi)
        npart=npart+1
        PPART(4,npart)=E3                     !  its  parameters in lab
        PTOT=DSQRT(E3**2-mpi**2)
        PPART(1,npart)=PTOT*p3(1)
        PPART(2,npart)=PTOT*p3(2)
        PPART(3,npart)=PTOT*p3(3)
        MPART(npart)=mpi
        TYPE(npart)=211
        MOTHER(npart)=ind
      endif  
      return
      END


      SUBROUTINE GGRET1D(ind,IMOD)
c*********************************************************************
c                                                                    *
c  Simulation of eta' decays                                         *
c                                                                    *
c  IMODE=decay mode:
c             1  - eta' -> 2gamma,
c             2  - eta' -> pi+ pi- eta (eta->2gamma)
c             3  - eta' -> pi+ pi- eta (eta ->3pi0)
c             4  - eta' -> pi+ pi- eta (eta ->pi+ pi- pi0)
c             5  - eta' -> pi+ pi- eta (eta ->pi+ pi- gamma)
c             6  - eta' -> pi0 pi0 eta (eta ->2gamma)
c             7  - eta' -> pi0 pi0 eta (eta ->3pi0)
c             8  - eta' -> pi0 pi0 eta (eta ->pi+ pi- pi0)
c             9  - eta' -> pi0 pi0 eta (eta ->pi+ pi- gamma)
c            10  - eta' -> rho0 gamma (rho0 ->pi+ pi-)
c                                                                    *
c   In /GGREV/ this subr. adds (i=4-8):                              *
c   PPART(1-3,i) - moments (GeV/c) of particles from resonance decay *
c   PPART(4,i)   - their total energy (GeV)                          *
c   PPART(5,i)   - their GEANT types                                 *
c                                                                    *
c*********************************************************************
      IMPLICIT none
c#include "AfkQed/ggrps.inc"
      include 'ggrps.inc'
      integer ind,IMOD,nn
c   eta' -> 2gamma
      IF(IMOD.EQ.1)THEN
        CALL GGRDEC1(IND)
c   eta' -> pi +pi- eta
      else if(imod.ge.2.and.imod.le.5)then
        CALL GGRDEC4(ind,1)
        call GGRETD(npart,imod-1)
c   eta' -> pi0 pi0 eta
      else if(imod.ge.6.and.imod.le.9)then
        CALL GGRDEC4(ind,2)
        nn=npart
        CALL GGRDEC1(nn-2)
        CALL GGRDEC1(nn-1)
        call GGRETD(nn,imod-5)
c   eta' -> rh0 gamma (rho->pi+pi-)
      else if(imod.eq.10)then
        CALL GGRDEC5(ind)
      endif  
      RETURN
      END
c
      SUBROUTINE GGRDEC5(ind)
c*********************************************************************
c
c  Decay eta' -> rh0 gamma
c
c    Input:
c
c    PRES(4) - lab.4-momentum of eta'
c
c    Output:
c
c    P1(4) - lab.4-momentum of pi+
c    P2(4) - lab.4-momentum of pi-
c    P2(4) - lab.4-momentum of gamma
c
c*********************************************************************
      IMPLICIT none
c#include "AfkQed/ggrps.inc"
      include 'ggrps.inc'
      integer ind,np
      REAL*8 RHOM,RHOG
      parameter (RHOM=0.7693)
      parameter (RHOG=0.1502)
      REAL*8 CC,R,P,RND
      COMMON /GGRARI/  CC(4),R(20)
      COMMON /GGRPRT/  P(4,20),NP
      real*8 FSPC,e0,p0,fmaj,P1TOT,P2TOT,P2PI(4),EPIPI,AMPIPI
      real*8 GAM,PPIPI,PREST,GG(3),C(3),E
      real*8 P1A(4),PTOT,P3A(4),PPI,CT,ST2,FAC,AMELSQ,F      
      integer i
c
      R(1)=mpi   ! pi+ mass
      R(2)=mpi   ! pi- mass
      R(3)=0.    ! gamma mass
      E0=mrho/2.
      P0=DSQRT(E0**2-mpi**2)   ! mom. of pi for M(pipi)=M(rho)
      FMAJ=0.5
c
      DO I=1,3
        CC(I)=0.
      END DO
      CC(4)=MPART(ind) ! inv. mass of system of particles (eta')
      NP=3             ! particles from eta' decay
c
c   At first - decay in eta' system
c
1     CALL GGRSPC(FSPC)
c
      P1TOT=DSQRT(P(4,1)**2-mpi**2)
      P2TOT=DSQRT(P(4,2)**2-mpi**2)
c
      DO I=1,3
	P2PI(I)=P1TOT*P(I,1)+P2TOT*P(I,2)  ! 3-momentum of pi+pi- system
      END DO
      P2PI(4)=P(4,1)+P(4,2)                ! energy of pi+pi- system
c
      EPIPI=P2PI(4)
      AMPIPI=DSQRT(EPIPI**2-P2PI(1)**2-P2PI(2)**2-P2PI(3)**2) ! mass pi+pi-
c
c   Lorentz transformation for pi+, gamma to pipi-system
c
      GAM=EPIPI/AMPIPI      ! Lorentz-factor of pipi in eta'-system
      PPIPI=DSQRT(P2PI(1)**2+P2PI(2)**2+P2PI(3)**2)
c
      DO I=1,3
        GG(I)=P2PI(I)/PPIPI ! direction cosines of pipi in eta'-system
	C(I)=P(I,1)         ! dir. cos of pi+
      END DO
      E=P(4,1)              ! energy of pi+
c
      CALL GGRLOR(GG,GAM,C,E,mpi)
c
      P1A(4)=E                     !  pi+ parameters in pipi-system
      PTOT=DSQRT(E**2-mpi**2)
      P1A(1)=PTOT*C(1)
      P1A(2)=PTOT*C(2)
      P1A(3)=PTOT*C(3)
c
      E=P(4,3)   ! energy of photon
c
      DO I=1,3
	C(I)=P(I,3)
      END DO
c
      CALL GGRLOR(GG,GAM,C,E,0.D0)
c
      P3A(1)=E*C(1)   !  photon parameters in pipi-system
      P3A(2)=E*C(2)
      P3A(3)=E*C(3)
      P3A(4)=E
c
c   |M|**2 calculation
c
      PPI=DSQRT(P1A(1)**2+P1A(2)**2+P1A(3)**2)
      CT=(P1A(1)*P3A(1)+P1A(2)*P3A(2)+P1A(3)*P3A(3))/PPI/E
      ST2=DSQRT(1.-CT**2)
      GAM=RHOG*(PPI/P0)**3*2*P0**2/(PPI**2+P0**2)
      FAC=1./((AMPIPI**2-RHOM**2)**2+(RHOM*GAM)**2)
      AMELSQ=PPI**2*E**2*AMPIPI**2*ST2*FAC ! |M|**2 for eta' -> pi+ pi- gamma
c
      CALL GRNDM(RND,1)
c
      F=FSPC*AMELSQ
      IF(FMAJ*RND.GT.F) GOTO 1
      IF(F.GT.FMAJ) PRINT 920, F,FMAJ
c
c   Transformation of parameters of pi+ pi- gamma to lab.
c
      GAM=PPART(4,ind)/MPART(ind) ! Lorentz-factor of eta' resonance in lab.
      PREST=DSQRT(PPART(1,ind)**2+PPART(2,ind)**2+PPART(3,ind)**2)
      GG(1)=PPART(1,ind)/PREST ! direction cosines of eta' in lab.
      GG(2)=PPART(2,ind)/PREST
      GG(3)=PPART(3,ind)/PREST
c
c  pi+
c
      DO I=1,3
        C(I)=P(I,1)
      END DO
      E=P(4,1)
c
      CALL GGRLOR(GG,GAM,C,E,mpi)
c
      npart=npart+1
      PPART(4,npart)=E                     !  pi+ parameters in lab
      PTOT=DSQRT(E**2-mpi**2)
      PPART(1,npart)=PTOT*C(1)
      PPART(2,npart)=PTOT*C(2)
      PPART(3,npart)=PTOT*C(3)
      MPART(npart)=mpi
      TYPE(npart)=211
      MOTHER(npart)=ind
c
c  pi-
c
      DO I=1,3
        C(I)=P(I,2)
      END DO
      E=P(4,2)
c
      CALL GGRLOR(GG,GAM,C,E,mpi)
c
      npart=npart+1
      PPART(4,npart)=E                     !  pi- parameters in lab
      PTOT=DSQRT(E**2-mpi**2)
      PPART(1,npart)=PTOT*C(1)
      PPART(2,npart)=PTOT*C(2)
      PPART(3,npart)=PTOT*C(3)
      MPART(npart)=mpi
      TYPE(npart)=-211
      MOTHER(npart)=ind
c
c  photon
c
      DO I=1,3
        C(I)=P(I,3)
      END DO
      E=P(4,3)
c
      CALL GGRLOR(GG,GAM,C,E,0.D0)
c
      npart=npart+1
      PPART(4,npart)=E                     !  photon parameters in lab
      PTOT=DSQRT(E**2-mpi**2)
      PPART(1,npart)=PTOT*C(1)
      PPART(2,npart)=PTOT*C(2)
      PPART(3,npart)=PTOT*C(3)
      MPART(npart)=0.
      TYPE(npart)=22
      MOTHER(npart)=ind
c
920   Format(' GGRDEC5: F,Fmaj=',2(1pE10.2))
      RETURN
      END
c
      SUBROUTINE GGRDEC4(ind,IMOD)
c*********************************************************************
c
c   Decays eta' -> pi+ pi- eta (IMODE=1)
c                  pi0 pi0 eta (IMODE=2)
c
c*********************************************************************
      IMPLICIT none
c#include "AfkQed/ggrps.inc"
      include 'ggrps.inc'
      integer ind,imod,np
      REAL*8 CC,R,P
      COMMON /GGRARI/  CC(4),R(20)
      COMMON /GGRPRT/  P(4,20),NP
      integer ty(3)
      real*8 FSPC,fmaj,ENER1,ENER2,ENER3,y,AMELSQ,f,rnd
      real*8 GAM,PREST,GG(3),PMASS,ENER,C(3),PTOT
      integer i
c
      DO I=1,3
        CC(I)=0.
      END DO
      CC(4)=MPART(ind) ! inv. mass of system of particles (eta')
      NP=3     ! particles from eta' decay
c
      IF(IMOD.eq.1)then
        R(1)=mpi  
        R(2)=mpi
        R(3)=meta
        ty(1)=211
        ty(2)=-211
        ty(3)=221
        FMAJ=6.D-2
      else if(IMOD.eq.2)then
        R(1)=mpi0  
        R(2)=mpi0
        R(3)=meta
        ty(1)=111
        ty(2)=111
        ty(3)=221
        FMAJ=6.D-2
      else
        return
      endif  
c
1     CALL GGRSPC(FSPC)
c
c  |M|**2 for eta' -> eta pi pi (G.R.Kalbfleisch, PR D10(1974)916), PDG2000
c
      ENER1=P(4,1)-R(1)   ! Ekin for pi
      ENER2=P(4,2)-R(2)
      ENER3=P(4,3)-R(3)   ! Ekin for eta
c
      y=(2*ENER3-ENER1-ENER2)/(ENER1+ENER2+ENER3)
      AMELSQ=(1-0.058*y)**2
c
      CALL GRNDM(RND,1)
c
      F=FSPC*AMELSQ
      IF(FMAJ*RND.GT.F) GOTO 1
      IF(F.GT.FMAJ) PRINT 920, F,FMAJ
c
      GAM=PPART(4,ind)/MPART(ind) ! Lorentz-factor of eta' resonance in lab.
      PREST=DSQRT(PPART(1,ind)**2+PPART(2,ind)**2+PPART(3,ind)**2)
      GG(1)=PPART(1,ind)/PREST ! direction cosines of eta' in lab.
      GG(2)=PPART(2,ind)/PREST
      GG(3)=PPART(3,ind)/PREST
c
      DO I=1,3
        C(1)=P(1,I)
        C(2)=P(2,I)
        C(3)=P(3,I)
        PMASS=R(I)
        ENER=P(4,I)
c
        CALL GGRLOR(GG,GAM,C,ENER,PMASS)
c
        npart=npart+1
        PPART(4,npart)=ENER
        PTOT=DSQRT(ENER**2-PMASS**2)
        PPART(1,npart)=PTOT*C(1)
        PPART(2,npart)=PTOT*C(2)
        PPART(3,npart)=PTOT*C(3)
        MPART(npart)=PMASS
        TYPE(npart)=ty(i)
        MOTHER(npart)=ind
c
      enddo
920   Format(' GGRDEC4: F,Fmaj=',2(1pE10.2))
      RETURN
      END
c
      SUBROUTINE GGRDEC1(IND)
c-------------------------------------------------------------------
c  Isotropic decay of PS-resonance to 2 gamma in lab.system.
c---------------------------------------------------------------------
      IMPLICIT NONE
c#include "AfkQed/ggrps.inc"
      include 'ggrps.inc'
      integer ind
      REAL*8 PRES(4),PGAM1(4),PGAM2(4),GG(3),C(3)
      REAl*8 RND
      REAl*8 C_T,S_T,FI,S_F,C_F,EPI,GAM,PREST,E
      INTEGER i
c
      DO I=1,4
        PRES(I)=PPART(I,IND)
      END DO

      CALL GRNDM(RND,1)
      C_T=-1.+2.*RND
      S_T=DSQRT(DABS(1.D0-C_T**2))
      CALL GRNDM(RND,1)
      FI=2.*PI*RND
      S_F=DSIN(FI)
      C_F=DCOS(FI)
c
      EPI=PRES(4)
      GAM=EPI/MPART(IND) 
      PREST=DSQRT(PRES(1)**2+PRES(2)**2+PRES(3)**2)
      GG(1)=PRES(1)/PREST 
      GG(2)=PRES(2)/PREST
      GG(3)=PRES(3)/PREST
c
      C(3)=C_T
      C(1)=S_T*C_F
      C(2)=S_T*S_F
      E=MPART(IND)/2.
c
      CALL GGRLOR(GG,GAM,C,E,0.D0)
c
      PGAM1(4)=E
      PGAM1(1)=E*C(1)
      PGAM1(2)=E*C(2)
      PGAM1(3)=E*C(3)
c
      C(3)=-C_T
      C(1)=-S_T*C_F
      C(2)=-S_T*S_F
      E=MPART(IND)/2.
c
      CALL GGRLOR(GG,GAM,C,E,0.D0)
c
      PGAM2(4)=E
      PGAM2(1)=E*C(1)
      PGAM2(2)=E*C(2)
      PGAM2(3)=E*C(3)

      DO I=1,4
        PPART(I,NPART+1)=PGAM1(I)
        PPART(I,NPART+2)=PGAM2(I)
      END DO
      MPART(NPART+1)=0.d0
      MPART(NPART+2)=0.d0
      type(NPART+1)=22
      type(NPART+2)=22
      mother(NPART+1)=ind
      mother(NPART+2)=ind
      NPART=NPART+2
c
      RETURN
      END
c
      SUBROUTINE GGRETD(IND,IMOD)
c*********************************************************************
c                                                                    *
c  Simulation of ETA0 decay                                          *
c                                                                    *
c   Input:                                                           *
c                                                                    *
c  IMODE=1: eta -> 2gamma      (PDG: 39.2%)                          *
c        2: eta -> 3pi0        (PDG: 32.2%)                          *
c        3: eta -> pi+pi-pi0   (PDG: 21.3%)                          *
c        4: eta -> pi+pi-gamma (PDG: 4.8%)                           *
c*********************************************************************
      IMPLICIT NONE
c#include "AfkQed/ggrps.inc"
      include 'ggrps.inc'
      integer IND,IMOD
      REAL*8 CC,R,P
      integer NP
      COMMON /GGRARI/ CC(4),R(20)
      COMMON /GGRPRT/ P(4,20),NP
      integer i,i1,i2,ty(3)
      real*8 FMAT,FMAJ,RND,F,Q,T0,Y
      real*8 P1MOM,P2MOM,PGMOM,P1P2,P1PG,P2PG
      real*8 GAM,PREST,GG(3),C(3),ENER,PTOT

      FMAT=1.
c
c   eta -> 2gamma
      IF(IMOD.EQ.1)THEN
        CALL GGRDEC1(IND)
c   eta -> 3 particles
      else if(IMOD.le.4)THEN
        DO I=1,3
          CC(I)=0.
        END DO
        CC(4)=MPART(IND)    !  inv. mass of system of particles
        NP=3
c
        IF(IMOD.eq.2)then
c   eta -> 3pi0             (PDG: 32.2 %)
          R(1)=mpi0
          R(2)=mpi0
          R(3)=mpi0
          ty(1)=111
          ty(2)=111
          ty(3)=111
          FMAT=1.   !  uniform distribution
          FMAJ=3.
        else IF(IMOD.eq.3)then
c  eta -> pi+pi-pi0         (PDG: 21.3%)
          R(1)=mpi
          R(2)=mpi
          R(3)=mpi0
          ty(1)=211
          ty(2)=-211
          ty(3)=111
          FMAJ=1.
        else IF(IMOD.eq.4)then
c  eta -> pi+pi-gamma       (PDG: 4.8%)
          R(1)=mpi
          R(2)=mpi
          R(3)=0.
          ty(1)=211
          ty(2)=-211
          ty(3)=22
          FMAJ=5.E-5
        else
          return
        endif
1       CALL GGRSPC(F)
c        print *,cc
c        print 88,(p(1,i),p(2,i),p(3,i),p(4,i),R(i),i=1,3)
c 88     format(5f12.5)
        CALL GRNDM(RND,1)
        IF(IMOD.EQ.3) THEN
c   Mat.element**2. for decay  eta -> pi+pi-pi0
c        /M.Gormley e.a., PR D2(1970)501/
          Q=MPART(IND)-R(1)-R(2)-R(3)
          T0=P(4,3)-mpi0         ! kin.energy of pi0
          Y=3.*(T0/Q-1.)
          FMAT=1.-1.15*Y+0.16*Y**2
        else IF(IMOD.EQ.4) THEN
c   Mat.element**2. for decay  eta -> pi+pi-gamma
c    /A.D.Bukin, S.I.Eidelman, Preprint 77-101/
          P1MOM=SQRT(P(4,1)**2-mpi**2)  ! pi+ momentum
          P2MOM=SQRT(P(4,2)**2-mpi**2)  ! pi- ........
          PGMOM=P(4,3)                       ! gamma ......
          P1P2=P(4,1)*P(4,2)-P1MOM*P2MOM*
     &      (P(1,1)*P(1,2)+P(2,1)*P(2,2)+P(3,1)*P(3,2))
          P1PG=P(4,1)*P(4,3)-P1MOM*PGMOM*
     &      (P(1,1)*P(1,3)+P(2,1)*P(2,3)+P(3,1)*P(3,3))
          P2PG=P(4,2)*P(4,3)-P2MOM*PGMOM*
     &      (P(1,2)*P(1,3)+P(2,2)*P(2,3)+P(3,2)*P(3,3))
          FMAT=2*P1P2*P1PG*P2PG-mpi**2*(P1PG**2+P2PG**2)
        END IF
c
        IF(FMAJ*RND.GT.F*FMAT ) GOTO 1
        IF(F*FMAT.GT.FMAJ) PRINT 920, F*FMAT,FMAJ
920     Format(' GGRDEC3: F,Fmaj=',2(1pE10.2))
c..........
        GAM=PPART(4,ind)/MPART(ind) ! Lorentz-factor of resonance in lab.
        PREST=
     &    DSQRT(PPART(1,ind)**2+PPART(2,ind)**2+PPART(3,ind)**2)
        GG(1)=PPART(1,ind)/PREST ! direction cosines in lab.
        GG(2)=PPART(2,ind)/PREST
        GG(3)=PPART(3,ind)/PREST
c
        DO I=1,3
          npart=npart+1
          MOTHER(npart)=ind
          TYPE(npart)=ty(i)
          MPART(npart)=R(I)
          C(1)=P(1,I)
          C(2)=P(2,I)
          C(3)=P(3,I)
          ENER=P(4,I)
          CALL GGRLOR(GG,GAM,C,ENER,R(I))
          PTOT=DSQRT(ENER**2-R(I)**2)
          PPART(4,npart)=ENER
          PPART(1,npart)=PTOT*C(1)
          PPART(2,npart)=PTOT*C(2)
          PPART(3,npart)=PTOT*C(3)
        enddo
c pi0 decays
        i1=npart-2
        i2=npart
        do i=i1,i2
          if(TYPE(i).eq.111)CALL GGRDEC1(I)
        enddo
      endif
c
      Return
      End
c
      SUBROUTINE GGRLOR(GG,GAM,C,E,R)
c
c   Lorentz transformation into the lab. system
c
c   Input:
c
c   GG(3) - cosines of moving system in lab.
c   GAM   - gamma-factor of moving system
c   C(3)  - cosines of particle in moving system
c   E     - energy of particle in moving system
c   R     - mass of particle
c
c   Output:
c
c   C(3)  - cosines of particle in lab.syst.
c   E     - energy of particle in lab.system
c
      IMPLICIT none
      real*8 GG(3),GAM,C(3),E,R
      real*8 P_LAB(3),E_CM,P_CM,F,E_LAB,P_LABT

      E_CM=E                 ! part. energy in cms
      P_CM=DSQRT(E**2-R**2)  ! part. mom in cms
      F=P_CM*DSQRT(GAM**2-1.)*(C(1)*GG(1)+C(2)*GG(2)+C(3)*GG(3))
      E_LAB=GAM*E_CM+F
      F=(E_LAB+E_CM)*DSQRT(GAM**2-1.)/(GAM+1.)
c
      P_LAB(1)=P_CM*C(1)+F*GG(1)
      P_LAB(2)=P_CM*C(2)+F*GG(2)
      P_LAB(3)=P_CM*C(3)+F*GG(3)
      P_LABT=DSQRT(P_LAB(1)**2+P_LAB(2)**2+P_LAB(3)**2)
c
c  Output variables
c
      C(1)=P_LAB(1)/P_LABT
      C(2)=P_LAB(2)/P_LABT
      C(3)=P_LAB(3)/P_LABT
      E=E_LAB
      RETURN
      END
c
      SUBROUTINE GGRLOR1(GG,GAM,P)
c
c   Lorentz transformation into the lab. system
c
c   Input:
c
c   GG(3) - cosines of moving system in lab.
c   GAM   - gamma-factor of moving system
c   P(4)  - p(1,2,3) - momentum, p(4) - energy of particle in moving system
c   R     - mass of particle
c
c   Output:
c
c   P(4)  -  in lab.system
c
      IMPLICIT none
      real*8 GG(3),GAM,P(4)
      real*8 E_CM,P_CM,C(3),F,E_LAB

      E_CM=P(4)              ! part. energy in cms
      F=DSQRT(GAM**2-1.)*(P(1)*GG(1)+P(2)*GG(2)+P(3)*GG(3))
      E_LAB=GAM*E_CM+F
      F=(E_LAB+E_CM)*DSQRT(GAM**2-1.)/(GAM+1.)
cc  Output variables
      P(1)=P(1)+F*GG(1)
      P(2)=P(2)+F*GG(2)
      P(3)=P(3)+F*GG(3)
      P(4)=E_LAB
      RETURN
      END

      Real*8 function ggpol(t)
      IMPLICIT NONE
c#include "AfkQed/ggrps.inc"
      include 'ggrps.inc'
      integer i,n1,n2
      real*8 t

      n1=3665
      n2=3665

      do i=1,17
        if(sets(n1).lt.t) then
          n2=abs(n2)/2
          if((t.gt.sets(n1+1)).and.(n2.eq.0))n2=1
        else if(sets(n1).gt.t) then
          n2=-abs(n2)/2
          if(n2.eq.0)n2=-1
        else if(sets(n1).eq.t) then
          n2=0
        end if
        n1=n1+n2
      enddo

      ggpol=setpol(n1)+(setpol(n1+1)-setpol(n1))/
     &  (sets(n1+1)-sets(n1))*(t-sets(n1))

      end
