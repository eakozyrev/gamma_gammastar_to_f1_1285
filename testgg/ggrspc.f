      SUBROUTINE GGRSPC(WSPC)
c*********************************************************************
c                                                                    *
c   Calculation of weight of phase volume for <= 20 particles        *
c   and moments of particles                                         *
c                                                                    *
c  Input:                                                            *
c                                                                    *
c      NPRT       - number of particles                              *
c      RM(1-NPRT) - masses of particles  (GeV/c**2)                  *
c      COMS(1-4)  - initial 4-momentum (P,E) of the system in lab.   *
c                                                                    *
c                                                                    *
c  Output: WSPC - weight, in /GGRPRT/:                               *
c          E(1-3,.) - direction cosines of part. in init. system     *
c          E(4,.) - energy of part. in init. system (GeV)            *
c                                                                    *
c                                                                    *
c*********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ND=201)       !  ND=1/SKSI+1
c
      COMMON /GGRARI/  COMS(4),RM(20)
      COMMON /GGRPRT/  E(4,20),NPRT
      COMMON /GGRAAA/  A(ND,18)
      DIMENSION COMPO(3),COMP(3),CP(3),EFM(20)
      REAL*8 WSPC
      SAVE SKSI,PI
c
      DATA IAITR/1/
c
      IF(IAITR.EQ.0) GOTO 25
c
      IAITR=0
      SKSI=1./(ND-1)  ! SKSI - step in ksi
      PI=ACOS(-1.)
c
      CALL GGRTABL         ! preparation of tables for simulation

c--- Preliminary calculations ---------------------------------------

25    COEF=1.               ! beginning of simulation
      ITRUE=0
      TP1=0.
c
      DO I=1,3
        CP(I)=COMS(I)         ! COMS(1-3) momentum of gg-system
        TP1=TP1+CP(I)*CP(I)   !
      END DO
c
      OM=COMS(4)            !  energy of gg-system in lab
      IF(TP1.GT.0.) GOTO 26 !
      EFM(NPRT)=COMS(4)     !  M_eff for gg system in rest
      GOTO 27
c
26    TP1=OM*OM-TP1         ! (E1+ E2)**2-(P1 + P2)**2=M_eff**2
c
      PRINT 80, TP1
      IF(TP1.LT.0.) THEN
        PRINT 80, TP1
        STOP
      END IF
      EFM(NPRT)=SQRT(TP1)   ! eff. mass for moving gg-system
27    RMU=RM(1)             ! 1-st particle mass
c
      DO I=2,NPRT
        RMU=RMU+RM(I)         ! RMU=sum of particle masses
      END DO
c
      TK=EFM(NPRT)-RMU      ! energy deposition in c.m.s.
      TP1=TK**((3.*NPRT-5.)/2.)/(2.*EFM(NPRT))  ! factor
c
      COEF=COEF*TP1
c
      N=NPRT  ! N=number of particles
      N1=N/2
      IF(2*N1.ne.N)then
        N1=3*(N-1)/2         ! N - odd
        T=PI**N1
        N1=N1-1
        T1=1.
        DO I=1,N1
          T1=T1*I
        END DO       
        T=T/T1               ! T - factor for odd N
      else
        T=PI**(3*(N-1)/2.)   ! N - even
        N1=3*N/2-2
        T1=SQRT(PI)/2.
        DO I=2,N1
          T1=T1*(2*I-1)/2.
        END DO
        T=T/T1                ! T - factor for even N
      endif
      TP1=T
c
      COEF=COEF*TP1         ! new COEF
3     I=NPRT-1              ! I = number of particles - 1
c
c--------------- Recurrence part ------------------------------------
c
      DO 12 I1=1,I          ! recurrence part of algorithm
        I2=NPRT+1-I1          ! I2=NPRT,NPRT-1,...,2
        TP1=TK                ! energy deposition
        I3=I2-1               ! I3=NPRT-1,..,1
c
        if(I2.le.2)then
          TK=0.                 
          GOTO 6              ! end of recurrence part. I2=2
        endif
c                           ! I2=3,... , I3=2,...
        K=I3
        K1=K-1                ! K1=NPRT-2,...,1
c
        IF(K1.LT.1.OR.K1.GT.18) PRINT 102, K
        IF(K1.LT.1.OR.K1.GT.18) STOP
c
        CALL GRNDM(R,1)
c
        DO 62 I=2,ND
          Z=A(I,K1)
          if(r.le.z)then
            Z=SKSI*(I-1-(R-Z)/(A(I-1,K1)-Z)) ! K4: Z=ksi_k-1
            goto 63
          endif
62      continue
        
63      PSI=Z
        TK=TK*PSI             ! T_k-1=T_k*ksi_k-1
6       RMU=RMU-RM(I2)        ! K5: rest mass of particles 1,2,...,k-1
        TP1=TP1-TK            !
c
        if(TP1.lt.1.E-11)goto 25
c
        TP1=SQRT(TP1)         ! factor for coefficient in weight
7       COEF=COEF/TP1
        EFM(I3)=RMU+TK        ! K6: eff.mass of compound-particle O_k-1
        TP1=RM(I2)*RM(I2)     ! mass**2 of k-th particle
        OMO=(EFM(I2)*EFM(I2)+TP1-EFM(I3)*EFM(I3))/(2.*EFM(I2)) ! K7:
        PIMPO=OMO**2-TP1      ! energy and momentum**2 of k-th particle in O_k
c
        if(PIMPO.lt.1.E-11)goto 25

        PIMPO=SQRT(PIMPO)     ! momentum of k-th particle in O_k
        COEF=COEF*PIMPO       ! factor for coefficient in weight
c
        CALL GRNDM(RAND,1)    ! K8:
c
        ETA=-1.+2.*RAND       ! cos(theta) simulation
c
        CALL GRNDM(RAND,1)    ! K9:
c
        FI=2.*PI*RAND         ! phi simulation
        TP2=1.-ETA**2         !
c
        IF(TP2.lt.1.E-11)then
          TP2=0.
        else
          TP2=SQRT(TP2)         ! K10:
        endif  
c
        TP2=TP2*PIMPO
        COMPO(3)=ETA*PIMPO    ! Px ???????????
        COMPO(1)=TP2*COS(FI)  ! Py
        COMPO(2)=TP2*SIN(FI)  ! Pz
c
        TP2=0.                ! K11:
c
        DO I4=1,3
          TP2=TP2+COMPO(I4)*CP(I4)
        END DO
c
        E(4,I2)=(OM*OMO+TP2)/EFM(I2)  ! k-th particle energy in init. system
        GR=(E(4,I2)+OMO)/(OM+EFM(I2)) ! factor for Lorentz transformation
c
        DO I4=1,3                     ! K12:
          COMP(I4)=COMPO(I4)+CP(I4)*GR  ! particle 3-momentum in init.system
        END DO
c
        TP2=E(4,I2)**2-TP1            ! mom.**2 of k-th part. in init.system
        if(TP2.lt.1.E-11)then
          TP2=0.
          E(3,I2)=1.
          E(1,I2)=0.
          E(2,I2)=0.
        else  
          TP2=SQRT(TP2)
          DO I4=1,3
            E(I4,I2)=COMP(I4)/TP2      ! direction cos. of k-th part. in init.s.
          END DO
        endif  
c
        OM=OM-E(4,I2)              ! K13: energy of part.system in init. r.f.
c
        DO I4=1,3
          CP(I4)=CP(I4)-COMP(I4) ! 3-momentum of 1,...,k-1 particles in init.r.f.
        END DO
c
        if(i2.le.2)goto 13
c
12    CONTINUE ! -------- end of recurrence procedure -----------
c
c
13    DO I4=1,3         ! K14: last particle
        COMP(I4)=CP(I4)   ! 3-momentum of particle in init. r.f.
      END DO
c
      E(4,1)=OM         ! energy of last particle
      TP3=0.
c
      DO I4=1,3
        TP3=TP3+COMP(I4)*COMP(I4)
      END DO
c
      IF(TP3.lt.1.E-11)then
        TP4=0.
        E(3,1)=1.
        E(1,1)=0.
        E(2,1)=0.
      else
        TP4=SQRT(TP3)
38      DO I4=1,3
          E(I4,1)=COMP(I4)/TP4   ! direction cosines of last particle
        END DO
      endif  
c
41    WSPC=COEF              ! weight of point
c
80    Format(/,' ****** GGRSPS: TP1=',E12.5,/)
102   FORMAT(' GGRSPC: Call with wrong K: K=',I5,' out of 2-19')
900   Format(/,' Version of GGRSPS of 15.03.2000',/)
      Return
      End
c
      SUBROUTINE GGRTABL
      IMPLICIT REAL*8 (A-H,O-Z)
c
c  Calculation of /GGRAAA/ with step in ksi=0.005
c
      PARAMETER (ND=201)       !  ND=1/SKSI+1
c
      COMMON /GGRAAA/ A1(ND,4),A2(ND,4),A3(ND,4),A4(ND,4),A5(ND,2)
      COMMON /GGRMX/  IDUM(9),K
      DIMENSION ALP(ND,4)
      REAL*8 ALP,ANORM,AMAX,GGRGAUS,GGRFUN
      EXTERNAL GGRFUN
c
      SKSI=1./(ND-1)    ! SKSI - step in ksi
c
c  Calculation of alpha_k
c
      DO 10 N=1,5
      K1MAX=4
      IF(N.EQ.5) K1MAX=2
c
      DO 30 K1=1,K1MAX
      K=K1+1+4*(N-1)
c
      ANORM=GGRGAUS(GGRFUN,0.D0,1.D0,1.D-4)  ! normalization
c
c      PRINT 33, K,ANORM
c
c      DO 31 I=1,201   ! cycle on the 1-st index
      DO 31 I=1,ND   ! cycle on the 1-st index
c
      ALP(I,K1)=0.D0
      IF(I.EQ.1) GOTO 31
c      AMAX=5.D-3*(I-1)
      AMAX=SKSI*(I-1)
      ALP(I,K1)=GGRGAUS(GGRFUN,0.D0,AMAX,1.D-4)/ANORM
31    CONTINUE
c
30    CONTINUE
c
c      PRINT 12, 4*(N-1)+2,4*(N-1)+5
c--------------------------------------
      do i=1,ND
c      PRINT 13, ( ALP(i,j),j=1,4)
c
      JMAX=4
      IF(N.EQ.5) JMAX=2
c
      do j=1,JMAX
      IF(N.EQ.1) A1(i,j)=alp(i,j)
      IF(N.EQ.2) A2(i,j)=alp(i,j)
      IF(N.EQ.3) A3(i,j)=alp(i,j)
      IF(N.EQ.4) A4(i,j)=alp(i,j)
      IF(N.EQ.5) A5(i,j)=alp(i,j)
      end do
c
      end do
c
10    CONTINUE
c--------------------------------------
c
12    Format (/,' ALP at k=',I2,'-',I2,/)
13    Format (1x,4(1pE15.5))
20    FORMAT(' z=',1PE12.3,/)
33    FORMAT( ' K=',I2,'  ANORM=',1pE12.3,/)
      Return
      End
c
      REAL*8 FUNCTION GGRGAUS(F,A,B,EPS)
c
c   Integration subroutine
c
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W(12),X(12)
      EXTERNAL F
c
      DATA CONST /1.0D-12/
c
      DATA W
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      DATA X
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/
c
      DELTA=CONST*DABS(A-B)
      GGRGAUS=0.
      AA=A
c
5     Y=B-AA
      IF (DABS(Y).LE.DELTA) RETURN
c
2     BB=AA+Y
      C1=0.5*(AA+BB)
      C2=C1-AA
      S8=0.
      S16=0.
c
      DO 1 I = 1,4
      U=X(I)*C2
1     S8=S8+W(I)*(F(C1+U)+F(C1-U))
c
      DO 3 I = 5,12
      U=X(I)*C2
3     S16=S16+W(I)*(F(C1+U)+F(C1-U))
c
      S8=S8*C2
      S16=S16*C2
      IF (DABS(S16-S8).GT.EPS*(1.0+DABS(S16))) GOTO 4
      GGRGAUS=GGRGAUS+S16
      AA=BB
      GO TO 5
c
4     Y=0.5*Y
      IF (DABS(Y).GT.DELTA) GOTO 2
c
      WRITE (2,7)
      GGRGAUS=0.
c
7     FORMAT('  GGRGAUS: ... TOO HIGH ACCURACY REQUIRED')
      RETURN
      END
c
      REAL*8 FUNCTION GGRFUN(X)
c
      COMMON /GGRMX/ IDUM(9),K
      REAL*8 X,P
c
      P=(3.*K-5.)/2.
      GGRFUN=X**P*DSQRT(1.-X)
c
      RETURN
      END
c
      Subroutine GRNDM(VEC,N)
      implicit none
      integer N,i
      REAL*8 vec(N)
      real*4 rvec
      do i=1,N
        call RANLUX(rvec,1)
        vec(i)=rvec
      enddo
c      call afkrndm(VEC,N)
      RETURN
      END
      

