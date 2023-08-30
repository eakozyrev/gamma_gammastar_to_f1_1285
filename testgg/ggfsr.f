      subroutine ggfsr(t2)
c--------------------------------------------------------
c Generate photon emitted by the final positron.
c The photon is added to list of particles.
c Positron and resonance momenta are modified.
c input parameters:
c t2 - squared momentum transfer for positron
c kmin - minimal photon energy (GeV)
c funm - majorant value for function discribing the photon
c        energy and angular distribution      
c ---------------------------------------------------------------------
      IMPLICIT NONE
c#include "AfkQed/ggrps.inc" 
      include 'ggrps.inc' 
      real*8 t2,kmin,funm,fmax1
      real*8 l,delta,prob
      real*8 rnd,xmin,xmax,ymin,ymax,x,y,ekf,fac1
      real*8 phik,beta,cost,sint,fac2
      real*8 v(4),tmp,cosp,sinp,kf(4)
      real*8 c(3),c1(3),e,e1,gam,mepi
      real*8 ac(3),dzerox,fun
      integer ng
      common /ggfuc/ac
      real*8 ggfu
      external ggfu
      parameter (kmin=0.001)
      parameter (funm=300.)

      if(PPART(4,2).le.max(0.01d0,kmin))return
c probability of FSR
      l=log(-t2/em**2)-1.
      delta=1.-alpha/pi*((log(1./rmax)-0.75)*l+0.25)
      xmin=kmin/PPART(4,2)
      prob=alpha/pi*((log(1./xmin)-0.75)*l+0.25)/delta
      Call GRNDM(rnd,1)
      if(rnd.gt.prob)return
c generate FSR
      fmax1=0.
      ng=0
  898 ng=ng+1
      if(ng.gt.1000)then
        print *,'FSR generation: too many iterations ( > 1000)'
        return
      endif
c generate photon energy according 1/x
      xmax=1.-0.001/PPART(4,2)
      ymin=log(xmin)
      ymax=log(xmax)
      Call GRNDM(rnd,1)
      y=(ymax-ymin)*rnd+ymin
      x=exp(y)
      ekf=PPART(4,2)*x
      fac1=(ymax-ymin)*x
c generate photon angles relative positron,
c beta*dcos(theta)/(1-beta*cos(theta))
      Call GRNDM(rnd,1)
      phik=2.*PI*RND
      Call GRNDM(rnd,1)
      beta=sqrt(1.-em**2/(PPART(4,2)-ekf)**2)
      ymin=log(1.-beta)
      ymax=log(1.+beta)
      y=(ymax-ymin)*rnd+ymin
      cost=(1-exp(y))/beta
      sint=sqrt(1.-cost**2)
      fac2=(ymax-ymin)*(1.-beta*cost)/beta
c photon momentum 
      v(1)=ekf*cos(phik)*sint
      v(2)=ekf*sin(phik)*sint
      v(3)=ekf*cost
      v(4)=ekf
      tmp=sqrt(PPART(1,2)**2+PPART(2,2)**2)
      cosp=PPART(1,2)/tmp
      sinp=PPART(2,2)/tmp
      tmp=sqrt(PPART(1,2)**2+PPART(2,2)**2+PPART(3,2)**2)
      cost=PPART(3,2)/tmp
      sint=sqrt(1.-cost**2)
      kf(1)= v(1)*cost*cosp-v(2)*sinp+v(3)*sint*cosp
      kf(2)= v(1)*cost*sinp+v(2)*cosp+v(3)*sint*sinp
      kf(3)=-v(1)*sint               +v(3)*cost
      kf(4)= v(4)
c energy and momentum positron+pi0 system
      v(1)=ppart(1,2)+ppart(1,3)
      v(2)=ppart(2,2)+ppart(2,3)
      v(3)=ppart(3,2)+ppart(3,3)
      v(4)=ppart(4,2)+ppart(4,3)
      tmp=sqrt(v(1)**2+v(2)**2+v(3)**2)
      v(1)=-v(1)/tmp
      v(2)=-v(2)/tmp
      v(3)=-v(3)/tmp
      mepi=sqrt(v(4)**2-tmp**2)
      gam=v(4)/mepi
c photon in positron+pi0 rest frame
      tmp=DSQRT(kf(1)**2+kf(2)**2+kf(3)**2)
      C1(1)=kf(1)/tmp
      C1(2)=kf(2)/tmp
      C1(3)=kf(3)/tmp
      e1=kf(4)
      CALL GGRLOR(v,gam,C1,e1,0.d0)
c pi0 direction in positron+pi0 rest frame
      tmp=DSQRT(PPART(1,3)**2+PPART(2,3)**2+PPART(3,3)**2)
      C(1)=PPART(1,3)/tmp
      C(2)=PPART(2,3)/tmp
      C(3)=PPART(3,3)/tmp
      e=PPART(4,3)
      CALL GGRLOR(v,gam,C,e,rm)
c new pi0 momentum
      ac(1)=e1
      ac(2)=(c1(1)*c(1)+c1(2)*c(2)+c1(3)*c(3))*e1
      ac(3)=mepi
      if(ggfu(0.d0).gt.0.)goto 898
      tmp=DZEROX(0.d0,mepi/2.d0,1.d-10,10000,ggfu,1)
      e=sqrt(tmp**2+rm**2)
c new positron momentum in lab
      C1(1)=-e1*c1(1)-tmp*c(1)
      C1(2)=-e1*c1(2)-tmp*c(2)
      C1(3)=-e1*c1(3)-tmp*c(3)
      tmp=sqrt(c1(1)**2+c1(2)**2+c1(3)**2)
      C1(1)=c1(1)/tmp
      C1(2)=c1(2)/tmp
      C1(3)=c1(3)/tmp
      e1=sqrt(tmp**2+em**2)
      v(1)=-v(1)
      v(2)=-v(2)
      v(3)=-v(3)
      CALL GGRLOR(v,gam,C1,e1,em)
c new pi0 momentum in lab
      CALL GGRLOR(v,gam,C,e,rm)
c correct photon energy and angular distribution
      beta=sqrt(1.-em**2/e1)
      cost=(c1(1)*kf(1)+c1(2)*kf(2)+c1(3)*kf(3))/kf(4)
      fun=fac1*fac2*(1.-x+x**2/2.)/x/(1.-beta*cost)
      if(fun.gt.fmax1)fmax1=fun
      if(fun.gt.funm)then
        write(6 ,*) 'FSR: Function value exceeds majorantvalue: ',
     &    fun,'>',funm
      endif
      Call GRNDM(rnd,1)
      if(funm*rnd.gt.fun)goto 898
c new momenta
      tmp=dsqrt(e1**2-em**2)
      PPART(1,2)=C1(1)*tmp
      PPART(2,2)=C1(2)*tmp
      PPART(3,2)=C1(3)*tmp
      PPART(4,2)=e1
      tmp=dsqrt(e**2-rm**2)
      PPART(1,3)=C(1)*tmp
      PPART(2,3)=C(2)*tmp
      PPART(3,3)=C(3)*tmp
      PPART(4,3)=e
      NPART=npart+1
      type(npart)=22
      mother(npart)=0
      MPART(npart)=0.
      PPART(1,npart)=kf(1)
      PPART(2,npart)=kf(2)
      PPART(3,npart)=kf(3)
      PPART(4,npart)=kf(4)

      END
      
      real*8 function ggfu(x)      
      implicit none
c#include "AfkQed/ggrps.inc" 
      include 'ggrps.inc' 
      real*8 ac(3),x
      common /ggfuc/ac
      ggfu=ac(1)+sqrt(x**2+ac(1)**2+2.*ac(2)*x+em**2)+
     &  sqrt(x**2+rm**2)-ac(3)
      end
