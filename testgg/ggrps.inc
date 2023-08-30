      REAL*8 sum,es,fm
      integer NOBR,Nact,ngt
      COMMON /GGRSTA/sum,es,fm,NOBR,Nact,ngt
      REAL*8 EB,RMAS,RWID,RM,RG,FMAX,rmax,rmin,Q2max
      integer IR,IMODE,KVDM,RAD
      COMMON /GGRPAR/EB,RMAS,RWID,RG,RM,FMAX,IR,IMODE,KVDM,RAD,
     & rmax,rmin,Q2max
      REAL*8 alpha,PI,EM,mpi0,mpi,meta,metap,mks,mkc,
     &  mrho,mjpsi,mups,ma0,wa0,mf0,wf0,wrho,g12
      COMMON /GGRCON/alpha,PI,EM,mpi0,mpi,meta,metap,mks,mkc,
     &  mrho,mjpsi,mups,ma0,wa0,mf0,wf0,wrho,g12
      real*8 PPART,MPART
      integer npart,TYPE,MOTHER
      COMMON /GGREV/PPART(4,25),MPART(25),TYPE(25),MOTHER(25),npart
      REAL*8 sets,setpol
      COMMON /GGRPOL/ sets(7330),setpol(7330)
