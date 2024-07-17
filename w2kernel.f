c      This program, w2kernel, converts weight function in frequency to 
c      kernel function in time domain. 
c      w2kernel works interactvely, just start the binary file, then
c      1. you need to give the number of beads in the PIMD simulation
c      2. the temperature of the simulation must be entered
c      3. the program asks for the timestep between two snapshots
c      The weight function is expected in the file of "gx_function.dat". 
c      The weight function should be given from 0 to 11999.95 with 0.05
c      step size. The weight function can be 
c       a) determined from QTB simulation in CP2K 
c       b) determined from a fortran code:
c          https://github.com/madaraszadam/PI_weight_functions
c          https://doi.org/10.26434/chemrxiv-2024-x36vm-v2
c       c) downloaded from this site:
c          https://zenodo.org/records/10702413
c      The resulting kernel function is written in the file of 
c      "kernel.dat"
c      The algorithm is based on Fourier transformation that can be
c      found in the following paper:
c      Dénes Berta, Dávid Ferenc, Imre Bakó and Ádám Madarász 
c      "Nuclear Quantum Effects from the Analysis of Smoothed Trajectories: 
c      Pilot Study for Water"
c      https://doi.org/10.1021/acs.jctc.9b00703
c      email: madarasz.adam@ttk.hu
c     

      program w2kernel

      implicit none

      include "fftw3.f"

      integer*8 plan

c These parameters determines the numerical integrals:
      DOUBLE PRECISION ftstep,kercut,precint

      parameter (ftstep=1.0d-2)    ! step size in the Fourier transform
      parameter (kercut=1.0d-10)    ! cutoff value of the kernel function
      parameter (precint=1.0d-12)  ! precision of the integral

      DOUBLE PRECISION avogadro,boltzmann,planck

      parameter (avogadro=6.02214129d+23)
      parameter (boltzmann=0.831446215d0)
      parameter (planck=6.62606957d-34)
      integer gnum        ! number of table for the weight function
      integer i,j,k,m,ipos,tav,midpos,nst,glim
      integer datn  
      integer ierror,iline
      integer pbead            ! P: number of replacs or beads
      DOUBLE PRECISION kernel(0:2000000),PI,temp
      DOUBLE PRECISION finker(0:2000000)
      DOUBLE PRECISION dnu,sum,koz,timestep,deltat,v,seg,ana,datf
      DOUBLE PRECISION, allocatable :: gx(:,:),gy(:),ftg(:)
      DOUBLE PRECISION, allocatable :: gynew(:),gyold(:)
      DOUBLE PRECISION gmax,dx!,loc,poc,relx,egyik,masik
      DOUBLE PRECISION absx,step
      common pbead
      DOUBLE PRECISION gweight
      DOUBLE PRECISION outdx,outgmax
      LOGICAL :: gx_file_exist
      character*120 gx_file_name,kernel_file_name
      INTEGER outgnum
      DOUBLE PRECISION tot,mindev,dev
      INTEGER devpos

      gx_file_name="gx_function.dat"

      kernel_file_name="kernel.dat"

      PI=4.D0*DATAN(1.D0)

      write (*,821)
  821 format (/,' Number of beads:',
     &             '',$)
c  921    format (i20.0)
         read (*,*) pbead
         write (*,*) pbead

      write (*,81)
   81 format (/,' Temperature for filtration',
     &             ' in Kelvin:  ',$)
   91    format (f20.0)
         read (*,91) temp
         write (*,*) temp

      write (*,80)
   80 format (/,' Timestep between frames',
     &              ' in femtoseconds:  ',$)
c   90 format (f20.0)
      read (*,*) datf

      gnum=240000

      gmax=24000.0d0

      timestep=avogadro*planck/boltzmann/temp*1E11/gmax/2.0d0      

      datn=max(2*floor(datf/2000.0d0/timestep),1)

      timestep=datf/1000.0d0/dble(datn)

      outgmax=avogadro*planck/boltzmann/temp*1E11/timestep/2.0d0

      write (*,*) "datn, timestep",datn,timestep

      outgnum=gnum

      dx=gmax/dble(gnum)

      step=1.0d0/dble(pbead)**2.0d0

      allocate (gx(0:gnum,0:pbead-1))
      allocate (gy(0:gnum))
      allocate (gynew(0:gnum))
      allocate (gyold(0:gnum))
      allocate (ftg(0:gnum))

c Read the 

      INQUIRE(FILE=gx_file_name, EXIST=gx_file_exist)

      if (gx_file_exist) then

        OPEN (51, file = gx_file_name)

         i = 0
         DO
           READ (51,*, END=13) gx(i,0),gy(i)

           gx(i,0)=gx(i,0)*2.0d0

           gyold(i)=gy(i)

           i = i + 1
         END DO
   13 CLOSE (51)

      endif

      outdx = outgmax/outgnum

      do i=0,outgnum

         absx = i*outdx

         seg=gweight(absx,dx,gyold,gnum) 

         gynew(i)=dsqrt(seg)-dsqrt(absx/2.0d0/dble(pbead))

      enddo

      write(*,*) 'Calculation of the kernel function:'

      deltat=2*PI/(avogadro*planck/boltzmann/temp/timestep*1E11)

c FFT

      call dfftw_plan_r2r_1d_ (plan,gnum+1,gynew,ftg,
     &  FFTW_REDFT00, FFTW_ESTIMATE )

      call dfftw_execute_ ( plan )

      seg=0.5d0/dble(gnum)

      sum=0.0d0

      do i = 0, gnum
        ftg(i)=ftg(i)*seg
c        write ( *, * ) i, ftg(i)

        sum=sum+ftg(i)

      end do

      call dfftw_destroy_plan_ ( plan )

c end of FFT

       write(*,*) "sum ftg"

       write(*,*) sum

       write(*,*) "FFT finished"

       finker=0.0d0

       j=0

       tot=0.0d0

       do j = 0, gnum

         ana=0.5/dsqrt(PI*(dble(j)+0.5)*deltat)
      ana=ana-sign(0.5,float(j)-0.5)/dsqrt(PI*(abs(dble(j)-0.5))*deltat)

         seg=ana/dsqrt(dble(pbead))

         sum=ftg(j)+seg

         finker(j)=sum

         if (j.eq.0) then

           tot=sum

           mindev=abs(tot-1.0d0)

         else

           tot=tot+2.0d0*sum

         endif

         dev=abs(tot-1.0d0)

         if (dev.lt.mindev) then

           devpos=j

           mindev=dev

         endif

      enddo

c determine the kernel function on less dense grid

      do j=devpos+1,gnum

        finker(j)=0.0d0

      enddo

     
      j=0
      sum=1.0d0

      seg=0.0d0

      do while (abs(sum) .gt. 0.0d0)

        sum=finker(abs(j*datn-datn/2))+finker(j*datn+datn/2)

        sum=sum/2

        do i=1-datn/2,datn/2-1

          sum=sum+finker(abs(j*datn+i))

        enddo

        kernel(j)=sum

        seg=seg+sum

        j=j+1

      enddo

      seg=seg*2.0d0-kernel(0)

      open (unit=53,file=kernel_file_name)

      do i=0,j-1

        kernel(i)=kernel(i)/seg

        write(53,*) kernel(i)

      enddo

      close(53)

c end determine the kernel function on less dense grid

      write(*,*) "qfilter terminated normally."

      end program w2kernel

      function gweight(x,dx,gy,gnum)
      implicit none
      integer maxsite
      parameter (maxsite=10000)
      integer i,j,k,gnum
      integer ixa,ixb
      integer pbead
      DOUBLE PRECISION x,gweight,dx,gy(0:gnum)
      DOUBLE PRECISION xa,xb,ya,yb,gmax,absx
      DOUBLE PRECISION loc,poc,relx,seg,sumu,sumb
      DOUBLE PRECISION kernel(0:100000),PI,temp
      DOUBLE PRECISION xk(0:1000),xksq(0:1000)
      DOUBLE PRECISION aex,bex,egyik,masik


      common pbead

      PI=4.D0*DATAN(1.D0)

      absx=abs(x)

      gmax = dx*gnum

      if (absx .eq. gmax) then

          gweight=gy(gnum)

          return

      endif

      if (absx .gt. gmax) then

           gweight = absx / 2 / pbead

      else

          ixa = int(absx/dx)
          ixb = ixa + 1
          xa= ixa * dx
          xb= ixb * dx

          ya= gy(ixa)
          yb= gy(ixb)

c linear interpolation:

          gweight = ((absx-xa)*yb+(xb-absx)*ya)/(xb-xa)

      endif

      return
      end
