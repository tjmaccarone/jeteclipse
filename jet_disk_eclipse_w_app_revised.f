      program eclipsejet
c     Copyright 2020 Thomas J. Maccarone
c     All code written by first author of paper
c     Conceptual input from A. Jakob van den Eijnden, Thomas D. Russell
c     and Nathalie Degenaar in developing the science
c     Paper by these authors under submission to MNRAS at the time
c     of posting of this code.
c     This code may only be modified with permission of the lead author
c     but may be download and used as is by anyone.
      implicit none
      real nu_min,nu_max
      integer nfreq,hbins,i,check_eclipse,eclipsed
      real h_min,n,alpha,f_flat,h_max,h,nu_break
      real jetspec(100000),jetspectot(100000),appjetspectot(100000)
      real appjetspec(100000)
      real phase, inclination, r,a,pi,m1,m2,beta_jet
      real delta_app,delta_counter,reldelta
      character*80 outfile
   
      pi=4.0*atan(1.0)
      call getinputs (alpha,nfreq,nu_min,nu_max,nu_break,h_min,h_max,
     &     hbins,f_flat,m1,m2,a,r,phase,inclination,outfile,beta_jet)
      delta_app=reldelta(beta_jet,inclination)
      delta_counter=reldelta(-beta_jet,inclination)
      write (*,*)'deltas=',delta_app,delta_counter


      open (50,status='unknown',file=outfile)
c      write(*,*)'phase, in units of fraction of orbital period=',phase


c     initialize the spectrum
      do i=1,nfreq
         jetspectot(i)=0.0
         appjetspectot(i)=0.0
      enddo

      call getjetspec(nfreq,h,h_min,h_max,alpha,f_flat,
     &     nu_min,nu_max,nu_break,phase, inclination,r,a,m1,m2,
     &     jetspectot,hbins)
c     Above does counter jet, below does the approaching jet
      call getappjetspec(nfreq,h,h_min,h_max,alpha,f_flat,
     &     nu_min,nu_max,nu_break,phase, inclination,r,a,m1,m2,
     &     appjetspectot,hbins)
c     Write out the summed jet spectrum
      call sumjetspecs(jetspectot,appjetspectot,nfreq,delta_app,
     &     delta_counter)
      call writespec(jetspectot,nfreq,nu_min,nu_max,50)
      end


      subroutine sumjetspecs(jetspectot,appjetspectot,nfreq,delta_app,
     &     delta_counter)
      real jetspectot(100000),appjetspectot(100000),relexp
      integer nfreq,i
      real delta_app,delta_counter
      relexp=2.22
c     This should be (7+3p)/(4+p) where p is the electron spectral index
      do i=1,nfreq
         jetspectot(i)=delta_counter**relexp *jetspectot(i)+
     &        delta_app**relexp*appjetspectot(i)
      enddo
      return
      end




      subroutine getjetspec(nfreq,h,h_min,h_max,alpha,f_flat,
     &     nu_min,nu_max,nu_break,phase, inclination,r,a,m1,m2,
     &     jetspectot,hbins)
      implicit none
      real nu_max,nu_min
      integer nfreq,hbins,i,check_eclipse,eclipsed,check_eclipse_disk
      real h_min,n,alpha,f_flat,h_max,h,nu_break
      real jetspec(100000),jetspectot(100000)
      real phase, inclination, r,a,pi,m1,m2,r_out,rcirc
      r_out=rcirc(a,m1,m2)
      do i=1,hbins
         h=h_min*(h_max/h_min)**(float(i)/float(hbins))
         call makespec(jetspec,nfreq,h,h_min,alpha,f_flat,nu_min,
     &     nu_max,nu_break)
         eclipsed=check_eclipse(phase,inclination,r,a,h,m1,m2)+
     &        check_eclipse_disk(h,r_out,inclination)
         if (eclipsed.ge.1) write(*,*)'eclipsed at',h
         if (eclipsed.eq.0) call addjetspec(jetspec,jetspectot,nfreq)
c         if (i.eq.100) call writespec(jetspec,nfreq,nu_min,nu_max,60)         
c         if (i.eq.100) then
c            write(*,*)'writeout of spectrum triggered'
c            write(*,*)'did I get here?'
c            call writespec(jetspectot,nfreq,nu_min,nu_max,60)
      enddo
      return
      end

      subroutine getappjetspec(nfreq,h,h_min,h_max,alpha,f_flat,
     &     nu_min,nu_max,nu_break,phase, inclination,r,a,m1,m2,
     &     jetspectot,hbins)
      implicit none
      real nu_max,nu_min
      integer nfreq,hbins,i,check_eclipse,eclipsed,check_eclipse_disk
      real h_min,n,alpha,f_flat,h_max,h,nu_break
      real jetspec(100000),jetspectot(100000)
      real phase, inclination, r,a,pi,m1,m2,r_out,rcirc
      do i=1,hbins
         h=h_min*(h_max/h_min)**(float(i)/float(hbins))
         call makespec(jetspec,nfreq,h,h_min,alpha,f_flat,nu_min,
     &     nu_max,nu_break)
         call addjetspec(jetspec,jetspectot,nfreq)

c         if (i.eq.100) call writespec(jetspec,nfreq,nu_min,nu_max,60)         
c         if (i.eq.100) then
c            write(*,*)'writeout of spectrum triggered'
c            write(*,*)'did I get here?'
c            call writespec(jetspectot,nfreq,nu_min,nu_max,60)
      enddo
      return
      end



      subroutine writespec(jetspectot,nfreq,nu_min,nu_max,nwrite)
c     This writes out the jet spectrum
      real jetspectot(100000),nu_min,nu_max
      integer nfreq,write
      do i=1,nfreq
         write(nwrite,*) i,frequency(i,nu_min,nu_max,nfreq)
     &        ,jetspectot(i)
      enddo
      return
      end



      subroutine addjetspec(jetspec,jetspectot,nfreq)
c     This adds the spectrum at the current height to the spectra at the
c     other heights
      real jetspec(100000),jetspectot(100000)
      integer i,nfreq
      do i=1,nfreq
         jetspectot(i)=jetspectot(i)+jetspec(i)
      enddo
      return
      end
      


      real function frequency(i,nu_min,nu_max,nfreq)
      integer i,nfreq
      real nu_min,nu_max
      frequency=nu_min*(nu_max/nu_min)**(float(i)/float(nfreq))
      return
      end

      subroutine makespec(jetspec,nfreq,h,h_min,alpha,f_flat,nu_min,
     &     nu_max,nu_break)
c     This subroutine computes the spectrum at a particular height
c     within the jet
      real jetspec(100000)
      real h,h_min,alpha,f_flat,nu_min,nu_max,nu,nu_break,jetflux
      integer i,nfreq
      do i=1,nfreq
         nu=frequency(i,nu_min,nu_max,nfreq)
         jetspec(i)=jetflux(nu,nu_break,h_min,h,alpha,f_flat)
      enddo
      return
      end
      

      real function jetflux(nu,nu_break,h_min,h,alpha,f_flat)
      real nu,nu_break,h_min,h,alpha,nu_break_min, nu_peak,f_flat
      nu_peak=nu_break*(h_min/h)
c     Above assumes that the spectrum peaks at a frequency inversely
c     prop to radius up the jet, and gives the flux density for a single
c     frequency at a single location
      if (nu.le.nu_peak) then 
         jetflux=f_flat*(nu/nu_peak)**2.5 
      else
         jetflux=f_flat*(nu/nu_peak)**(-alpha)
      endif
      return
      end

      integer function check_eclipse(phase,inclination,r,a,h,m1,m2)
      implicit none
      real phase, pi,inclination,r,a,h
      real z_min, z_max,z,a1,a2,m1,m2
      real x1,y1,x2,y2,xprime,yprime,zprime,z_jet,deltaz
      check_eclipse=0
      pi=4.0*datan(1.0d0)

      if (phase.ge.(0.5*pi).and.phase.le.(1.5*pi)) goto 37

c     we need to: (1) for the given phase calculate the x,y,z coordinates 
c     for each object (2) translate to put the accretor back at the origin
c     (3) rotate the plane by the inclination angle
c     (4) determine if the jet is eclipsed, taking into account 
c     the inclination effects of where along the jet is being probed
c      write(*,*)'phase=',phase
      a1=a*m2/(m1+m2)
      a2=a*m1/(m1+m2)
      x1=a1*sin(phase+pi)
      y1=a1*cos(phase+pi)
      x2=a2*sin(phase)
      y2=a2*cos(phase)
c      write(*,*)'x1,x2,y1,y2=',x1,x2,y1,y2
c     now translate and rotate the system
      xprime=(x2-x1)
      yprime=(y2-y1)*cos(inclination)
      zprime=a*cos(inclination)
c     now, the h that is passed is the height along the jet in its own coordinate frame.
c     the z coordinate in the jet is h cos i, hypotenuse is the jet-star distance 
      z_jet = h*sin(inclination)
c     now, find the endpoints of the chords for the star at x=0
      if (r.ge.xprime) then
         deltaz = sqrt(r**2.0-xprime**2.0)
c***  Need to think about the line above -- does inclination affect it at all?
         z_min=zprime-deltaz
         z_max=zprime+deltaz
         if ((z_jet.le.z_max).and.(z_jet.ge.z_min)) check_eclipse=1
      endif
c      if (check_eclipse.eq.1) 
c      write(*,*) xprime,yprime,zprime,z_jet,h,z_min,z_max
 37   continue
      return
      end


      integer function check_eclipse_disk(h,r_out,inclination)
      real h,r_out,inclination
      check_eclipse_disk=0
      if (h.le.r_out/tan(inclination)) check_eclipse_disk=1
      return
      end

      
      subroutine getinputs (alpha,nfreq,nu_min,nu_max,nu_break,h_min
     &     ,h_max,
     &     hbins,f_flat,m1,m2,a,r,phase,inclination,outfile,beta_jet)
      implicit none
      real alpha,nu_min,nu_max,nu_break,h_min,h_max,f_flat,m1,m2,a,r
      real phase, inclination,pi,beta_jet
      integer nfreq,hbins
      character*80 outfile
      pi=4.0*atan(1.0)
c     alpha gives the spectral index of the synchrotron spectrum at each height
c      alpha=0.7 
      read(*,*) alpha 
c     nfreq gives the number of frequencies to put into the jet spectrum
c      nfreq=1000
      read(*,*) nfreq
c     nu_min gives the lowest frequency to consider 
c      nu_min=1e8
      read(*,*) nu_min
c     nu_max gives the highest frequency to consider
c      nu_max=1e17
      read(*,*)nu_max
c     nu_break gives the break frequency of the jet
      read(*,*)nu_break
c      nu_break=1e14
c     h_min gives the base height
      read(*,*) h_min
c      h_min=1e8
c     h_max gives the highest height considered
      read(*,*) h_max
c      h_max=1e18
c     sets the number of height bins
      read(*,*)hbins
c      hbins=1000
c     f_flat gives the flux of the flat part of the spectrum
      read(*,*) f_flat
c      f_flat=1.0
c      m1,m2,phase,inclination,r,a
c      choosing values for 4U 1820-30 -- a is in AU times the AU
      read(*,*) m1,m2
c      m1=0.07
c      m2=1.4
      read(*,*)a,r
c      a=0.0136*1.5e13
c      r=0.167*a
c     the mass may be a bit higher for the neutron star, which would also
c     change the mass ratio
      read(*,*) phase
      phase=phase*2*pi
c      phase=0.0
c     phase is defined incorrectly here
      read(*,*) inclination
      inclination=inclination*pi/180.0
c      inclination=30.0*pi/180.0
      read(*,*) beta_jet
      read(*,'(a)')outfile
      return
      end

      real function rcirc(a,m1,m2)
      real a,m1,m2,q
      q=m1/m2     
      rcirc=a*(1.0+q)*(0.5-0.227*log10(q))**4.0
      return
      end


      real function reldelta(beta,theta)
      real gamma,beta,delta
      gamma=1.0/sqrt(1.0-beta**(2.0))
      write(*,*)'gamma=',gamma
      reldelta=1.0/(gamma*(1.0-beta*cos(theta)))
      return
      end
