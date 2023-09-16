      program cyl_builder
c
c     This program fills a cylinder with counter-ions and solvent
c     so as to neutralize a specified wall charge and meet a density spec.
c     
c     The fluid atoms are laid down in lines parallel to the axis
c     and packed hexagonally from the center out. The ions are put
c     in the innermost line.
c     
      implicit none
      real*8 rc,xc,pi,rx,ry,rz
      real*8 rsphere,qw
      real*8 sigma,areac
      real*8 rhof_targ,rhof,rf,vf,rh,rl
      real*8 hex(2,6),half,one,rt32,zero,qf,qfion
      integer imol,itag
      integer nf,nh,nfmax,nlayer,il,ih,ihex,istep,nfion
      integer itypef,itypefion
c
c     Calculate cylinder dimensions
c      
      pi = 4.0*datan(1.0d0)

      imol = 0.0
      itypef = 1                       !! Laamps type of solvent
      itypefion = 2                    !! Laamps type of ion
      qf = 0.0                         !! Charge of solvent
      qfion = 1.0                      !! Charge of Ion (assumed to be Na+) 
c
c     Specify surface charge in C/m^2, then convert to 
c     Angstrom units
c
      qw = -0.1                        !! Total charge on wall in C/m^2          
      qw = qw/(1.6022e-19*1.0e20)      !! Charge in Angstroms  
c
c     Now we will work in Angstrom units
c
      rhof_targ = 0.033456             !! target fluid density
      sigma = 3.2                      !! inner radius for ions only
c      rc = 40.0
c      xc = 80.0
      rc = 20.0                        !! Radius of cylinder
      xc = 40.0                        !! length of cylinder
      areac = 2.0*pi*rc*xc             !! volume of cylinder
      nfion = -nint(areac*qw/qfion)    !! number of solvent ions

      if (nfion.le.0) then
         write(6,*) 'qw and qfion have same sign'
         stop
      endif

c     Now the real raidus is equalt to the original radius - (sigma/2), since only
c     Ions will be placed within sigma of the walls
      areac = -qfion*nfion/qw          !! recalculate the volume of the cylinder
      xc = areac/(2.0*pi*rc)           !! recalculate the length of the cylinder
      rhof = rhof_targ                 !! real fluid density
      rf = rc-0.5*sigma                !! real cylinder radius
      vf = pi*rf**2*xc                 !! real volume of fluid
      nf = nint(vf*rhof)               !! real number of fluid atoms (mass)
      rhof = nf/vf                     !! recalculate the real fluid density
c
c     Echo key dimensions to screen
c
      write(6,*) 'Cylinder radius is ',rc
      write(6,*) 'Cylinder length is ',xc
      write(6,*) 'Wall charge density in C/m^2 is ',
     $     qw*(1.6022e-19*1.0e20)
      write(6,*) 'Wall charge density in Angstrom units is ',qw
      write(6,*) 'Number of fluid ions is ',nfion
c
c     Find rank of smallest hexagonal column 
c     that has at least nf sites.
c
      nh = 0                           !! Rank of hexagon
      nfmax = 0                        !! fluid atom counter
      do while (nfmax.lt.nf)
         nh = nh+1                    
         rh = rf/nh                    !! Spacing of hexagon site
         nlayer = nint(xc/rh)          !! Number of hexagon layers
         rl = xc/nlayer                !! Spacing of hexagon layers
         nfmax = (1+6*nh*(nh+1)/2)*nlayer
         write(6,*) nh,nlayer,nfmax,nf
      end do

      write(6,*) 'Target fluid density is ',rhof_targ
      write(6,*) 'Actual fluid density is ',rhof
      write(6,*) 'Number of fluid atoms is ',nf
      write(6,*) 'Maximum number of fluid atoms is ',nfmax
      write(6,*) 'Rank of hexagon is ',nh
      write(6,*) 'Number of hexagon layers is ',nlayer
      write(6,*) 'Spacing of hexagon sites is ',rh
      write(6,*) 'Spacing of hexagon layers is ',rl
c
c     Now generate atom positions
c
      open (unit=10,file = 'cyl_builder_fluid.r3d',status = 'unknown')
      open (unit=11,file = 'cyl_builder_fluid.lmp',status = 'unknown')

      one  = rh
      zero = rh*0.0
      half = rh*0.5
      rt32 = rh*sqrt(3.0)/2.0

      hex(1,1) = half
      hex(2,1) = rt32
      hex(1,2) = one
      hex(2,2) = zero
      hex(1,3) = half
      hex(2,3) =-rt32
      hex(1,4) =-half
      hex(2,4) =-rt32
      hex(1,5) =-one
      hex(2,5) = zero
      hex(1,6) =-half
      hex(2,6) = rt32

      rx = 0.0
      ry = 0.0
      rz = 0.0
      itag = 0
      rsphere = sigma/2.0

      do il = 1,nlayer
         if (itag.lt.nf) then
            itag = itag+1
            if (itag.gt.nfion) then
               write(11,3) itag,imol,itypef,qf,rx,ry,rz
 3             format(3i8,f10.4,3f15.6)
               write(10,4) rx,ry,rz,rsphere
 4             format('2'/4e15.8,' 0.0 0.0 1.0')
            else
               write(11,3) itag,imol,itypefion,qfion,rx,ry,rz
               write(10,5) rx,ry,rz,rsphere
 5             format('2'/4e15.8,' 0.0 1.0 1.0')
            endif
            rx = rx+rl
         end if
      end do

      do ih = 1,nh
         ry = ry+hex(1,5)
         rz = rz+hex(2,5)
         do ihex = 1,6
            do istep = 1,ih
               ry = ry+hex(1,ihex)
               rz = rz+hex(2,ihex)

               rx = 0.0
               do il = 1,nlayer
                  if (itag.lt.nf) then
                     itag = itag+1
                     if (itag.gt.nfion) then
                        write(11,3) itag,imol,itypef,qf,rx,ry,rz
                        write(10,4) rx,ry,rz,rsphere
                     else
                        write(11,3) itag,imol,itypefion,qfion,rx,ry,rz
                        write(10,5) rx,ry,rz,rsphere
                     endif
                     rx = rx+rl
                  end if
               end do

            end do
         end do
      end do
               
      stop
      end




