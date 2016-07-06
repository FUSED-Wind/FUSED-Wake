c ----------------------------------------------------------------------
c
c GC Larsen wake model applied to offshore wind farms (WindFarm object)
c by Juan P. Murcia <jumu@dtu.dk>
c
c ----------------------------------------------------------------------


c ----------------------------------------------------------------------
c get_R96(D,CT,TI)
c ----------------------------------------------------------------------
c Computes the wake radius at 9.6D downstream location of a turbine
c
c
c Inputs
c ----------
c D (float): Wind turbine diameter
c CT (float): Outputs WindTurbine object's thrust coefficient
c TI (float): Ambient turbulence intensity
c
c
c Outputs
c ----------
c R96 (float): Wake radius at 9.6D downstream location
c
      subroutine get_R96(D,CT,TI,a1,a2,a3,a4,b1,b2,R96)

      implicit none
      real(kind=8) D,CT,TI,a1,a2,a3,a4,b1,b2,R96
cf2py real(kind=8) intent(in) :: D, CT, TI
cf2py real(kind=8) optional,intent(in) :: a1=0.435449861
cf2py real(kind=8) optional,intent(in) :: a2=0.797853685
cf2py real(kind=8) optional,intent(in) :: a3=-0.124807893
cf2py real(kind=8) optional,intent(in) :: a4=0.136821858
cf2py real(kind=8) optional,intent(in) :: b1=15.6298
cf2py real(kind=8) optional,intent(in) :: b2=1.0
cf2py real(kind=8) intent(out) :: R96

      R96=a1*(EXP(a2*CT*CT+a3*CT+a4))*(b1*TI+b2)*D
      end subroutine get_R96


c ----------------------------------------------------------------------
c get_RW(x,D,CT,TI)
c ----------------------------------------------------------------------
c Computes the wake radius at a location
c
c Inputs
c ----------
c x (float): Distance between turbines in the stream-wise direction
c D (float): Wind turbine diameter
c TI (float): Ambient turbulence intensity
c CT (float): Outputs WindTurbine object's thrust coefficient
c
c Outputs
c ----------
c Rw (float): Wake radius at a location
c xT_st (float): Distance from rotor disc to origin of wake expansion
c                (wake becomes zero)
c c1(float): Integration constant
      subroutine get_RW(n,x,D,CT,TI,a1,a2,a3,a4,b1,b2,RW,xT_st,c1)

      implicit none
      integer :: n
      real(kind=8) :: x(n),D,CT,TI,a1,a2,a3,a4,b1,b2
      real(kind=8) :: RW(n),xT_st,c1
cf2py integer intent(hide),depend(x) :: n=len(x)
cf2py real(kind=8) intent(in),dimension(n) :: x
cf2py real(kind=8) intent(in) :: D,CT,TI
cf2py real(kind=8) optional,intent(in) :: a1=0.435449861
cf2py real(kind=8) optional,intent(in) :: a2=0.797853685
cf2py real(kind=8) optional,intent(in) :: a3=-0.124807893
cf2py real(kind=8) optional,intent(in) :: a4=0.136821858
cf2py real(kind=8) optional,intent(in) :: b1=15.6298
cf2py real(kind=8) optional,intent(in) :: b2=1.0
cf2py real(kind=8) intent(out), depend(n),dimension(n) :: RW
cf2py real(kind=8) intent(out) :: xT_st, c1
      ! internal variables
      real(kind=8), parameter :: pi=3.1415926535897932384626433832795d0
      integer :: i
      real(kind=8) :: Area,a,k,R96,tm1,tm2,tm3,tm4,tm5

      Area=pi*D*D/4.0d0
      ! axial velocity deficit at rotor disc
      a=(1.0d0-(sqrt(1.0d0-CT)))/2.0d0
      ! Near wake expansion ratio: k = D0/D from momentum balance
      k=sqrt((1.0d0-a)/(1.0d0-2.0d0*a))
      call get_R96(D,CT,TI,a1,a2,a3,a4,b1,b2,R96)
      ! Distance from rotor disc to origin of wake expansion, where the
      ! wake becomes zero; assuming wake expansion RW ~ (x+xT_st)**(1/3)
      xT_st=(9.6d0*D)/(((2.0d0*R96/(k*D))**3.0d0)-1.0d0)
      do i = 1, n
        if ((x(i)+xT_st) <= 0) then
          ! Null wake radius for locations in front of the rotor disc
          RW(i) = 0.0d0
        else
          ! Integration constant from first order asymptotic solution of
          ! axis-symmetrical momentum balance in the wake
          tm1=(k*D/2.0d0)**(2.5d0)
          tm2=(105.0d0/(2.0d0*pi))**(-0.5d0)
          tm3=(CT*Area*xT_st)**(-5.0d0/6.0d0)
          c1=tm1*tm2*tm3

          tm4=(105.0d0*c1*c1/(2.0d0*pi))**(0.2d0)
          tm5=(CT*Area*(x(i)+xT_st))**(1.0d0/3.0d0)
          RW(i)=tm4*tm5
        end if
      end do
      end subroutine get_RW

c ----------------------------------------------------------------------
c get_dU(x,r,D,CT,TI)
c ----------------------------------------------------------------------
c Computes the wake velocity deficit at a location
c
c Inputs
c ----------
c x (float): Distance between turbines in the stream-wise direction
c r (array): Radial distance between the turbine and the location
c D (float): Wind turbine diameter
c TI (float): Ambient turbulence intensity [-]
c CT (float): Outputs WindTurbine object's thrust coefficient
c
c Outputs
c ----------
c dU (float): Wake velocity deficit at a location normalized
c             by rotor averaged (equivalent) inflow velocity
      subroutine get_dU(n,x,r,D,CT,TI,a1,a2,a3,a4,b1,b2,dU)

      implicit none
      integer :: n
      real(kind=8) :: x(n),r(n),D,CT,TI,a1,a2,a3,a4,b1,b2,dU(n)
cf2py integer intent(hide),depend(x) :: n=len(x)
cf2py real(kind=8) intent(in),dimension(n) :: x
cf2py real(kind=8) intent(in),depend(n),dimension(n) :: r
cf2py real(kind=8) intent(in) :: D,CT,TI
cf2py real(kind=8) optional,intent(in) :: a1=0.435449861
cf2py real(kind=8) optional,intent(in) :: a2=0.797853685
cf2py real(kind=8) optional,intent(in) :: a3=-0.124807893
cf2py real(kind=8) optional,intent(in) :: a4=0.136821858
cf2py real(kind=8) optional,intent(in) :: b1=15.6298
cf2py real(kind=8) optional,intent(in) :: b2=1.0
cf2py real(kind=8) intent(out),depend(n),dimension(n) :: dU
      ! internal variables
      real(kind=8), parameter :: pi=3.1415926535897932384626433832795d0
      integer :: i
      real(kind=8) :: Area
      real(kind=8) :: tm10,tm20,tm30,tm31,tm32,tm40,tm41,tm42,xT_st,c1
      real(kind=8), dimension(n) :: RW

      Area=pi*D*D/4.0d0

      call get_RW(n,x,D,CT,TI,a1,a2,a3,a4,b1,b2,RW,xT_st,c1)

      do i = 1, n
        if ((x(i) <= 0).or.(r(i) >= RW(i))) then
          ! Null deficit for locations in front of the rotor disc
          ! and outside the wake radius
          dU(i) = 0.0d0
        else
          tm10=1.0d0/9.0d0
          tm20=(CT*Area*(x(i)+xT_st)**(-2.0d0))**(1.0d0/3.0d0)
          tm31=(r(i)**(1.5d0))
          tm32=(3.0d0*c1*c1*CT*Area*(x(i)+xT_st))**(-0.5d0)
          tm30=tm31*tm32
          tm41=(35.0d0/(2.0d0*pi))**(3.0d0/10.0d0)
          tm42=(3.0d0*c1*c1)**(-0.2d0)
          tm40=tm41*tm42
          dU(i)=-tm10*tm20*(tm30-tm40)**(2.0d0)
        end if
      end do
      end subroutine get_dU

c ----------------------------------------------------------------------
c get_dUeq(x,y,z,DT,D,CT,TI)
c ----------------------------------------------------------------------
c Computes the rotor averaged (equivalent) wake velocity deficit at
c different turbine locations and diameters
c
c Inputs
c ----------
c x (array): Distance between turbines in the stream-wise direction
c y (array): Distance between turbines in the cross-flow direction
c z (array): Distance between turbines in the vertical direction
c DT (array): Wake operating turbines diameter
c D (float): Wake generating turbine diameter
c TI (float): Ambient turbulence intensity [-]
c CT (float): Outputs WindTurbine object's thrust coefficient
c
c Outputs
c ----------
c dUeq (float): Wake velocity deficit at a location normalized by
c               inflow velocity
      subroutine get_dUeq(n,x,y,z,DT,D,CT,TI,a1,a2,a3,a4,b1,b2,dUeq)

      implicit none
      integer :: n
      real(kind=8) :: x(n),y(n),z(n),DT(n),D,CT,TI
      real(kind=8) :: a1,a2,a3,a4,b1,b2,dUeq(n)
cf2py integer intent(hide),depend(x) :: n = len(x)
cf2py real(kind=8) intent(in),dimension(n) :: x
cf2py real(kind=8) intent(in),depend(n),dimension(n) :: y,z,DT
cf2py real(kind=8) intent(in) :: D,CT,TI
cf2py real(kind=8) optional,intent(in) :: a1=0.435449861
cf2py real(kind=8) optional,intent(in) :: a2=0.797853685
cf2py real(kind=8) optional,intent(in) :: a3=-0.124807893
cf2py real(kind=8) optional,intent(in) :: a4=0.136821858
cf2py real(kind=8) optional,intent(in) :: b1=15.6298
cf2py real(kind=8) optional,intent(in) :: b2=1.0
cf2py real(kind=8) intent(out),depend(n),dimension(n) :: dUeq
      ! internal variables
      real(kind=8), parameter :: pi=3.1415926535897932384626433832795d0
      integer :: i,j
      real(kind=8), dimension(40) :: r_node,th_node,weight
      real(kind=8), dimension(1) :: x_e,y_e,z_e,r_e,dU
      real(kind=8), dimension(n) :: RT

      RT = DT/2.0d0
      ! 2D Gauss quadrature nodes and weights for polar integration
      ! int_(0,R) int_(0,2pi) f(r,theta) r dtheta dr =
      ! sum weight f(r_node, theta_node)
      r_node(1) = 0.26349922998554242692
      th_node(1) = 4.79436403870179805864
      weight(1) = 0.00579798753740115753

      r_node(2) = 0.26349922998554242692
      th_node(2) = 5.13630491629471475079
      weight(2) = 0.01299684397858970851

      r_node(3) = 0.26349922998554242692
      th_node(3) = 5.71955352542765460555
      weight(3) = 0.01905256317618122044

      r_node(4) = 0.26349922998554242692
      th_node(4) = 0.20924454049880022999
      weight(4) = 0.02341643323656225281

      r_node(5) = 0.26349922998554242692
      th_node(5) = 1.10309379714216659885
      weight(5) = 0.02569988335562909190

      r_node(6) = 0.26349922998554242692
      th_node(6) = 2.03849885644762496284
      weight(6) = 0.02569988335562912660

      r_node(7) = 0.26349922998554242692
      th_node(7) = 2.93234811309099407950
      weight(7) = 0.02341643323656214179

      r_node(8) = 0.26349922998554242692
      th_node(8) = 3.70522443534172518653
      weight(8) = 0.01905256317618119616

      r_node(9) = 0.26349922998554242692
      th_node(9) = 4.28847304447466459720
      weight(9) = 0.01299684397858971198

      r_node(10) = 0.26349922998554242692
      th_node(10) = 4.63041392206758217753
      weight(10) = 0.00579798753740114539

      r_node(11) = 0.57446451431535072718
      th_node(11) = 4.79436403870179805864
      weight(11) = 0.01086984853977092380

      r_node(12) = 0.57446451431535072718
      th_node(12) = 5.13630491629471475079
      weight(12) = 0.02436599330905551281

      r_node(13) = 0.57446451431535072718
      th_node(13) = 5.71955352542765460555
      weight(13) = 0.03571902745281423097

      r_node(14) = 0.57446451431535072718
      th_node(14) = 0.20924454049880022999
      weight(14) = 0.04390024659093685194

      r_node(15) = 0.57446451431535072718
      th_node(15) = 1.10309379714216659885
      weight(15) = 0.04818117282305908744

      r_node(16) = 0.57446451431535072718
      th_node(16) = 2.03849885644762496284
      weight(16) = 0.04818117282305915683

      r_node(17) = 0.57446451431535072718
      th_node(17) = 2.93234811309099407950
      weight(17) = 0.04390024659093664378

      r_node(18) = 0.57446451431535072718
      th_node(18) = 3.70522443534172518653
      weight(18) = 0.03571902745281418240

      r_node(19) = 0.57446451431535072718
      th_node(19) = 4.28847304447466459720
      weight(19) = 0.02436599330905552321

      r_node(20) = 0.57446451431535072718
      th_node(20) = 4.63041392206758217753
      weight(20) = 0.01086984853977089951

      r_node(21) = 0.81852948743000586429
      th_node(21) = 4.79436403870179805864
      weight(21) = 0.01086984853977090992

      r_node(22) = 0.81852948743000586429
      th_node(22) = 5.13630491629471475079
      weight(22) = 0.02436599330905548505

      r_node(23) = 0.81852948743000586429
      th_node(23) = 5.71955352542765460555
      weight(23) = 0.03571902745281418934

      r_node(24) = 0.81852948743000586429
      th_node(24) = 0.20924454049880022999
      weight(24) = 0.04390024659093679643

      r_node(25) = 0.81852948743000586429
      th_node(25) = 1.10309379714216659885
      weight(25) = 0.04818117282305903193

      r_node(26) = 0.81852948743000586429
      th_node(26) = 2.03849885644762496284
      weight(26) = 0.04818117282305909438

      r_node(27) = 0.81852948743000586429
      th_node(27) = 2.93234811309099407950
      weight(27) = 0.04390024659093658826

      r_node(28) = 0.81852948743000586429
      th_node(28) = 3.70522443534172518653
      weight(28) = 0.03571902745281413383

      r_node(29) = 0.81852948743000586429
      th_node(29) = 4.28847304447466459720
      weight(29) = 0.02436599330905549199

      r_node(30) = 0.81852948743000586429
      th_node(30) = 4.63041392206758217753
      weight(30) = 0.01086984853977088737

      r_node(31) = 0.96465960618086743494
      th_node(31) = 4.79436403870179805864
      weight(31) = 0.00579798753740116100

      r_node(32) = 0.96465960618086743494
      th_node(32) = 5.13630491629471475079
      weight(32) = 0.01299684397858971545

      r_node(33) = 0.96465960618086743494
      th_node(33) = 5.71955352542765460555
      weight(33) = 0.01905256317618123432

      r_node(34) = 0.96465960618086743494
      th_node(34) = 0.20924454049880022999
      weight(34) = 0.02341643323656226669

      r_node(35) = 0.96465960618086743494
      th_node(35) = 1.10309379714216659885
      weight(35) = 0.02569988335562910925

      r_node(36) = 0.96465960618086743494
      th_node(36) = 2.03849885644762496284
      weight(36) = 0.02569988335562914394

      r_node(37) = 0.96465960618086743494
      th_node(37) = 2.93234811309099407950
      weight(37) = 0.02341643323656215567

      r_node(38) = 0.96465960618086743494
      th_node(38) = 3.70522443534172518653
      weight(38) = 0.01905256317618120657

      r_node(39) = 0.96465960618086743494
      th_node(39) = 4.28847304447466459720
      weight(39) = 0.01299684397858972065

      r_node(40) = 0.96465960618086743494
      th_node(40) = 4.63041392206758217753
      weight(40) = 0.00579798753740114886

      do i = 1, n
        dUeq(i) = 0.0d0
        if (x(i) > 0.0) then
          ! Evaluation of wake and sum of quadrature
          do j = 1, 40
            ! Location of evaluation points in the local rotor coordinates
            x_e = x(i)
            y_e = y(i) + RT(i)*r_node(j)*cos(th_node(j))
            z_e = z(i) + RT(i)*r_node(j)*sin(th_node(j))
            r_e  = (y_e**(2.0d0) + z_e**(2.0d0))**(0.5d0)
            call get_dU(1,x_e,r_e,D,CT,TI,a1,a2,a3,a4,b1,b2,dU)
            dUeq(i)=dUeq(i)+dU(1)*weight(j)
          end do
        end if
      end do

      end subroutine get_dUeq


c ----------------------------------------------------------------------
c gcl_s(x,y,z,DT,P_c,CT_c,WS,TI)
c ----------------------------------------------------------------------
c SINGLE FLOW CASE
c Computes the WindFarm flow and Power using G. C. Larsen model:
c Larsen, G. C., and P. E. Re'thore'. A simple stationary semi-analytical
c wake model. Technical Report Risoe, 2009.
c
c Inputs
c ----------
c x_g (array): Distance between turbines in the global coordinates
c y_g (array): Distance between turbines in the global coordinates
c z_g (array): Distance between turbines in the global coordinates
c DT (array): Turbines diameter
c P_c (array): Power curves
c CT_c (array): Thrust coefficient curves
c WS (float): Undisturbed rotor averaged (equivalent) wind speed at hub
c             height [m/s]
c WD (float): Undisturbed wind direction at hub height [deg.]
c             Meteorological coordinates (N=0,E=90,S=180,W=270)
c TI (float): Ambient turbulence intensity [-]
c
c rho (float): Air density at which the power curve is valid [kg/m^3]
c WS_CI (array): Cut in wind speed [m/s] for each turbine
c WS_CO (array): Cut out wind speed [m/s] for each turbine
c CT_idle (array): Thrust coefficient at rest [-] for each turbine
c
c Outputs
c ----------
c P (array): Power production of the wind turbines (nWT,1) [W]
c T (array): Thrust force of the wind turbines (nWT,1) [N]
c U (array): Rotor averaged (equivalent) Wind speed at hub height
c            (nWT,1) [m/s]
      subroutine gcl_s(n,nP,nCT,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,TI,
     &a1,a2,a3,a4,b1,b2,rho,WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),WS,WD,TI,a1,a2,a3,a4,b1,b2
      real(kind=8) :: rho,WS_CI(n),WS_CO(n),CT_idle(n),P(n),T(n),U(n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in) :: WS,WD,TI
cf2py real(kind=8) optional,intent(in) :: a1=0.435449861
cf2py real(kind=8) optional,intent(in) :: a2=0.797853685
cf2py real(kind=8) optional,intent(in) :: a3=-0.124807893
cf2py real(kind=8) optional,intent(in) :: a4=0.136821858
cf2py real(kind=8) optional,intent(in) :: b1=15.6298
cf2py real(kind=8) optional,intent(in) :: b2=1.0
cf2py real(kind=8) optional,intent(in) :: rho=1.225
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CI=4.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CO=25.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: CT_idle=0.053
cf2py real(kind=8) intent(out),depend(n),dimension(n) :: P,T,U
      ! internal variables
      real(kind=8), parameter :: pi=3.1415926535897932384626433832795d0
      integer :: i,j,k,nDownstream(n),idT(n)
      real(kind=8) :: x(n),y(n),z(n),D,CT,dUeq(n),angle
      real(kind=8) :: x_l(n,n),y_l(n,n)

      ! Rotates the global coordinates to local flow coordinates
      angle = pi*(270.0d0-WD)/180.0d0
      do i=1,n
        do j=1,n
          x_l(i,j) = cos(angle)*x_g(i,j)+sin(angle)*y_g(i,j)
          y_l(i,j) = -sin(angle)*x_g(i,j)+cos(angle)*y_g(i,j)
        end do
        ! counts the number of turbines in front of turbine
        nDownstream(i) = count(x_l(i,:).lt.0)
      end do
      ! Indexes of ordered turbines from most upstream turbine
      call order_id(n,nDownstream,idT)
      ! Initializes the rotor averaged (equivalent) velocity
      U = WS
      ! Computes the rotor averaged (equivalent) velocity deficit
      do j=1,n
        i=idT(j)
        x = x_l(i,:)
        y = y_l(i,:)
        z = z_g(i,:)
        D = DT(i)
        if ((U(i) >= WS_CI(i)).and.(U(i) <= WS_CO(i))) then
          call interp_l(CT_c(i,:,1),CT_c(i,:,2),nCT,U(i),CT)
        else
          CT = CT_idle(i)
        end if
        call get_dUeq(n,x,y,z,DT,D,CT,TI,a1,a2,a3,a4,b1,b2,dUeq)
        U = U + U(i)*dUeq
      end do
      ! Calculates the power and thrust
      do k=1,n
        if ((U(k) >= WS_CI(k)).and.(U(k) <= WS_CO(k))) then
          call interp_l(P_c(k,:,1),P_c(k,:,2),nP,U(k),P(k))
          call interp_l(CT_c(k,:,1),CT_c(k,:,2),nCT,U(k),CT)
        else
          P(k)=0.0d0
          CT = CT_idle(k)
        end if
        T(k) = CT*0.5d0*rho*U(k)*U(k)*pi*DT(k)*DT(k)/4.0d0
      end do

      end subroutine gcl_s

c ----------------------------------------------------------------------
c gcl(x,y,z,DT,P_c,CT_c,WS,WD,TI)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES
c Computes the WindFarm flow and Power using G. C. Larsen model:
c Larsen, G. C., and P. E. Re'thore'. A simple stationary semi-analytical
c wake model. Technical Report Risoe, 2009.
c
c Inputs
c ----------
c x_g (array): Distance between turbines in the global coordinates
c y_g (array): Distance between turbines in the global coordinates
c z_g (array): Distance between turbines in the global coordinates
c DT (array): Turbines diameter
c P_c (array): Power curves
c CT_c (array): Thrust coefficient curves
c WS (array): Undisturbed rotor averaged (equivalent) wind speed at hub
c             height [m/s]
c WD (array): Undisturbed wind direction at hub height [deg.]
c             Meteorological coordinates (N=0,E=90,S=180,W=270)
c TI (array): Ambient turbulence intensity [-]
c
c rho (float): Air density at which the power curve is valid [kg/m^3]
c WS_CI (array): Cut in wind speed [m/s] for each turbine
c WS_CO (array): Cut out wind speed [m/s] for each turbine
c CT_idle (array): Thrust coefficient at rest [-] for each turbine
c
c Outputs
c ----------
c P (array): Power production of the wind turbines (nWT,1) [W]
c T (array): Thrust force of the wind turbines (nWT,1) [N]
c U (array): Rotor averaged (equivalent) Wind speed at hub height
c            (nWT,1) [m/s]
      subroutine gcl(n,nP,nCT,nF,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,TI,a1,a2,
     &a3,a4,b1,b2,rho,WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT,nF
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),rho,WS_CI(n),WS_CO(n),CT_idle(n)
      real(kind=8),dimension(nF) :: WS,WD,TI,a1,a2,a3,a4,b1,b2
      real(kind=8) :: P(nF,n),T(nF,n),U(nF,n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py integer intent(hide),depend(WS) :: nF = len(WS)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in),dimension(nF) :: WS,WD,TI
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a1=0.435449861
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a2=0.797853685
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a3=-0.124807893
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a4=0.136821858
cf2py real(kind=8) optional,intent(in),dimension(nF) :: b1=15.6298
cf2py real(kind=8) optional,intent(in),dimension(nF) :: b2=1.0
cf2py real(kind=8) optional,intent(in) :: rho=1.225
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CI=4.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CO=25.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: CT_idle=0.053
cf2py real(kind=8) intent(out),depend(n),dimension(nF,n) :: P,T,U
      ! internal variables
      integer :: i

      do i=1,nF
        call gcl_s(n,nP,nCT,x_g,y_g,z_g,DT,P_c,CT_c,WS(i),WD(i),TI(i),
     &       a1(i),a2(i),a3(i),a4(i),b1(i),b2(i),rho,WS_CI,WS_CO,
     &       CT_idle,P(i,:),T(i,:),U(i,:))
      end do

      end subroutine gcl

c ----------------------------------------------------------------------
c gcl_av(x,y,z,DT,P_c,CT_c,WS,WD,TI,AV)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES with wt available
c Computes the WindFarm flow and Power using G. C. Larsen model:
c Larsen, G. C., and P. E. Re'thore'. A simple stationary semi-analytical
c wake model. Technical Report Risoe, 2009.
c
c Inputs
c ----------
c x_g (array): Distance between turbines in the global coordinates
c y_g (array): Distance between turbines in the global coordinates
c z_g (array): Distance between turbines in the global coordinates
c DT (array): Turbines diameter
c P_c (array): Power curves
c CT_c (array): Thrust coefficient curves
c WS (array): Undisturbed rotor averaged (equivalent) wind speed at hub
c             height [m/s]
c WD (array): Undisturbed wind direction at hub height [deg.]
c             Meteorological coordinates (N=0,E=90,S=180,W=270)
c TI (array): Ambient turbulence intensity [-]
c AV (array): Wind turbine available per flow [nF,n]
c
c rho (float): Air density at which the power curve is valid [kg/m^3]
c WS_CI (array): Cut in wind speed [m/s] for each turbine
c WS_CO (array): Cut out wind speed [m/s] for each turbine
c CT_idle (array): Thrust coefficient at rest [-] for each turbine
c
c Outputs
c ----------
c P (array): Power production of the wind turbines (nWT,1) [W]
c T (array): Thrust force of the wind turbines (nWT,1) [N]
c U (array): Rotor averaged (equivalent) Wind speed at hub height
c            (nWT,1) [m/s]
      subroutine gcl_av(n,nP,nCT,nF,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,TI,AV,
     &a1,a2,a3,a4,b1,b2,rho,WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT,nF,AV(nf,n)
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),rho,WS_CI(n),WS_CO(n),CT_idle(n)
      real(kind=8),dimension(nF) :: WS,WD,TI,a1,a2,a3,a4,b1,b2
      real(kind=8) :: P(nF,n),T(nF,n),U(nF,n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py integer intent(hide),depend(WS) :: nF = len(WS)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in),dimension(nF) :: WS,WD,TI
cf2py integer intent(in),dimension(nF,n) :: AV
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a1=0.435449861
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a2=0.797853685
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a3=-0.124807893
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a4=0.136821858
cf2py real(kind=8) optional,intent(in),dimension(nF) :: b1=15.6298
cf2py real(kind=8) optional,intent(in),dimension(nF) :: b2=1.0
cf2py real(kind=8) optional,intent(in) :: rho=1.225
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CI=4.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CO=25.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: CT_idle=0.053
cf2py real(kind=8) intent(out),depend(n),dimension(nF,n) :: P,T,U
      ! internal variables
      integer :: i,j
      real(kind=8) :: CT_c_AV(n,nCT,2), P_c_AV(n,nCT,2)

      do i=1,nF
        CT_c_AV = CT_c
        P_c_AV  = P_c
        ! Re-defines the trust curve for non available turbines
        do j=1,n
          if (AV(i,j)==0) then
            CT_c_AV(j,:,2) = CT_idle(j)
            P_c_AV(j,:,2) = 0.0d0
          end if
        end do
        call gcl_s(n,nP,nCT,x_g,y_g,z_g,DT,P_c_AV,CT_c_AV,WS(i),WD(i),
     &       TI(i),a1(i),a2(i),a3(i),a4(i),b1(i),b2(i),rho,WS_CI,
     &       WS_CO,CT_idle,P(i,:),T(i,:),U(i,:))
      end do

      end subroutine gcl_av

c ----------------------------------------------------------------------
c New Multiple turbine types (height, P and CT curves) implementation
c ----------------------------------------------------------------------


c ----------------------------------------------------------------------
c gclm_s(x,y,z,DT,P_c,CT_c,WS,TI)
c ----------------------------------------------------------------------
c SINGLE FLOW CASE
c Computes the WindFarm flow and Power using G. C. Larsen model:
c Larsen, G. C., and P. E. Re'thore'. A simple stationary semi-analytical
c wake model. Technical Report Risoe, 2009.
c
c Inputs
c ----------
c x_g (array): Distance between turbines in the global coordinates
c y_g (array): Distance between turbines in the global coordinates
c z_g (array): Distance between turbines in the global coordinates
c DT (array): Turbines diameter
c P_c (array): Power curves
c CT_c (array): Thrust coefficient curves
c WS (array): Undisturbed rotor averaged (equivalent) wind speed at hub
c             height [m/s] for each wind turbine
c WD (array): Undisturbed wind direction at hub height [deg.] for each
c             turbine. Meteorological coordinates (N=0,E=90,S=180,W=270)
c TI (array): Turbulence intensity [-] for each turbine
c
c rho (float): Air density at which the power curve is valid [kg/m^3]
c WS_CI (array): Cut in wind speed [m/s] for each turbine
c WS_CO (array): Cut out wind speed [m/s] for each turbine
c CT_idle (array): Thrust coefficient at rest [-] for each turbine
c
c Outputs
c ----------
c P (array): Power production of the wind turbines (nWT,1) [W]
c T (array): Thrust force of the wind turbines (nWT,1) [N]
c U (array): Rotor averaged (equivalent) Wind speed at hub height
c            (nWT,1) [m/s]
      subroutine gclm_s(n,nP,nCT,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,TI,
     &a1,a2,a3,a4,b1,b2,rho,WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),WS(n),WD(n),TI(n),a1,a2,a3,a4,b1,b2
      real(kind=8) :: rho,WS_CI(n),WS_CO(n),CT_idle(n),P(n),T(n),U(n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in),dimension(n) :: WS,WD,TI
cf2py real(kind=8) optional,intent(in) :: a1=0.435449861
cf2py real(kind=8) optional,intent(in) :: a2=0.797853685
cf2py real(kind=8) optional,intent(in) :: a3=-0.124807893
cf2py real(kind=8) optional,intent(in) :: a4=0.136821858
cf2py real(kind=8) optional,intent(in) :: b1=15.6298
cf2py real(kind=8) optional,intent(in) :: b2=1.0
cf2py real(kind=8) optional,intent(in) :: rho=1.225
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CI=4.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CO=25.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: CT_idle=0.053
cf2py real(kind=8) intent(out),depend(n),dimension(n) :: P,T,U
      ! internal variables
      real(kind=8), parameter :: pi=3.1415926535897932384626433832795d0
      integer :: i,j,k,nDownstream(n),idT(n)
      real(kind=8) :: x(n),y(n),z(n),D,CT,dUeq(n),angle
      real(kind=8) :: x_l(n,n),y_l(n,n)

      ! Rotates the global coordinates to local flow coordinates
      do i=1,n
        angle = pi*(270.0d0-WD(i))/180.0d0
        do j=1,n
          x_l(i,j) = cos(angle)*x_g(i,j)+sin(angle)*y_g(i,j)
          y_l(i,j) = -sin(angle)*x_g(i,j)+cos(angle)*y_g(i,j)
        end do
        ! counts the number of turbines in front of turbine
        nDownstream(i) = count(x_l(i,:).lt.0)
      end do
      ! Indexes of ordered turbines from most upstream turbine
      call order_id(n,nDownstream,idT)
      ! Initializes the rotor averaged (equivalent) velocity
      U = WS
      ! Computes the rotor averaged (equivalent) velocity deficit
      do j=1,n
        i=idT(j)
        x = x_l(i,:)
        y = y_l(i,:)
        z = z_g(i,:)
        D = DT(i)
        if ((U(i) >= WS_CI(i)).and.(U(i) <= WS_CO(i))) then
          call interp_l(CT_c(i,:,1),CT_c(i,:,2),nCT,U(i),CT)
        else
          CT = CT_idle(i)
        end if
        call get_dUeq(n,x,y,z,DT,D,CT,TI(i),a1,a2,a3,a4,b1,b2,dUeq)
        U = U + U(i)*dUeq
      end do
      ! Calculates the power and thrust
      do k=1,n
        if ((U(k) >= WS_CI(k)).and.(U(k) <= WS_CO(k))) then
          call interp_l(P_c(k,:,1),P_c(k,:,2),nP,U(k),P(k))
          call interp_l(CT_c(k,:,1),CT_c(k,:,2),nCT,U(k),CT)
        else
          P(k)=0.0d0
          CT = CT_idle(k)
        end if
        T(k) = CT*0.5d0*rho*U(k)*U(k)*pi*DT(k)*DT(k)/4.0d0
      end do

      end subroutine gclm_s

c ----------------------------------------------------------------------
c gclm(x,y,z,DT,P_c,CT_c,WS,WD,TI)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES
c Computes the WindFarm flow and Power using G. C. Larsen model:
c Larsen, G. C., and P. E. Re'thore'. A simple stationary semi-analytical
c wake model. Technical Report Risoe, 2009.
c
c Inputs
c ----------
c x_g (array): Distance between turbines in the global coordinates
c y_g (array): Distance between turbines in the global coordinates
c z_g (array): Distance between turbines in the global coordinates
c DT (array): Turbines diameter
c P_c (array): Power curves
c CT_c (array): Thrust coefficient curves
c WS (array): Undisturbed rotor averaged (equivalent) wind speed at hub
c             height [m/s] for each wind turbine x flow [nF,n]
c WD (array): Undisturbed wind direction at hub height [deg.] for each
c             turbine x flow [nF,n]
c             Meteorological coordinates (N=0,E=90,S=180,W=270)
c TI (array): Turbulence intensity [-] for each turbine x flow [nF,n]
c
c rho (float): Air density at which the power curve is valid [kg/m^3]
c WS_CI (array): Cut in wind speed [m/s] for each turbine
c WS_CO (array): Cut out wind speed [m/s] for each turbine
c CT_idle (array): Thrust coefficient at rest [-] for each turbine
c
c Outputs
c ----------
c P (array): Power production of the wind turbines (nWT,1) [W]
c T (array): Thrust force of the wind turbines (nWT,1) [N]
c U (array): Rotor averaged (equivalent) Wind speed at hub height
c            (nWT,1) [m/s]
      subroutine gclm(n,nP,nCT,nF,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,TI,a1,
     &a2,a3,a4,b1,b2,rho,WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT,nF
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),rho,WS_CI(n),WS_CO(n),CT_idle(n)
      real(kind=8),dimension(nF) :: a1,a2,a3,a4,b1,b2
      real(kind=8) :: WS(nF,n),WD(nF,n),TI(nF,n),P(nF,n),T(nF,n),U(nF,n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py integer intent(hide),depend(WS) :: nF = size(WS,1)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in),dimension(nF,n) :: WS,WD,TI
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a1=0.435449861
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a2=0.797853685
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a3=-0.124807893
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a4=0.136821858
cf2py real(kind=8) optional,intent(in),dimension(nF) :: b1=15.6298
cf2py real(kind=8) optional,intent(in),dimension(nF) :: b2=1.0
cf2py real(kind=8) optional,intent(in) :: rho=1.225
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CI=4.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CO=25.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: CT_idle=0.053
cf2py real(kind=8) intent(out),depend(nF,nF,n),dimension(nF,n) :: P,T,U
      ! internal variables
      integer :: i

      do i=1,nF
        call gclm_s(n,nP,nCT,x_g,y_g,z_g,DT,P_c,CT_c,WS(i,:),WD(i,:),
     &       TI(i,:),a1(i),a2(i),a3(i),a4(i),b1(i),b2(i),rho,WS_CI,
     &       WS_CO,CT_idle,P(i,:),T(i,:),U(i,:))
      end do

      end subroutine gclm

c ----------------------------------------------------------------------
c gclm_av(x,y,z,DT,P_c,CT_c,WS,WD,TI,AV)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES with wt available
c Computes the WindFarm flow and Power using G. C. Larsen model:
c Larsen, G. C., and P. E. Re'thore'. A simple stationary semi-analytical
c wake model. Technical Report Risoe, 2009.
c
c Inputs
c ----------
c x_g (array): Distance between turbines in the global coordinates
c y_g (array): Distance between turbines in the global coordinates
c z_g (array): Distance between turbines in the global coordinates
c DT (array): Turbines diameter
c P_c (array): Power curves
c CT_c (array): Thrust coefficient curves
c WS (array): Undisturbed rotor averaged (equivalent) wind speed at hub
c             height [m/s] for each wind turbine x flow [nF,n]
c WD (array): Undisturbed wind direction at hub height [deg.] for each
c             turbine x flow [nF,n]
c             Meteorological coordinates (N=0,E=90,S=180,W=270)
c TI (array): Turbulence intensity [-] for each turbine x flow [nF,n]
c AV (array): Wind turbine available per flow [nF,n]
c
c rho (float): Air density at which the power curve is valid [kg/m^3]
c WS_CI (array): Cut in wind speed [m/s] for each turbine
c WS_CO (array): Cut out wind speed [m/s] for each turbine
c CT_idle (array): Thrust coefficient at rest [-] for each turbine
c
c Outputs
c ----------
c P (array): Power production of the wind turbines (nWT,1) [W]
c T (array): Thrust force of the wind turbines (nWT,1) [N]
c U (array): Rotor averaged (equivalent) Wind speed at hub height
c            (nWT,1) [m/s]
      subroutine gclm_av(n,nP,nCT,nF,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,TI,
     &AV,a1,a2,a3,a4,b1,b2,rho,WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT,nF,AV(nf,n)
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),rho,WS_CI(n),WS_CO(n),CT_idle(n)
      real(kind=8),dimension(nF) :: a1,a2,a3,a4,b1,b2
      real(kind=8) :: WS(nF,n),WD(nF,n),TI(nF,n),P(nF,n),T(nF,n),U(nF,n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py integer intent(hide),depend(WS) :: nF = size(WS,1)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in),dimension(nF,n) :: WS,WD,TI
cf2py integer intent(in),dimension(nF,n) :: AV
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a1=0.435449861
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a2=0.797853685
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a3=-0.124807893
cf2py real(kind=8) optional,intent(in),dimension(nF) :: a4=0.136821858
cf2py real(kind=8) optional,intent(in),dimension(nF) :: b1=15.6298
cf2py real(kind=8) optional,intent(in),dimension(nF) :: b2=1.0
cf2py real(kind=8) optional,intent(in) :: rho=1.225
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CI=4.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CO=25.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: CT_idle=0.053
cf2py real(kind=8) intent(out),depend(nF,n),dimension(nF,n) :: P,T,U
      ! internal variables
      integer :: i,j
      real(kind=8) :: CT_c_AV(n,nCT,2), P_c_AV(n,nCT,2)

      do i=1,nF
        CT_c_AV = CT_c
        P_c_AV  = P_c
        ! Re-defines the trust curve for non available turbines
        do j=1,n
          if (AV(i,j)==0) then
            CT_c_AV(j,:,2) = CT_idle(j)
            P_c_AV(j,:,2) = 0.0d0
          end if
        end do
        call gclm_s(n,nP,nCT,x_g,y_g,z_g,DT,P_c_AV,CT_c_AV,WS(i,:),
     &       WD(i,:),TI(i,:),a1(i),a2(i),a3(i),a4(i),b1(i),b2(i),rho,
     &       WS_CI,WS_CO,CT_idle,P(i,:),T(i,:),U(i,:))
      end do

      end subroutine gclm_av


c ----------------------------------------------------------------------
c Additional routines
c ----------------------------------------------------------------------

c ----------------------------------------------------------------------
c Find the lower location index of a point in an ordered array
c ----------------------------------------------------------------------
      subroutine loc_in_list(n,xa,x,j)
      implicit none
      integer j,n
      real(kind=8) x,xa(n)
cf2py integer intent(hide),depend(xa) :: n = len(xa)
cf2py real(kind=8) intent(in),dimension(n) :: xa
cf2py real(kind=8) intent(in) :: x
cf2py integer intent(out) :: j
      ! internal variables
      integer j_low,j_mid,j_up
      j_low=0
      j_up=n+1
      ! bisecting loop
      do while (j_up-j_low.gt.1)
        j_mid=(j_up+j_low)/2
        if((xa(n).ge.xa(1)).eqv.(x.ge.xa(j_mid)))then
          j_low=j_mid
        else
          j_up=j_mid
        end if
      end do
      ! check for cases
      if(x.eq.xa(1))then
        j=1
      else if(x.eq.xa(n))then
        j=n-1
      else
        j=j_low
      end if
      end subroutine loc_in_list

c ----------------------------------------------------------------------
c Linear interpolation
c ----------------------------------------------------------------------
      subroutine interp_l(xa,ya,n,x,y)
      implicit none
      integer n
      real(kind=8) x,y,xa(n),ya(n)
cf2py integer intent(hide),depend(xa) :: n = len(xa)
cf2py real(kind=8) intent(in),dimension(n) :: xa
cf2py real(kind=8) intent(in),depend(n),dimension(n) :: ya
cf2py real(kind=8) intent(in) :: x
cf2py real(kind=8) intent(out) :: y
      ! internal variables
      integer j
      call loc_in_list(n,xa,x,j)
      y = ((ya(j+1)-ya(j))/(xa(j+1)-xa(j)))*(x-xa(j)) + ya(j)
      end subroutine interp_l

c ----------------------------------------------------------------------
c Finds the index that order an array from lower to higher values
c only for integers
c
c Knuth, D. E. "Sorting and Searching, vol. 3 of The Art of Computer
c Programming, section 6.2. 2." (1997): 430-31.
c ----------------------------------------------------------------------
      subroutine order_id(n,a,id)
      integer n,a(n),id(n)
cf2py integer intent(hide),depend(a) :: n = len(a)
cf2py integer intent(in),dimension(n) :: a
cf2py integer intent(out),depend(n),dimension(n) :: id
      ! internal variables
      integer i,j,inc,v,w
      do i=1,n
        id(i)=i
      end do
      ! Determine starting increment
      inc=1
      do while (inc.le.n)
        inc=3*inc+1
      end do
      ! Partial sorts loop
      do while (inc.gt.1)
        inc=inc/3
        do i=inc+1,n
          v=a(i)
          w=id(i)
          j=i
          do while (a(j-inc).gt.v)
            a(j)=a(j-inc)
            id(j)=id(j-inc)
            j=j-inc
            if(j.le.inc) exit
          end do
          a(j)=v
          id(j)=w
        end do
      end do
      end subroutine order_id