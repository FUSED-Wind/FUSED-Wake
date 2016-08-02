c ----------------------------------------------------------------------
c
c Gaussian wake model applied to offshore wind farms (WindFarm object)
c by Juan P. Murcia <jumu@dtu.dk>
c
c Bastankhah, M., & Porte'-Agel, F. (2014). A new analytical model for
c wind-turbine wakes. Renewable Energy, 70, 116-123.
c
c ----------------------------------------------------------------------


c ----------------------------------------------------------------------
c get_RW(x,D,CT,ks)
c ----------------------------------------------------------------------
c Computes the wake radius at a location
c
c Inputs
c ----------
c x (float): Distance between turbines in the stream-wise direction
c D (float): Wind turbine diameter
c CT (float): Outputs WindTurbine object's thrust coefficient
c ks (float): Wake (linear) expansion coefficient
c
c Outputs
c ----------
c Rw (float): Wake radius at a location
      subroutine get_RW(n,x,D,CT,RW,ks)

      implicit none
      integer :: n
      real(kind=8) :: x(n),D,CT,RW(n),ks
cf2py integer intent(hide),depend(x) :: n=len(x)
cf2py real(kind=8) intent(in),dimension(n) :: x
cf2py real(kind=8) intent(in) :: D,CT
cf2py real(kind=8) optional,intent(in) :: ks = 0.030 
cf2py real(kind=8) intent(out), depend(n),dimension(n) :: RW
      ! internal variables
      real(kind=8), parameter :: pi=3.1415926535897932384626433832795d0
      integer :: i
      real(kind=8) :: Area,a,k

      Area=pi*D*D/4.0d0
      ! axial velocity deficit at rotor disc
      a=(1.0d0-(sqrt(1.0d0-CT)))/2.0d0
      ! Near wake expansion ratio: k = D0/D from momentum balance
      k=sqrt((1.0d0-a)/(1.0d0-2.0d0*a))

      do i = 1, n
        ! Linear wake expansion:
        ! RW is defined as the sigma of the Gaussian wake deficit
        RW(i)=0.2d0*k*D + ks*x(i)
      end do
      end subroutine get_RW

c ----------------------------------------------------------------------
c get_dU(x,r,D,CT,ks)
c ----------------------------------------------------------------------
c Computes the wake velocity deficit at a location
c Bastankhah, M., & Porte'-Agel, F. (2014). A new analytical model for
c wind-turbine wakes. Renewable Energy, 70, 116-123.
c
c Inputs
c ----------
c x (float): Distance between turbines in the stream-wise direction
c r (array): Radial distance between the turbine and the location
c D (float): Wind turbine diameter
c ks (float): Wake (linear) expansion coefficient [-]
c CT (float): Outputs WindTurbine object's thrust coefficient
c
c Outputs
c ----------
c dU (float): Wake velocity deficit at a location normalized
c             by rotor averaged (equivalent) inflow velocity
      subroutine get_dU(n,x,r,D,CT,ks,dU)

      implicit none
      integer :: n
      real(kind=8) :: x(n),r(n),D,CT,ks,dU(n)
cf2py integer intent(hide),depend(x) :: n=len(x)
cf2py real(kind=8) intent(in),dimension(n) :: x
cf2py real(kind=8) intent(in),depend(n),dimension(n) :: r
cf2py real(kind=8) intent(in) :: D,CT,ks
cf2py real(kind=8) intent(out),depend(n),dimension(n) :: dU
      ! internal variables
      integer :: i
      real(kind=8), dimension(n) :: RW, tm0, tm10,tm20

      call get_RW(n,x,D,CT,RW,ks)

      do i = 1, n
        tm0(i) = (RW(i)/D)**(2.0d0)
        if ( (x(i)<=2.d0*D).or.(1.0d0<=(CT/(8.0d0*tm0(i)))) ) then
          ! Null wake radius for locations in front of the rotor disc
          dU(i) = 0.0d0
        else
          tm10(i)=1.0d0 - ( (1.0d0 - (CT/(8.0d0*tm0(i))))**(0.5d0) )
          tm20(i)=exp(-(1.0d0/(2.0d0*tm0(i)))*((r(i)/D)**(2.0d0)))
          dU(i)=-tm10(i)*tm20(i)
        end if
      end do

      end subroutine get_dU

c ----------------------------------------------------------------------
c get_dUeq(x,y,z,DT,D,CT,ks)
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
c ks (float): Wake (linear) expansion coefficient [-]
c CT (float): Outputs WindTurbine object's thrust coefficient
c Ng (int): Polynomial order for Gauss-Legendre quadrature integration
c           in both radial and angular positions
c
c Outputs
c ----------
c dUeq (float): Wake velocity deficit at a location normalized by
c               inflow velocity
      subroutine get_dUeq(n,x,y,z,DT,D,CT,ks,Ng,dUeq)

      implicit none
      integer :: n,Ng
      real(kind=8) :: x(n),y(n),z(n),DT(n),D,CT,ks
      real(kind=8) :: dUeq(n)
cf2py integer intent(hide),depend(x) :: n = len(x)
cf2py real(kind=8) intent(in),dimension(n) :: x
cf2py real(kind=8) intent(in),depend(n),dimension(n) :: y,z,DT
cf2py real(kind=8) intent(in) :: D,CT,ks
cf2py integer optional intent(in) :: Ng = 4
cf2py real(kind=8) intent(out),depend(n),dimension(n) :: dUeq
      ! internal variables
      real(kind=8), parameter :: pi=3.1415926535897932384626433832795d0
      integer :: i,j,k
      real(kind=8) :: tm1,tm2,tm3,tm4
      real(kind=8), dimension(Ng) :: root,weight,r_pr,th_pr
      real(kind=8), dimension(1) :: x_e,r_e,dU
      real(kind=8), dimension(n) :: RT,r_R,th_R

      RT = DT/2.0d0
      ! Gauss-Legendre quadrature points and weights
      select case (ng)
       case ( 4 )
          root(1) = -0.3399810435848563d0
          root(2) = 0.3399810435848563d0
          root(3) = -0.8611363115940526d0
          root(4) = 0.8611363115940526d0

          weight(1) = 0.6521451548625461d0
          weight(2) = 0.6521451548625461d0
          weight(3) = 0.3478548451374538d0
          weight(4) = 0.3478548451374538d0
       case (5)
          root(1) = 0.0000000000000000d0
          root(2) = -0.5384693101056831d0
          root(3) = 0.5384693101056831d0
          root(4) = -0.9061798459386640d0
          root(5) = 0.9061798459386640d0

          weight(1) = 0.5688888888888889d0
          weight(2) = 0.4786286704993665d0
          weight(3) = 0.4786286704993665d0
          weight(4) = 0.2369268850561891d0
          weight(5) = 0.2369268850561891d0
       case (6)
          root(1) = 0.6612093864662645d0
          root(2) = -0.6612093864662645d0
          root(3) = -0.2386191860831969d0
          root(4) = 0.2386191860831969d0
          root(5) = -0.9324695142031521d0
          root(6) = 0.9324695142031521d0

          weight(1) = 0.3607615730481386d0
          weight(2) = 0.3607615730481386d0
          weight(3) = 0.4679139345726910d0
          weight(4) = 0.4679139345726910d0
          weight(5) = 0.1713244923791704d0
          weight(6) = 0.1713244923791704d0
       case (7)
          root(1) = 0.0000000000000000d0
          root(2) = 0.4058451513773972d0
          root(3) = -0.4058451513773972d0
          root(4) = -0.7415311855993945d0
          root(5) = 0.7415311855993945d0
          root(6) = -0.9491079123427585d0
          root(7) = 0.9491079123427585d0

          weight(1) = 0.4179591836734694d0
          weight(2) = 0.3818300505051189d0
          weight(3) = 0.3818300505051189d0
          weight(4) = 0.2797053914892766d0
          weight(5) = 0.2797053914892766d0
          weight(6) = 0.1294849661688697d0
          weight(7) = 0.1294849661688697d0
        case ( 8 )
          root(1) = -0.960289856497536d0
          root(2) = -0.796666477413627d0
          root(3) = -0.525532409916329d0
          root(4) = -0.183434642495650d0
          root(5) = 0.183434642495650d0
          root(6) = 0.525532409916329d0
          root(7) = 0.796666477413627d0
          root(8) = 0.960289856497536d0

          weight(1) = 0.101228536290374d0
          weight(2) = 0.222381034453375d0
          weight(3) = 0.313706645877888d0
          weight(4) = 0.362683783378363d0
          weight(5) = 0.362683783378363d0
          weight(6) = 0.313706645877888d0
          weight(7) = 0.222381034453375d0
          weight(8) = 0.101228536290374d0
      end select

      ! Location of the turbines in wake coordinates
      r_R  = (y**(2.0d0) + z**(2.0d0))**(0.5d0)
      th_R = modulo(atan2(z,y),2.0d0*pi)

      do i = 1, n
        dUeq(i) = 0.0d0
        if (x(i) > 0.0) then
          ! Location of evaluation points in the local rotor coordinates
          r_pr  = RT(i)*(root+1d0)/2.0d0!uniform in [0, RT]
          !th_pr = pi*(root+1d0)        !uniform in [0,2*pi]
          th_pr = pi*(root+1d0)-pi/2.d0 !uniform in [-pi/2,3/2*pi]
          ! Location of evaluation points in wake coordinates
          ! Evaluation of wake and sum of quadrature
          do j = 1, Ng
            do k = 1, Ng
              x_e = x(i)
              tm1 = (r_R(i))**(2.0d0)
              tm2 = (r_pr(k))**(2.0d0)
              tm3 = 2d0*r_R(i)*r_pr(k)*cos(th_R(i) - th_pr(j))
              r_e = sqrt( tm1+tm2+tm3)
              call get_dU(1,x_e,r_e,D,CT,ks,dU)
              tm4 = weight(j)*weight(k)*dU(1)*(root(k)+1d0)/4d0
              dUeq(i)=dUeq(i)+tm4
            end do
          end do
        end if
      end do

      end subroutine get_dUeq


c ----------------------------------------------------------------------
c gau_s(x,y,z,DT,P_c,CT_c,WS,ks)
c ----------------------------------------------------------------------
c SINGLE FLOW CASE
c Bastankhah, M., & Porte'-Agel, F. (2014). A new analytical model for
c wind-turbine wakes. Renewable Energy, 70, 116-123.
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
c ks (float): Wake (linear) expansion coefficient [-]
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
c            (nWT,1) [m/s]
      subroutine gau_s(n,nP,nCT,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,ks,Ng,rho,
     &WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT,Ng
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),WS,WD,ks
      real(kind=8) :: rho,WS_CI(n),WS_CO(n),CT_idle(n),P(n),T(n),U(n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in) :: WS,WD,ks
cf2py integer optional intent(in) :: Ng = 4
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
        call get_dUeq(n,x,y,z,DT,D,CT,ks,Ng,dUeq)
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

      end subroutine gau_s

c ----------------------------------------------------------------------
c gau(x,y,z,DT,P_c,CT_c,WS,WD,ks)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES
c Bastankhah, M., & Porte'-Agel, F. (2014). A new analytical model for
c wind-turbine wakes. Renewable Energy, 70, 116-123.
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
c ks (array): Linear wake expansion [-]
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
      subroutine gau(n,nP,nCT,nF,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,ks,
     &Ng,rho,WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT,nF,Ng
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),rho,WS_CI(n),WS_CO(n),CT_idle(n)
      real(kind=8),dimension(nF) :: WS,WD,ks
      real(kind=8) :: P(nF,n),T(nF,n),U(nF,n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py integer intent(hide),depend(WS) :: nF = len(WS)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in),dimension(nF) :: WS,WD,ks
cf2py integer optional intent(in) :: Ng = 4
cf2py real(kind=8) optional,intent(in) :: rho=1.225
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CI=4.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CO=25.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: CT_idle=0.053
cf2py real(kind=8) intent(out),depend(n),dimension(nF,n) :: P,T,U
      ! internal variables
      integer :: i

      do i=1,nF
        call gau_s(n,nP,nCT,x_g,y_g,z_g,DT,P_c,CT_c,WS(i),WD(i),ks(i),
     &       Ng,rho,WS_CI,WS_CO,CT_idle,P(i,:),T(i,:),U(i,:))
      end do

      end subroutine gau

c ----------------------------------------------------------------------
c gau_av(x,y,z,DT,P_c,CT_c,WS,WD,ks,AV)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES with wt available
c Bastankhah, M., & Porte'-Agel, F. (2014). A new analytical model for
c wind-turbine wakes. Renewable Energy, 70, 116-123.
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
c ks (array): Linear wake expansion [-]
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
      subroutine gau_av(n,nP,nCT,nF,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,ks,AV,
     &Ng,rho,WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT,nF,Ng,AV(nf,n)
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),rho,WS_CI(n),WS_CO(n),CT_idle(n)
      real(kind=8),dimension(nF) :: WS,WD,ks
      real(kind=8) :: P(nF,n),T(nF,n),U(nF,n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py integer intent(hide),depend(WS) :: nF = len(WS)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in),dimension(nF) :: WS,WD,ks
cf2py integer intent(in),dimension(nF,n) :: AV
cf2py integer optional intent(in) :: Ng = 4
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
        call gau_s(n,nP,nCT,x_g,y_g,z_g,DT,P_c_AV,CT_c_AV,WS(i),WD(i),
     &       ks(i),Ng,rho,WS_CI,WS_CO,CT_idle,P(i,:),T(i,:),U(i,:))
      end do

      end subroutine gau_av

c ----------------------------------------------------------------------
c gau_GA(x,y,z,DT,P_c,CT_c,WS,WD,ks,STD_WD,Nga)
c ----------------------------------------------------------------------
c GAUSSIAN AVERAGE
c Gauss-Hermite quadrature for normally distributed wind direction
c uncertainty (inside Reynolds averaging time), for a global/unique wind
c direction uncertainty
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
c ks (array): Linear wake expansion [-]
c STD_WD(array): Standard deviation of wind direction uncertainty
c Nga (int): Number of quadrature points for Gaussian averaging
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
      subroutine gau_GA(n,nP,nCT,nF,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,
     &ks,STD_WD,Nga,Ng,rho,WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT,nF,Ng,Nga
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),rho,WS_CI(n),WS_CO(n),CT_idle(n)
      real(kind=8),dimension(nF) :: WS,WD,ks,STD_WD
      real(kind=8) :: P(nF,n),T(nF,n),U(nF,n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py integer intent(hide),depend(WS) :: nF = len(WS)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in),dimension(nF) :: WS,ks,WD,STD_WD
cf2py integer optional,intent(in) :: Nga = 4
cf2py integer optional intent(in) :: Ng = 4
cf2py real(kind=8) optional,intent(in) :: rho=1.225
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CI=4.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CO=25.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: CT_idle=0.053
cf2py real(kind=8) intent(out),depend(n),dimension(nF,n) :: P,T,U
      ! internal variables
      real(kind=8), parameter :: pi=3.1415926535897932384626433832795d0
      integer :: i,j
      real(kind=8) :: WD_aux
      real(kind=8), dimension(n) :: P_aux,T_aux,U_aux
      real(kind=8), dimension(Nga) :: root,weight

      ! Gauss-Hermite quadrature points and weights
      select case (Nga)
        case ( 1 )
          root(1) = 0.0d0
          weight(1) = sqrt(pi)
        case ( 4 )
          root(1) = -1.650680123885785d0
          root(2) = -0.524647623275290d0
          root(3) = 0.524647623275290d0
          root(4) = 1.650680123885785d0

          weight(1) = 0.081312835447245d0
          weight(2) = 0.804914090005513d0
          weight(3) = 0.804914090005513d0
          weight(4) = 0.081312835447245d0
        case ( 5 )
          root(1) = -2.020182870456086d0
          root(2) = -0.958572464613819d0
          root(3) = 0.000000000000000d0
          root(4) = 0.958572464613819d0
          root(5) = 2.020182870456086d0

          weight(1) = 0.019953242059046d0
          weight(2) = 0.393619323152241d0
          weight(3) = 0.945308720482942d0
          weight(4) = 0.393619323152241d0
          weight(5) = 0.019953242059046d0
        case ( 6 )
          root(1) = -2.350604973674492d0
          root(2) = -1.335849074013697d0
          root(3) = -0.436077411927617d0
          root(4) = 0.436077411927617d0
          root(5) = 1.335849074013697d0
          root(6) = 2.350604973674492d0

          weight(1) = 0.004530009905509d0
          weight(2) = 0.157067320322856d0
          weight(3) = 0.724629595224392d0
          weight(4) = 0.724629595224392d0
          weight(5) = 0.157067320322856d0
          weight(6) = 0.004530009905509d0
        case ( 7 )
          root(1) = -2.651961356835233d0
          root(2) = -1.673551628767471d0
          root(3) = -0.816287882858965d0
          root(4) = 0.000000000000000d0
          root(5) = 0.816287882858965d0
          root(6) = 1.673551628767471d0
          root(7) = 2.651961356835233d0

          weight(1) = 0.000971781245100d0
          weight(2) = 0.054515582819127d0
          weight(3) = 0.425607252610128d0
          weight(4) = 0.810264617556807d0
          weight(5) = 0.425607252610128d0
          weight(6) = 0.054515582819127d0
          weight(7) = 0.000971781245100d0
        case ( 8 )
          root(1) = -2.930637420257244d0
          root(2) = -1.981656756695843d0
          root(3) = -1.157193712446780d0
          root(4) = -0.381186990207322d0
          root(5) = 0.381186990207322d0
          root(6) = 1.157193712446780d0
          root(7) = 1.981656756695843d0
          root(8) = 2.930637420257244d0

          weight(1) = 0.000199604072211d0
          weight(2) = 0.017077983007413d0
          weight(3) = 0.207802325814892d0
          weight(4) = 0.661147012558241d0
          weight(5) = 0.661147012558241d0
          weight(6) = 0.207802325814892d0
          weight(7) = 0.017077983007413d0
          weight(8) = 0.000199604072211d0
      end select

      do i=1,nF
        P(i,:)=0.0d0
        T(i,:)=0.0d0
        U(i,:)=0.0d0
        do j=1,Nga
          WD_aux = WD(i)+sqrt(2.0d0)*STD_WD(i)*root(j)
          call gau_s(n,nP,nCT,x_g,y_g,z_g,DT,P_c,CT_c,WS(i),WD_aux,
     &     ks(i),Ng,rho,WS_CI,WS_CO,CT_idle,P_aux,T_aux,U_aux)
          P(i,:)=P(i,:)+weight(j)*P_aux*(1.0d0/sqrt(pi))
          T(i,:)=T(i,:)+weight(j)*T_aux*(1.0d0/sqrt(pi))
          U(i,:)=U(i,:)+weight(j)*U_aux*(1.0d0/sqrt(pi))
        end do
      end do

      end subroutine gau_GA

c ----------------------------------------------------------------------
c gau_mult_wd(x,y,z,DT,P_c,CT_c,WS,WD,ks)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES with individual wind direction for each turbine
c and optional individual Gaussian averaging
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
c ks (array): Linear wake expansion [-]
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
      subroutine gau_mult_wd(n,nP,nCT,nF,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,
     &ks,STD_WD,Nga,Ng,rho,WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT,nF,Ng,Nga
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),rho,WS_CI(n),WS_CO(n),CT_idle(n)
      real(kind=8),dimension(nF) :: WS,ks
      real(kind=8),dimension(nF,n) :: WD,STD_WD
      real(kind=8) :: P(nF,n),T(nF,n),U(nF,n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py integer intent(hide),depend(WS) :: nF = len(WS)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in),dimension(nF) :: WS,ks
cf2py real(kind=8) intent(in),dimension(nF,n) :: WD
cf2py real(kind=8) optional,intent(in),dimension(nF,n) :: STD_WD = 0.0
cf2py integer optional intent(in) :: Nga = 1
cf2py integer optional intent(in) :: Ng = 4
cf2py real(kind=8) optional,intent(in) :: rho=1.225
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CI=4.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CO=25.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: CT_idle=0.053
cf2py real(kind=8) intent(out),depend(n),dimension(nF,n) :: P,T,U
      ! internal variables
      integer :: i,j
      real(kind=8), dimension(n) :: P_aux,T_aux,U_aux

      do j=1,n
        do i=1,nF
          call gau_GA(n,nP,nCT,1,x_g,y_g,z_g,DT,P_c,CT_c,WS(i),WD(i,j),
     &ks(i),STD_WD(i,j),Nga,Ng,rho,
     &WS_CI,WS_CO,CT_idle,P_aux,T_aux,U_aux)

          P(i,j)=P_aux(j)
          T(i,j)=T_aux(j)
          U(i,j)=U_aux(j)
        end do
      end do

      end subroutine gau_mult_wd

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
