c ----------------------------------------------------------------------
c
c N.O.Jensen wake model
c by Juan P. Murcia <jumu@dtu.dk>
c
c ----------------------------------------------------------------------

c ----------------------------------------------------------------------
c get_RW(x,D,CT,kj)
c ----------------------------------------------------------------------
c Computes the wake radius at a location
c
c Inputs
c ----------
c x (float): Distance between turbines in the stream-wise direction
c D (float): Wind turbine diameter
c CT (float): Outputs WindTurbine object's thrust coefficient
c kj (float): Wake (linear) expansion coefficient
c
c Outputs
c ----------
c Rw (float): Wake radius at a location
      subroutine get_RW(n,x,D,CT,RW,kj)

      implicit none
      integer :: n
      real(kind=8) :: x(n),D,CT,RW(n),kj
cf2py integer intent(hide),depend(x) :: n=len(x)
cf2py real(kind=8) intent(in),dimension(n) :: x
cf2py real(kind=8) intent(in) :: D,CT
cf2py real(kind=8) optional,intent(in) :: kj = 0.050
cf2py real(kind=8) intent(out), depend(n),dimension(n) :: RW
      ! internal variables
      real(kind=8), parameter :: pi=3.1415926535897932384626433832795d0
      integer :: i
      real(kind=8) :: Area,a,k

      Area=pi*D*D/4.0d0
      ! axial velocity deficit at rotor disc
      a=(1.0d0-(sqrt(1.0d0-CT)))/2.0d0
      ! Near wake expansion ratio: k = D0/D from momentum balance
      k=1.0d0!sqrt((1.0d0-a)/(1.0d0-2.0d0*a))

      do i = 1, n
        if (x(i)<= 0) then !3.0d0*D ) then !
          ! Null wake radius for locations in front of the rotor disc
          RW(i) = 0.0d0
        else
          ! Linear wake expansion x0 = 3D
          !RW(i)=0.5d0*k*D + kj*( x(i) - 3.0d0*D )
          RW(i)=0.5d0*k*D + kj*( x(i) )
        end if
      end do
      end subroutine get_RW

c ----------------------------------------------------------------------
c get_dU(x,r,D,CT,kj)
c ----------------------------------------------------------------------
c Computes the wake velocity deficit at a location
c
c Inputs
c ----------
c x (float): Distance between turbines in the stream-wise direction
c r (array): Radial distance between the turbine and the location
c D (float): Wind turbine diameter
c CT (float): Outputs WindTurbine object's thrust coefficient
c kj (float): Wake (linear) expansion coefficient
c
c Outputs
c ----------
c dU (float): Wake velocity deficit at a location normalized
c             by rotor averaged (equivalent) inflow velocity
      subroutine get_dU(n,x,r,D,CT,kj,dU)

      implicit none
      integer :: n
      real(kind=8) :: x(n),r(n),D,CT,kj,dU(n)
cf2py integer intent(hide),depend(x) :: n=len(x)
cf2py real(kind=8) intent(in),dimension(n) :: x
cf2py real(kind=8) intent(in),depend(n),dimension(n) :: r
cf2py real(kind=8) intent(in) :: D,CT
cf2py real(kind=8) optional,intent(in) :: kj = 0.050
cf2py real(kind=8) intent(out),depend(n),dimension(n) :: dU
      ! internal variables
      real(kind=8), parameter :: pi=3.1415926535897932384626433832795d0
      integer :: i
      real(kind=8) :: a
      real(kind=8), dimension(n) :: RW

      ! axial velocity deficit at rotor disc
      a=(1.0d0-(sqrt(1.0d0-CT)))/2.0d0
      ! wake radius at x locations
      call get_RW(n,x,D,CT,RW,kj)

      do i = 1, n
        if ((x(i) <= 0).or.(r(i) >= RW(i))) then
          ! Null deficit for locations in front of the rotor disc
          ! and outside the wake radius
          dU(i) = 0.0d0
        else
          dU(i)=-(2.0d0*a*(D**(2.0d0)))/((2.0d0*RW(i))**(2.0d0))
        end if
      end do
      end subroutine get_dU

c ----------------------------------------------------------------------
c get_dUeq(x,y,z,DT,D,CT,kj)
c ----------------------------------------------------------------------
c Computes the rotor averaged (equivalent) wake velocity deficit at
c different turbine locations and diameters.
c Takes into account partial wake conditions as presented in:
c Wan, C., Wang, J., Yang, G., Gu, H., & Zhang, X. (2012).
c Wind farm micro-siting by Gaussian particle swarm optimization with
c local search strategy. Renewable Energy, 48, 276-286.
c
c Inputs
c ----------
c x (array): Distance between turbines in the stream-wise direction
c y (array): Distance between turbines in the cross-flow direction
c z (array): Distance between turbines in the vertical direction
c DT (array): Wake operating turbines diameter
c D (float): Wake generating turbine diameter
c CT (float): Outputs WindTurbine object's thrust coefficient
c kj (float): Wake (linear) expansion coefficient
c
c Outputs
c ----------
c dUeq (float): Wake velocity deficit at a location normalized by
c               inflow velocity
      subroutine get_dUeq(n,x,y,z,DT,D,CT,kj,dUeq)

      implicit none
      integer :: n
      real(kind=8) :: x(n),y(n),z(n),DT(n),D,CT,kj,dUeq(n)
cf2py integer intent(hide),depend(x) :: n = len(x)
cf2py real(kind=8) intent(in),dimension(n) :: x
cf2py real(kind=8) intent(in),depend(n),dimension(n) :: y,z,DT
cf2py real(kind=8) intent(in) :: D,CT
cf2py real(kind=8) optional,intent(in) :: kj = 0.050
cf2py real(kind=8) intent(out),depend(n),dimension(n) :: dUeq
      ! internal variables
      real(kind=8), parameter :: pi=3.1415926535897932384626433832795d0
      integer :: i
      real(kind=8) :: Area,alpha1,alpha2,q,tm1,tm2,tm3,tm4,tm5,tm6
      real(kind=8), dimension(1) :: x_e,r_e,dU
      real(kind=8), dimension(n) :: r,RT,RW

      Area=pi*D*D/4.0d0
      RT = DT/2.0d0
      ! Location of the turbines in wake coordinates
      r  = (y**(2.0d0) + z**(2.0d0))**(0.5d0)
      ! wake radius at turbine x locations
      call get_RW(n,x,D,CT,RW,kj)
      do i = 1, n
        if (x(i) > 0.0) then
          if (((x(i)==0.0d0).and.(r(i)==0.0d0)).or.(RW(i)==0.0d0)) then
            dUeq(i)=0.0d0
          ! turbine is totally inside the wake
          else if (r(i) <= (RW(i)-RT(i))) then
            x_e = x(i)
            r_e = r(i)
            call get_dU(1,x_e,r_e,D,CT,kj,dU)
            dUeq(i)=dU(1)
          ! Partial wakes
          else if ((r(i)>(RW(i)-RT(i))).and.(r(i)<(RW(i)+RT(i)))) then
            x_e = x(i)
            r_e = 0.5d0*(RW(i)+r(i)-RT(i))
            call get_dU(1,x_e,r_e,D,CT,kj,dU)
            tm1 = RT(i)**(2.0d0)+r(i)**(2.0d0)-RW(i)**(2.0d0)
            tm2 = tm1/(2.0d0*RT(i)*r(i))
            alpha1 = 2.0d0*acos(tm2)

            tm3 = RW(i)**(2.0d0)+r(i)**(2.0d0)-RT(i)**(2.0d0)
            tm4 = tm3/(2.0d0*RW(i)*r(i))
            alpha2 = 2.0d0*acos(tm4)

            tm5 = 0.5d0*((RT(i))**(2.0d0))*(alpha1 - sin(alpha1))
            tm6 = 0.5d0*((RW(i))**(2.0d0))*(alpha2 - sin(alpha2))
            q = (tm5 + tm6)/Area

            dUeq(i)=dU(1)*q
          else
            dUeq(i)=0.0d0
          end if
        else
          dUeq(i)=0.0d0
        end if
      end do

      end subroutine get_dUeq


c ----------------------------------------------------------------------
c noj_s(x,y,z,DT,P_c,CT_c,WS,kj)
c ----------------------------------------------------------------------
c SINGLE FLOW CASE
c Computes the WindFarm flow and Power using N. O. Jensen model:
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
c kj (float): Wake (linear) expansion coefficient
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
      subroutine noj_s(n,nP,nCT,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,kj,
     &rho,WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),WS,WD,kj
      real(kind=8) :: rho,WS_CI(n),WS_CO(n),CT_idle(n),P(n),T(n),U(n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in) :: WS,WD
cf2py real(kind=8) optional,intent(in) :: kj = 0.050
cf2py real(kind=8) optional,intent(in) :: rho = 1.225
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CI = 4.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CO = 25.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: CT_idle = 0.053
cf2py real(kind=8) intent(out),depend(n),dimension(n) :: P,T,U
      ! internal variables
      real(kind=8), parameter :: pi=3.1415926535897932384626433832795d0
      integer :: i,j,k,nDownstream(n),idT(n)
      real(kind=8) :: x(n),y(n),z(n),D,CT,angle,dUeq(n),dUsq(n)
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
      dUsq = 0d0
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
        call get_dUeq(n,x,y,z,DT,D,CT,kj,dUeq)
        !
        !Wake deficits are normalized by the local velocity
        !Linear sum wake deficits superposition
        !dUsq = (U(i)*dUeq)**(2.0d0)
        !U = U - dUsq**(0.5d0)
        !
        !All deficits are normalized by the inflow velocity
        !Squared root of the sum of squares wake deficits superposition
        dUsq = dUsq + (WS*dUeq)**(2.0d0)
        U = WS - dUsq**(0.5d0)
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

      end subroutine noj_s

c ----------------------------------------------------------------------
c noj(x,y,z,DT,P_c,CT_c,WS,WD,kj)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES
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
c kj (float): Wake (linear) expansion coefficient
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
      subroutine noj(n,nP,nCT,nF,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,kj,
     &rho,WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT,nF
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),rho,WS_CI(n),WS_CO(n),CT_idle(n)
      real(kind=8),dimension(nF) :: WS,WD,kj
      real(kind=8) :: P(nF,n),T(nF,n),U(nF,n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py integer intent(hide),depend(WS) :: nF = len(WS)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in),dimension(nF) :: WS,WD
cf2py real(kind=8) optional,intent(in),dimension(nF)::kj = 0.050
cf2py real(kind=8) optional,intent(in) :: rho = 1.225
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CI = 4.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CO = 25.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: CT_idle = 0.053
cf2py real(kind=8) intent(out),depend(n),dimension(nF,n) :: P,T,U
      ! internal variables
      integer :: i

      do i=1,nF
        call noj_s(n,nP,nCT,x_g,y_g,z_g,DT,P_c,CT_c,WS(i),WD(i),kj(i),
     &            rho,WS_CI,WS_CO,CT_idle,P(i,:),T(i,:),U(i,:))
      end do

      end subroutine noj

c ----------------------------------------------------------------------
c noj_av(x,y,z,DT,P_c,CT_c,WS,WD,kj,AV)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES with wt available
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
c kj (float): Wake (linear) expansion coefficient
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
      subroutine noj_av(n,nP,nCT,nF,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,kj,AV,
     &rho,WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT,nF,AV(nf,n)
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),rho,WS_CI(n),WS_CO(n),CT_idle(n)
      real(kind=8),dimension(nF) :: WS,WD,kj
      real(kind=8) :: P(nF,n),T(nF,n),U(nF,n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py integer intent(hide),depend(WS) :: nF = len(WS)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in),dimension(nF) :: WS,WD
cf2py integer intent(in),dimension(nF,n) :: AV
cf2py real(kind=8) optional,intent(in),dimension(nF)::kj = 0.050
cf2py real(kind=8) optional,intent(in) :: rho = 1.225
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CI = 4.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CO = 25.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: CT_idle = 0.053
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
        call noj_s(n,nP,nCT,x_g,y_g,z_g,DT,P_c_AV,CT_c_AV,WS(i),WD(i),
     &            kj(i),rho,WS_CI,WS_CO,CT_idle,P(i,:),T(i,:),U(i,:))
      end do

      end subroutine noj_av

c ----------------------------------------------------------------------
c noj_GA(x,y,z,DT,P_c,CT_c,WS,WD,STD_WD,Nga)
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
      subroutine noj_GA(n,nP,nCT,nF,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,
     &STD_WD,Nga,kj,rho,WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT,nF,Nga
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),rho,WS_CI(n),WS_CO(n),CT_idle(n)
      real(kind=8),dimension(nF) :: WS,WD,STD_WD,kj
      real(kind=8) :: P(nF,n),T(nF,n),U(nF,n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py integer intent(hide),depend(WS) :: nF = len(WS)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in),dimension(nF) :: WS,WD,STD_WD
cf2py integer optional,intent(in) :: Nga = 4
cf2py real(kind=8) optional,intent(in),dimension(nF)::kj = 0.050
cf2py real(kind=8) optional,intent(in) :: rho = 1.225
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CI = 4.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CO = 25.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: CT_idle = 0.053
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
          call noj_s(n,nP,nCT,x_g,y_g,z_g,DT,P_c,CT_c,WS(i),WD_aux,
     &      kj(i),rho,WS_CI,WS_CO,CT_idle,P_aux,T_aux,U_aux)
          P(i,:)=P(i,:)+weight(j)*P_aux*(1.0d0/sqrt(pi))
          T(i,:)=T(i,:)+weight(j)*T_aux*(1.0d0/sqrt(pi))
          U(i,:)=U(i,:)+weight(j)*U_aux*(1.0d0/sqrt(pi))
        end do
      end do

      end subroutine noj_GA

c ----------------------------------------------------------------------
c noj_mult_wd(x,y,z,DT,P_c,CT_c,WS,WD,TI)
c ----------------------------------------------------------------------
c MULTIPLE FLOW CASES with individual wind direction for each turbine
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
      subroutine noj_mult_wd(n,nP,nCT,nF,x_g,y_g,z_g,DT,P_c,CT_c,WS,WD,
     &STD_WD,Nga,kj,rho,WS_CI,WS_CO,CT_idle,P,T,U)

      implicit none
      integer :: n,nP,nCT,nF,Nga
      real(kind=8) :: x_g(n,n),y_g(n,n),z_g(n,n),DT(n),P_c(n,nP,2)
      real(kind=8) :: CT_c(n,nCT,2),rho,WS_CI(n),WS_CO(n),CT_idle(n)
      real(kind=8),dimension(nF) :: WS,kj
      real(kind=8),dimension(nF,n) :: WD, STD_WD
      real(kind=8) :: P(nF,n),T(nF,n),U(nF,n)
cf2py integer intent(hide),depend(DT) :: n = len(DT)
cf2py integer intent(hide),depend(P_c) :: nP = size(P_c,2)
cf2py integer intent(hide),depend(CT_c) :: nCT = size(CT_c,2)
cf2py integer intent(hide),depend(WS) :: nF = len(WS)
cf2py real(kind=8) intent(in),dimension(n) :: DT
cf2py real(kind=8) intent(in),depend(n),dimension(n,n) :: x_g,y_g,z_g
cf2py real(kind=8) intent(in),dimension(n,nP,2) :: P_c
cf2py real(kind=8) intent(in),dimension(n,nCT,2) :: CT_c
cf2py real(kind=8) intent(in),dimension(nF) :: WS
cf2py real(kind=8) intent(in),dimension(nF,n) :: WD
cf2py real(kind=8) optional,intent(in),dimension(nF,n) :: STD_WD = 0.0
cf2py integer optional intent(in) :: Nga = 1
cf2py real(kind=8) optional,intent(in),dimension(nF)::kj = 0.050
cf2py real(kind=8) optional,intent(in) :: rho = 1.225
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CI = 4.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: WS_CO = 25.0
cf2py real(kind=8) optional,intent(in),dimension(n) :: CT_idle = 0.053
cf2py real(kind=8) intent(out),depend(n),dimension(nF,n) :: P,T,U
      ! internal variables
      integer :: i,j
      real(kind=8), dimension(n) :: P_aux,T_aux,U_aux

      do j=1,n
        do i=1,nF
          call noj_GA(n,nP,nCT,1,x_g,y_g,z_g,DT,P_c,CT_c,WS(i),WD(i,j),
     &STD_WD(i,j),Nga,kj(i),rho,WS_CI,WS_CO,CT_idle,P_aux,T_aux,U_aux)
!          call noj_s(n,nP,nCT,x_g,y_g,z_g,DT,P_c,CT_c,WS(i),WD(i,j),
!     &      kj(i),rho,WS_CI,WS_CO,CT_idle,P_aux,T_aux,U_aux)
          P(i,j)=P_aux(j)
          T(i,j)=T_aux(j)
          U(i,j)=U_aux(j)
        end do
      end do

      end subroutine noj_mult_wd

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
