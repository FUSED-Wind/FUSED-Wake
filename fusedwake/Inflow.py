"""An inflow wind model

@moduleauthor:: Juan P. Murcia <jumu@dtu.dk>

"""
import numpy as np

def get_ABL_U(z,Ur,zr,type='Log',**kwargs):
    """Function of undisturbed inflow wind speed

    Inputs
    ----------
    z:  np.array, float
        Heights to evaluate the streamwise wind speed
    Ur: float
        Reference wind speed
    zr: float
        Reference height

    type: str, optional
        type of ABL model.

        type='log': Log-law
        .. math::
            U = \\frac{u^*}{\\kappa} \\ln \\left( \\frac{ z }{ z0 } \\right)

            kwargs
            ------
            z0      Characteristic roughness length of the terrain

        type='pow': Power law
        .. math::
            U = Ur \\left(\\frac{ z }{ zr } \\right)^\\alpha

            kwargs
            ------
            alpha   Shear exponent

        type='MOB': Monin-Obukhov from [1]
        .. math::
            U = \\frac{u^*}{\\kappa} \\ln \\left( \\frac{ z }{ z0 } - \\phi \\left \\frac{z}{L} \\right  f_s \\right)

            kwargs
            ------
            z0      Characteristic roughness length of the terrain
            L       Obukhov length
            zi      BLH, Boundary layer height


    Returns
    -------
    U:  np.array, float
        Axial wind speed [m/s]

    References
    ----------
    [1] A. Penña, T. Mikkelsen, S.-E. Gryning, C.B. Hasager,
        A.N. Hahmann, M. Badger, et al., O shore vertical wind shear,
        DTU Wind Energy-E-Report-0005(EN), Technical University of Denmark,
        2012.

    """
    if  type=='log':
        kappa = 0.4 # Kappa: Von Karman constant

        if  kwargs is not None and 'z0' in kwargs.keys():
            z0 = kwargs['z0']
        else: # Default
            print("Using default LogLaw characteristic roughness ",
                  "length of the terrain: z0 = 0.0002 (Offshore)")
            z0 = 0.0002

        us = Ur * kappa/np.log(zr/z0)  # friction velocity

        return us / kappa * np.log(z/z0)

    elif type=='pow':
        if  kwargs is not None and 'alpha' in kwargs.keys():
            alpha = kwargs['alpha']
        else: # Default
            print("Using default ABL power law shear coefficient: ",
                  "alpha = 0.143 (Offshore)")
            alpha = 0.143

        return Ur * (z/zr)**alpha

    elif type=='MOB':

        if  kwargs is not None:
            # Roughness Length
            if 'z0' in kwargs.keys():
                z0 = kwargs['z0']
            else:
                print("Using default MOB characteristic roughness ",
                "length of the terrain: z0 = 0.0002 (Offshore)")
                z0 = 0.0002

            # Stability
            if 'L' in kwargs.keys():
                L = kwargs['L']
            else:
                print("Using default MOB Monin-Obukhov length: ",
                "L = -1000 (Neutral from unstable asymptote)")
                L = -1000.

            # ABL height
            if 'zi' in kwargs.keys():
                zi = kwargs['zi']
            else:
                print("Using default MOB ABL height: ",
                "zi = 400 [m]")
                zi = 400.

        else: # Default
            print("Using default MOB ABL coefficients: ",
                "z0 = 0.0002 (Offshore), ",
                "L = -1000 (Neutral from unstable asymptote), ",
                "zi = 400 [m]")
            z0 = 0.0002
            L = -1000.
            zi = 400.
            alpha = 0.143

        if  L>0: #   Stable atmospheric conditions
            phi_m = -4.7*z/L
            phi_m_r = -4.7*zr/L
        else:   #L<0 Unstable atmospheric conditions
            x = (1. - 12.*z/L)**(1./3.)
            phi_m = (3./2.)*np.log( (1. + x + x**2.)/3. ) - \
                    np.sqrt(3.)*np.arctan( (2.*x+1)/np.sqrt(3.) ) + \
                    np.pi/np.sqrt(3.)

            x_r = (1. - 12.*zr/L)**(1./3.)
            phi_m_r = (3./2.)*np.log( (1. + x_r + x_r**2.)/3. ) - \
                    np.sqrt(3.)*np.arctan( (2.*x_r+1)/np.sqrt(3.) ) + \
                    np.pi/np.sqrt(3.)

        return Ur * np.log(z/z0 - phi_m)/np.log(zr/z0 - phi_m_r)

def RotorAvg(f,H,R,dep='z',**kwargs):
    """Rotor averaging function

    .. math::
        Feq = \\int_0^{2\\pi} \\int_0^R f(r_i,\\theta_i) r  dr d\\theta
        Feq = sum w_i f(r_i,\\theta_i)

    Gaussian quadrature using:
    .. math::
        \\theta \sim \\text{Uniform}(0,2\\pi)
        r \sim \\text{Triangular}(0,R) = r/C = c*2/R**2.


    Inputs
    ----------
    f:  python function
        function to be rotor averaged over multiple rotors
    H:  array
        Multiple rotors hub height
    R:  array
        Multiple rotors radii

    dep: str, optional
        type of function dependency

        type='z':
            f is only a function of the height

        type='r':
            f is an axis-symmetrical function of the radius

        type='yz'
            f is a function of the height and horizontal position

    Returns
    -------
    Ueq:  np.array, float
          Rotor averaged axial wind speed [m/s]

    """
    # New improved quadrature rule for wake deficit rotor averaging
    node_R, node_th, weight = np.array([[ 0.26349922998554242692 ,  4.79436403870179805864 ,  0.00579798753740115753 ],
    [ 0.26349922998554242692 ,  5.13630491629471475079 ,  0.01299684397858970851 ],
    [ 0.26349922998554242692 ,  5.71955352542765460555 ,  0.01905256317618122044 ],
    [ 0.26349922998554242692 ,  0.20924454049880022999 ,  0.02341643323656225281 ],
    [ 0.26349922998554242692 ,  1.10309379714216659885 ,  0.02569988335562909190 ],
    [ 0.26349922998554242692 ,  2.03849885644762496284 ,  0.02569988335562912660 ],
    [ 0.26349922998554242692 ,  2.93234811309099407950 ,  0.02341643323656214179 ],
    [ 0.26349922998554242692 ,  3.70522443534172518653 ,  0.01905256317618119616 ],
    [ 0.26349922998554242692 ,  4.28847304447466459720 ,  0.01299684397858971198 ],
    [ 0.26349922998554242692 ,  4.63041392206758217753 ,  0.00579798753740114539 ],
    [ 0.57446451431535072718 ,  4.79436403870179805864 ,  0.01086984853977092380 ],
    [ 0.57446451431535072718 ,  5.13630491629471475079 ,  0.02436599330905551281 ],
    [ 0.57446451431535072718 ,  5.71955352542765460555 ,  0.03571902745281423097 ],
    [ 0.57446451431535072718 ,  0.20924454049880022999 ,  0.04390024659093685194 ],
    [ 0.57446451431535072718 ,  1.10309379714216659885 ,  0.04818117282305908744 ],
    [ 0.57446451431535072718 ,  2.03849885644762496284 ,  0.04818117282305915683 ],
    [ 0.57446451431535072718 ,  2.93234811309099407950 ,  0.04390024659093664378 ],
    [ 0.57446451431535072718 ,  3.70522443534172518653 ,  0.03571902745281418240 ],
    [ 0.57446451431535072718 ,  4.28847304447466459720 ,  0.02436599330905552321 ],
    [ 0.57446451431535072718 ,  4.63041392206758217753 ,  0.01086984853977089951 ],
    [ 0.81852948743000586429 ,  4.79436403870179805864 ,  0.01086984853977090992 ],
    [ 0.81852948743000586429 ,  5.13630491629471475079 ,  0.02436599330905548505 ],
    [ 0.81852948743000586429 ,  5.71955352542765460555 ,  0.03571902745281418934 ],
    [ 0.81852948743000586429 ,  0.20924454049880022999 ,  0.04390024659093679643 ],
    [ 0.81852948743000586429 ,  1.10309379714216659885 ,  0.04818117282305903193 ],
    [ 0.81852948743000586429 ,  2.03849885644762496284 ,  0.04818117282305909438 ],
    [ 0.81852948743000586429 ,  2.93234811309099407950 ,  0.04390024659093658826 ],
    [ 0.81852948743000586429 ,  3.70522443534172518653 ,  0.03571902745281413383 ],
    [ 0.81852948743000586429 ,  4.28847304447466459720 ,  0.02436599330905549199 ],
    [ 0.81852948743000586429 ,  4.63041392206758217753 ,  0.01086984853977088737 ],
    [ 0.96465960618086743494 ,  4.79436403870179805864 ,  0.00579798753740116100 ],
    [ 0.96465960618086743494 ,  5.13630491629471475079 ,  0.01299684397858971545 ],
    [ 0.96465960618086743494 ,  5.71955352542765460555 ,  0.01905256317618123432 ],
    [ 0.96465960618086743494 ,  0.20924454049880022999 ,  0.02341643323656226669 ],
    [ 0.96465960618086743494 ,  1.10309379714216659885 ,  0.02569988335562910925 ],
    [ 0.96465960618086743494 ,  2.03849885644762496284 ,  0.02569988335562914394 ],
    [ 0.96465960618086743494 ,  2.93234811309099407950 ,  0.02341643323656215567 ],
    [ 0.96465960618086743494 ,  3.70522443534172518653 ,  0.01905256317618120657 ],
    [ 0.96465960618086743494 ,  4.28847304447466459720 ,  0.01299684397858972065 ],
    [ 0.96465960618086743494 ,  4.63041392206758217753 ,  0.00579798753740114886 ]]).T

    H_msh, node_R_msh = np.meshgrid(H,node_R)
    R_msh, node_th_msh = np.meshgrid(R,node_th)
    _, weight_msh = np.meshgrid(H,weight)

    if  dep=='z':
        ze = H_msh + R_msh*node_R_msh*np.sin(node_th_msh)
        f_msh = f(ze,**kwargs)

    elif dep=='r':
        ye = R_msh*node_R_msh*np.cos(node_th_msh)
        ze = H_msh + R_msh*node_R_msh*np.sin(node_th_msh)
        re = np.sqrt( ye**2. + ze**2. )

        f_msh = f(re,**kwargs)

    elif dep=='yz':
        ye = R_msh*node_R_msh*np.cos(node_th_msh)
        ze = H_msh + R_msh*node_R_msh*np.sin(node_th_msh)

        f_msh = f(ye,ze,**kwargs)

    elif dep=='xyz':
        xe, weight_msh = np.meshgrid(x,weight)
        ye = R_msh*node_R_msh*np.cos(node_th_msh)
        ze = H_msh + R_msh*node_R_msh*np.sin(node_th_msh)

        f_msh = f(xe,ye,ze,**kwargs)

    return np.sum(weight_msh*f_msh,axis=0)
