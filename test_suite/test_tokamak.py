import numpy as np
from io import StringIO

from .. import tokamak
from .. import _geqdsk

def test_tokamak_interpolations():
    """Test interpolations and derivatives"""
    
    # Define 2D (R,Z) grid
    r1d = np.linspace(1.0, 2.0, 65)
    z1d = np.linspace(-1.0, 1.0, 65)
    r2d, z2d = np.meshgrid(r1d, z1d, indexing='ij')

    # A poloidal flux function
    r0 = 1.5
    z0 = 0.0

    def psi_func(R,Z):
        return np.exp(-((R - r0)**2 + (Z - z0)**2)/0.3**2)

    def dpsi_dr(R, Z):
        "Derivative of psi in R"
        return - (2/0.3**2) * (R - r0) * psi_func(R, Z)
    
    def dpsi_dz(R, Z):
        "Derivative of psi in Z"
        return - (2/0.3**2) * (Z - z0) * psi_func(R, Z)
    
    psi2d = psi_func(r2d, z2d)
    
    def fpol_func(psi):
        "Define a simple profile for poloidal current function f = R * Bt"
        return 1. - 0.1*psi**2
    
    def fpolprime_func(psi):
        "Derivative of fpol"
        return - 0.2 * psi
    
    psi1d = np.linspace(0, 1, 65)
    fpol1d = fpol_func(psi1d)
    
    eq = tokamak.TokamakEquilibrium(r1d, z1d, psi2d, psi1d, fpol1d)

    # Check interpolation of psi, f and f' at some locations
    for r, z in [(1.2, 0.1), (1.6, -0.4), (1.8, 0.9)]:
        psi = psi_func(r,z)

        # Poloidal flux
        assert np.isclose(eq.psi(r,z), psi, atol=1e-5)

        # Poloidal current function
        assert np.isclose(eq.fpol(psi), fpol_func(psi))
        assert np.isclose(eq.fpolprime(psi), fpolprime_func(psi))

        # Poloidal magnetic field
        assert np.isclose(eq.Bp_R(r, z), -dpsi_dz(r, z)/r, atol=1e-3)
        assert np.isclose(eq.Bp_Z(r, z), dpsi_dr(r, z)/r, atol=1e-3)

        # vector Grad(psi)/|Grad(psi)|**2
        assert np.isclose(eq.f_R(r, z), dpsi_dr(r, z) / np.sqrt(dpsi_dr(r, z)**2 +
                                                                dpsi_dz(r, z)**2),
                          atol=1e-3)
        assert np.isclose(eq.f_Z(r, z), dpsi_dz(r, z) / np.sqrt(dpsi_dr(r, z)**2 +
                                                                dpsi_dz(r, z)**2),
                          atol=1e-3)
        

def test_read_geqdsk():
    # Number of mesh points
    nx = 65
    ny = 65

    # Limits of the domain
    Rmin = 1.0
    Rmax = 2.0
    Zmin = -1.0
    Zmax = 1.0

    # Centre of "plasma"
    r0 = 1.1
    z0 = 0.2
    
    # A poloidal flux function
    def psi_func(R,Z):
        return - 1.5 * np.exp(-((R - r0)**2 + (Z - z0)**2)/0.3**2)
    
    def dpsi_dr(R, Z):
        "Derivative of psi in R"
        return - (2/0.3**2) * (R - r0) * psi_func(R, Z)
    
    def dpsi_dz(R, Z):
        "Derivative of psi in Z"
        return - (2/0.3**2) * (Z - z0) * psi_func(R, Z)

    psi_boundary = psi_func(1.5, z0)
    
    def fpol_func(psi):
        "Define a simple profile for poloidal current function f = R * Bt"
        psi = np.clip(psi, None, psi_boundary)
        return 1. - 0.1*psi**2
    
    def fpolprime_func(psi):
        "Derivative of fpol"
        return np.where(psi < psi_boundary, - 0.2 * psi, 0.0)
    
    psi1d = np.linspace(psi_func(r0, z0), psi_boundary, nx)
    
    data = {"nx":nx, "ny":ny,
            "rdim":Rmax - Rmin,
            "zdim":Zmax - Zmin,
            "rleft": Rmin,
            "rcentr": r0,
            "bcentr": 1.0,
            "zmid": 0.5*(Zmax + Zmin),
            "rmagx": r0,
            "zmagx": z0,
            "simagx": psi_func(r0, z0),
            "sibdry": psi_func(1.5, z0),
            "cpasma": 1234521,
            "fpol": fpol_func(psi1d),
            "pres": np.zeros(nx),
            "qpsi": np.zeros(nx),
            "psi": psi_func(*np.meshgrid(np.linspace(Rmin, Rmax, nx),
                                         np.linspace(Zmin, Zmax, ny),
                                         indexing='ij'))}
    
    # Write to string
    output = StringIO()
    _geqdsk.write(data, output)

    # Move to the beginning of the buffer
    output.seek(0)

    # Read from string
    eq = tokamak.read_geqdsk(output)
    
    # Check interpolation of psi, f and f' at some locations
    for r, z in [(1.2, 0.1), (1.6, -0.4), (1.8, 0.9)]:
        psi = psi_func(r,z)

        # Poloidal flux
        assert np.isclose(eq.psi(r,z), psi, atol=1e-5)

        # Poloidal current function
        assert np.isclose(eq.fpol(psi), fpol_func(psi))
        assert np.isclose(eq.fpolprime(psi), fpolprime_func(psi))

        # Poloidal magnetic field
        assert np.isclose(eq.Bp_R(r, z), -dpsi_dz(r, z)/r, atol=1e-3)
        assert np.isclose(eq.Bp_Z(r, z), dpsi_dr(r, z)/r, atol=1e-3)

        # vector Grad(psi)/|Grad(psi)|**2
        assert np.isclose(eq.f_R(r, z), dpsi_dr(r, z) / np.sqrt(dpsi_dr(r, z)**2 +
                                                                dpsi_dz(r, z)**2),
                          atol=1e-3)
        assert np.isclose(eq.f_Z(r, z), dpsi_dz(r, z) / np.sqrt(dpsi_dr(r, z)**2 +
                                                                dpsi_dz(r, z)**2),
                          atol=1e-3)
