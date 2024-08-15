#
# COCOS convention class
#
# Coordinate conventions based on
#     COCOS paper:
#       O. Sauter and S. Yu. Medvevdev, "Tokamak Coordinate Conventions: COCOS"
#       Comput. Phys. Commun. 184 (2013) 293
#
#      and cocos_module.f90 in CHEASE code: https://gitlab.epfl.ch/spc/chease
#
# written by Haruki SETO seto.haruki@qst.go.jp
#

import numpy as np


class Cocos:
    """
    cocos_index: +1/+11, +2/+12, +3/+13, +4/+14, +5/+15, +6/+16, +7/+17, +8/+18

        exp_Bp: 0 for poloidal flux devided by 2pi, and +1 for full poloidal flux
        sign_Bp: sign of dpsi/drho (+1 or -1) for Ip > 0 and Bt > 0
        sign_cylind: +1 for (R, phi, Z), and -1 for (R, Z, phi)
        sign_poloid: +1 for (rho, theta, phi), and -1 for (rho, phi, theta)
        sign_qsafe: sign of qsafe (+1/-1) for Ip > 0 and Bt > 0
        sign_pprime: sign of pprime (+1/-1) for Ip > 0 and Bt > 0
        theta_clockwize: True for clockwize and False for cnt-clockwize
    """

    def __init__(self, cocos_index):

        sign_cocos_index = np.sign(cocos_index)

        if sign_cocos_index < 0:
            raise ValueError("negative inputs are not supported yet:", cocos_index)

        cocos_index_2nd_digit = cocos_index // 10

        if cocos_index_2nd_digit == 0:
            self.exp_Bp = 0.0

        elif cocos_index_2nd_digit == 1:
            self.exp_Bp = 1.0

        else:
            raise ValueError("Invalid COCOS input:", cocos_index)

        cocos_index_1st_digit = cocos_index % 10

        if cocos_index_1st_digit == 1:
            self.sign_Bp = +1.0
            self.sign_cylind = +1.0
            self.sign_poloid = +1.0
            self.sign_qsafe = +1.0
            self.sign_pprime = -1.0

        elif cocos_index_1st_digit == 2:
            self.sign_Bp = +1.0
            self.sign_cylind = -1.0
            self.sign_poloid = +1.0
            self.sign_qsafe = +1.0
            self.sign_pprime = -1.0

        elif cocos_index_1st_digit == 3:
            self.sign_Bp = -1.0
            self.sign_cylind = +1.0
            self.sign_poloid = -1.0
            self.sign_qsafe = -1.0
            self.sign_pprime = +1.0

        elif cocos_index_1st_digit == 4:
            self.sign_Bp = -1.0
            self.sign_cylind = -1.0
            self.sign_poloid = -1.0
            self.sign_qsafe = -1.0
            self.sign_pprime = +1.0

        elif cocos_index_1st_digit == 5:
            self.sign_Bp = +1.0
            self.sign_cylind = +1.0
            self.sign_poloid = -1.0
            self.sign_qsafe = -1.0
            self.sign_pprime = -1.0

        elif cocos_index_1st_digit == 6:
            self.sign_Bp = +1.0
            self.sign_cylind = -1.0
            self.sign_poloid = -1.0
            self.sign_qsafe = -1.0
            self.sign_pprime = -1.0

        elif cocos_index_1st_digit == 7:
            self.sign_Bp = -1.0
            self.sign_cylind = +1.0
            self.sign_poloid = +1.0
            self.sign_qsafe = +1.0
            self.sign_pprime = +1.0

        elif cocos_index_1st_digit == 8:
            self.sign_Bp = -1.0
            self.sign_cylind = -1.0
            self.sign_poloid = +1.0
            self.sign_qsafe = +1.0
            self.sign_pprime = +1.0

        else:
            raise ValueError("Invalid COCOS input:", cocos_index)

        self.theta_clockwize = self.sign_cylind * self.sign_poloid > 0
        self.cocos_index = cocos_index


class CocosConversion(Cocos):

    def __init__(
        self, cocos_index, geqdsk_dict, positive_definite_qsafe=False, show_signs=False
    ):

        super().__init__(cocos_index)

        self.check_consistency(
            cocos_index, geqdsk_dict, positive_definite_qsafe, show_signs
        )

    def check_consistency(
        self, cocos_index, geqdsk_dict, positive_definite_qsafe, show_signs
    ):
        """
        check consistency between input GEQDSK and COCOS convension based on
        Section 5 in COCOS paper.

        geqdsk_dict - dictionary of GEQDSK read by _geqdsk.py

        positive_definite_qsafe - always use |q| rather than q (e.g. JET)
        """

        xpos = geqdsk_dict["nx"] // 2  # radial index for checking sign of 1D quantities

        sign_Ip_geqdsk = np.sign(geqdsk_dict["cpasma"])
        sign_B0_geqdsk = np.sign(geqdsk_dict["fpol"][xpos])
        sign_dpsi_geqdsk = np.sign(geqdsk_dict["sibdry"] - geqdsk_dict["simagx"])
        sign_pprime_geqdsk = np.sign(geqdsk_dict["pprime"][xpos])
        sign_qsafe_geqdsk = np.sign(geqdsk_dict["qpsi"][xpos])

        if show_signs:
            print("sign_Ip_geqdsk %d" % sign_Ip_geqdsk)
            print("sign_B0_geqdsk %d" % sign_B0_geqdsk)
            print("sign_dpsi_geqdsk %d" % sign_dpsi_geqdsk)
            print("sign_pprime_geqdsk %d" % sign_pprime_geqdsk)
            print("sign_qsafe_geqdsk %d" % sign_qsafe_geqdsk)

        if sign_dpsi_geqdsk == np.sign(sign_Ip_geqdsk * self.sign_Bp):
            self.sign_dpsi = sign_dpsi_geqdsk

        else:
            raise ValueError("geqdsk is not consistent with COCOS (dpsi)")

        if sign_pprime_geqdsk == np.sign(-1 * sign_Ip_geqdsk * self.sign_Bp):
            self.sign_pprime = sign_pprime_geqdsk

        else:
            raise ValueError("geqdsk is not consistent with COCOS (pprime)")

        sign_qsafe_cocos = np.sign(sign_Ip_geqdsk * sign_B0_geqdsk * self.sign_poloid)

        if sign_qsafe_geqdsk == sign_qsafe_cocos:
            self.sign_qsafe = sign_qsafe_cocos

        elif positive_definite_qsafe and sign_qsafe_cocos < 0:
            print("geqdsk has positive definite qsafe but qsafe should be negative")
            self.sign_qsafe = sign_qsafe_cocos
            geqdsk_dict["qpsi"] *= -1

        else:
            raise ValueError("geqdsk is not consistent with COCOS input (qsafe)")

        self.sign_B0 = sign_B0_geqdsk
        self.sign_Ip = sign_Ip_geqdsk
        self.geqdsk_dict = geqdsk_dict

    def convert_geqdsk_to_cocos_out(
        self, cocos_index_out=1, sign_B0_out=None, sign_Ip_out=None
    ):
        """
        convert geqdsk file to cocos_index_out convension.

        Based on Appendix C in COCOS paper but, l_d = 1[m], l_b = 1[T], exp_mu = 1
        are assumed for both input and output.
        """

        twopi = np.pi * 2.0

        cocos_out = Cocos(cocos_index_out)

        sign_cylind_eff = np.sign(cocos_out.sign_cylind * self.sign_cylind)
        sign_poloid_eff = np.sign(cocos_out.sign_poloid * self.sign_poloid)
        sign_Bp_eff = np.sign(cocos_out.sign_Bp * self.sign_Bp)
        exp_Bp_eff = cocos_out.exp_Bp - self.exp_Bp

        if sign_B0_out is not None:
            sign_B0_eff = np.sign(sign_B0_out * self.sign_B0)
        else:
            sign_B0_eff = sign_cylind_eff

        if sign_Ip_out is not None:
            sign_Ip_eff = np.sign(sign_Ip_out * self.sign_Ip)
        else:
            sign_Ip_eff = sign_cylind_eff

        # convert self.geqdsk_dict to COCOS_out convension rule
        self.geqdsk_dict["cpasma"] *= sign_Ip_eff
        self.geqdsk_dict["bcentr"] *= sign_B0_eff

        factor_psi = sign_Ip_eff * sign_Bp_eff * np.power(twopi, exp_Bp_eff)
        self.geqdsk_dict["simagx"] *= factor_psi
        self.geqdsk_dict["sibdry"] *= factor_psi
        self.geqdsk_dict["psi"] *= factor_psi

        self.geqdsk_dict["fpol"] *= sign_B0_eff

        factor_dpsi = sign_Ip_eff * sign_Bp_eff / np.power(twopi, exp_Bp_eff)
        self.geqdsk_dict["pprime"] *= factor_dpsi
        self.geqdsk_dict["ffprime"] *= factor_dpsi

        self.geqdsk_dict["qpsi"] *= sign_Ip_eff * sign_B0_eff * sign_poloid_eff

        return self.geqdsk_dict
