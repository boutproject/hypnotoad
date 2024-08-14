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
    sigma_Bp: sign of dpsi/drho (+1 or -1) for Ip > 0 and Bt > 0
    sigma_cylind: +1 for (R, phi, Z), and -1 for (R, Z, phi)
    sigma_poloid: +1 for (rho, theta, phi), and -1 for (rho, phi, theta)
    sigma_qsafe: sign of qsafe (+1/-1) for Ip > 0 and Bt > 0
    sigma_pprime: sign of pprime (+1/-1) for Ip > 0 and Bt > 0
    theta_clockwize: True for clockwize and False for cnt-clockwize
    """

    def __init__(self, cocos_index):

        cocos_index_sign = np.sign(cocos_index)

        if cocos_index_sign < 0:
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
            self.sigma_Bp = +1.0
            self.sigma_cylind = +1.0
            self.sigma_poloid = +1.0
            self.sigma_qsafe = +1.0
            self.sigma_pprime = -1.0

        elif cocos_index_1st_digit == 2:
            self.sigma_Bp = +1.0
            self.sigma_cylind = -1.0
            self.sigma_poloid = +1.0
            self.sigma_qsafe = +1.0
            self.sigma_pprime = -1.0

        elif cocos_index_1st_digit == 3:
            self.sigma_Bp = -1.0
            self.sigma_cylind = +1.0
            self.sigma_poloid = -1.0
            self.sigma_qsafe = -1.0
            self.sigma_pprime = +1.0

        elif cocos_index_1st_digit == 4:
            self.sigma_Bp = -1.0
            self.sigma_cylind = -1.0
            self.sigma_poloid = -1.0
            self.sigma_qsafe = -1.0
            self.sigma_pprime = +1.0

        elif cocos_index_1st_digit == 5:
            self.sigma_Bp = +1.0
            self.sigma_cylind = +1.0
            self.sigma_poloid = -1.0
            self.sigma_qsafe = -1.0
            self.sigma_pprime = -1.0

        elif cocos_index_1st_digit == 6:
            self.sigma_Bp = +1.0
            self.sigma_cylind = -1.0
            self.sigma_poloid = -1.0
            self.sigma_qsafe = -1.0
            self.sigma_pprime = -1.0

        elif cocos_index_1st_digit == 7:
            self.sigma_Bp = -1.0
            self.sigma_cylind = +1.0
            self.sigma_poloid = +1.0
            self.sigma_qsafe = +1.0
            self.sigma_pprime = +1.0

        elif cocos_index_1st_digit == 8:
            self.sigma_Bp = -1.0
            self.sigma_cylind = -1.0
            self.sigma_poloid = +1.0
            self.sigma_qsafe = +1.0
            self.sigma_pprime = +1.0

        else:
            raise ValueError("Invalid COCOS input:", cocos_index)

        self.cocos_index = cocos_index
        self.theta_clockwize = self.sigma_cylind * self.sigma_poloid > 0
