# coding:utf-8
from osgeo import ogr
from scipy.special.cython_special import gdtrix
import numpy

import setting
import CalHydrological.Common as HC

"""
Calculate the design rainstorm for the basin
"""


def DR(P, H, Cv):
    """
    Computational design of Rainstorms
    param H: Average rainfall
    :param P: Frequency
    param Cv: Coefficient of variation
    :return: Design a rainstorm
    """
    Cs = 3.5 * Cv
    a = 4.0 / (Cs * Cs)
    b = 2.0 / (H * Cs * Cv)
    a_0 = H - 2 * H * Cv / Cs
    tp = gdtrix(1, a, 1 - P)
    Kp = 1 + Cv * (Cs / 2 * tp - 2 / Cs)
    Hp = Kp * H
    return Hp


def RainPow(H1, H2):
    """
    Calculate the rainstorm intensity and rainstorm attenuation index
    param H1: Such as a 6-hour design rainstorm
    param H2: Such as 24-hour designed rainstorm
    return:Sp and n
    """
    n = 1 + 1.661 * (numpy.log10(H1) - numpy.log10(H2))
    Sp = H2 * pow(24, n - 1)
    return Sp, n


def main(inP, inFc):
    """
    param inP: Frequency
    param inFc: Basin surface

    """
    print("================= the current frequency for", inP, " ======================================")
    print("-----------------the calculation design rainstorm began to perform-----------------")
    print(" Obtain data for 6 and 24 hours...")
    listdata6 = []  # 6-hour data
    listdata24 = []  # 24-hour data
    layer_unit = inFc.GetLayer()
    for i in range(layer_unit.GetFeatureCount()):
        fc_unit = layer_unit.GetFeature(i)
        listdata6.append([fc_unit.GetField(setting.unit_fields['H_6']), fc_unit.GetField(setting.unit_fields['Cv_6'])])
        listdata24.append([fc_unit.GetField(setting.unit_fields['H_24']), fc_unit.GetField(setting.unit_fields['Cv_24'])])

    print(" Calculate and save the design rainstorm...")
    # 添加字段
    HC.CreateNewField(layer_unit, setting.unit_fields['Sp'], ogr.OFTReal)
    HC.CreateNewField(layer_unit, setting.unit_fields['n'], ogr.OFTReal)

    for i in range(layer_unit.GetFeatureCount()):
        fc_unit = layer_unit.GetFeature(i)
        H6 = DR(inP, listdata6[i][0], listdata6[i][1])
        H24 = DR(inP, listdata24[i][0], listdata24[i][1])
        temp = RainPow(H6, H24)
        HC.UpdateField(layer_unit, fc_unit, setting.unit_fields['Sp'], temp[0])
        HC.UpdateField(layer_unit, fc_unit, setting.unit_fields['n'], temp[1])
        # print("h6 h24",H6,H24)
        # print("sp,n",temp[0],temp[1])
    return temp[0],temp[1]  # Sp,n
