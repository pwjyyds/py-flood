# coding:utf-8
import math
import copy

from osgeo import ogr

import setting
import CalHydrological.Common as HC

"""
Calculate the peak flood water level using Manning's formula
    Traverse the cross-sectional points
    Calculate the cross-sectional area and wet perimeter
    Substitute the formula to calculate Q
    If the difference between Q and the Qm of the basin is too large, increase the water level by 0.1
"""


def Mnll(n, S, L, I):
    """
    Calculate the flow rate according to Manning's formula
    :param n: Roughness
    param S: Cross-sectional area for water passage
    param L: Wet perimeter
    :param I: Descending
    :return: Flow
    """
    # print(n,S,L,I)
    if n == 0:
        n = 0.025
    Q = pow(n, -1) * pow(S, 5 / 3.0) * pow(L, -2 / 3.0) * pow(I, 1 / 2.0)
    return Q


def zhouchang(points, inH, length):
    """
    Calculate the circumference of the wet week
    :param points: Coordinate data
    return: Wet circumference L
    """
    sum = 0

    for i in range(len(points) - 1):
        if inH - points[i][2] >= 0:
            Xc = setting.raster_pixel_width / length
            Yc = abs(points[i][2] - points[i + 1][2])
            sum += math.sqrt(Xc * Xc + Yc * Yc)
    if sum <= 0:
        return setting.raster_pixel_width / length
    else:
        return abs(sum)


def square(points, inH, length):
    """
    Calculate the irregular area
    param length: Unit length
    param inH: Current water level
    :param points: Coordinate data
    :return: Area
    """
    sum0 = 0.0
    # for i in range(len(p) - 1):
    #     sum0 += (p[i][0] * p[i + 1][1] - p[i + 1][0] * p[i][1])
    # sum = (abs(sum0 + (p[len(p) - 1][0] * p[0][1]) - (p[0][0] * p[len(p) - 1][1]))) / 2
    for point in points:
        # print(setting.raster_pixel_width,inH-i[1])
        if inH - point[2] >= 0:
            sum0 += (setting.raster_pixel_width / length) * (inH - point[2])

    return sum0


def main(inPoint, inLine, field_z, field_h, fields_Qm, fields_L, fields_S):
    """
    :param fields_Qm:
    :param inPoint: Breakpoint of the cross-section line
    :param inLine: Cross-section line
    :param field_z: Field names for different frequencies
    """
    print("-----------------Qm-----------------")

    print(" Get key attributes...")
    layer_dmx = inLine.GetLayer()
    HC.CreateNewField(layer_dmx, field_z, ogr.OFTReal)
    HC.CreateNewField(layer_dmx, field_h, ogr.OFTReal)
    HC.CreateNewField(layer_dmx, fields_L, ogr.OFTReal)
    HC.CreateNewField(layer_dmx, fields_S, ogr.OFTReal)
    dmxPointDs = ogr.Open(inPoint, 1)
    layer_point = dmxPointDs.GetLayer()
    tempUnitDs = ogr.Open(setting.temp_units, 1)
    layer_tempUnit = tempUnitDs.GetLayer()
    count_lyrPoints = layer_point.GetFeatureCount()
    count_layer_dmx = layer_dmx.GetFeatureCount()
    thisDmxPoint_dmxID = layer_point.GetFeature(0).GetField(setting.dmxPoints_field['DmxID'])  # 当前断面点属于的断面线ID
    sameDmxPoints = []  # Points on the same cross-section
    h_result = 0

    for i in range(count_lyrPoints + 1):  # Traverse the cross-sectional points
        if i < count_lyrPoints:
            fc_point = layer_point.GetFeature(i)
            dmxID = fc_point.GetField(setting.dmxPoints_field['DmxID'])
        else:
            fc_point = layer_point.GetFeature(count_lyrPoints - 1)
            dmxID = count_lyrPoints + 1

        if thisDmxPoint_dmxID == dmxID:
            # print("ififififififif", thisDmxPoint_dmxID, dmxID)
            # Points on the same cross-sectional line
            tempList = []
            tempList.append(fc_point.GetField(setting.dmxPoints_field['ObjectID']))

            tempList.append(dmxID)
            pointDemValue = fc_point.GetField(setting.dmxPoints_field['DemValue'])
            if pointDemValue is None or pointDemValue <= 0:
                tempList.append(500)
            else:
                tempList.append(pointDemValue)
            tempList.append(fc_point.GetField(setting.dmxPoints_field['n0']))
            tempList.append(fc_point.GetField(setting.dmxPoints_field['J']))

            filter = "basin=" + str(dmxID + 1)
            layer_tempUnit.SetAttributeFilter(filter)  # Set the property filter
            unitQm = None
            for fc_ in layer_tempUnit:
                unitQm = fc_.GetField(fields_Qm)
                break
            if unitQm is None:
                tempList.append(100)
            else:
                tempList.append(unitQm)
            sameDmxPoints.append(tempList)

        else:
            # print('thisDmxPoint_dmxID != dmxID', thisDmxPoint_dmxID, dmxID)
            # Calculate the points of the previous cross-section
            if dmxID == 0:
                aaaaaaaa = 1
            p_h = []
            for p in sameDmxPoints:
                p_h.append(p[2])
            h_0 = (sum(p_h) - max(p_h) - min(p_h)) / len(p_h)

            sameDmxPoints_more = []
            for sdp in sameDmxPoints:
                sameDmxPoints_more.append(sdp)
            n = 0
            num_insert = 30  # Number of insertions
            for i_sdp in range(len(sameDmxPoints_more) - 1):
                for j_ni in range(num_insert):
                    eveValue = sameDmxPoints_more[n + j_ni + 1][2] + (
                            sameDmxPoints_more[n + j_ni + 1][2] - sameDmxPoints_more[n + j_ni][2]) / num_insert
                    insertValue = [sameDmxPoints[i_sdp][0],
                                   sameDmxPoints[i_sdp][1],
                                   eveValue,
                                   sameDmxPoints[i_sdp][3],
                                   sameDmxPoints[i_sdp][4],
                                   sameDmxPoints[i_sdp][5],
                                   ]  # The value to be inserted
                    sameDmxPoints_more.insert(n + j_ni + 1, insertValue)
                n += num_insert + 1
            count_step = 1
            step = 0.35
            index = 1
            max_count = 85
            # Obtain the peak flood water level at each cross-sectional point
            beCalPoints = []
            mid = int(len(sameDmxPoints_more) / 2)
            if len(sameDmxPoints_more) % 2 == 0:
                mid -= 1
                beCalPoints = [sameDmxPoints_more[mid], sameDmxPoints_more[mid + 1]]

            else:
                beCalPoints = [sameDmxPoints_more[mid]]

            biggerTime = 0

            temp_history = []  # Historical Iteration data
            S_result = 0
            L_result = 0

            while True:  # Expand from the middle to both sides
                S = square(beCalPoints, h_0,
                           len(sameDmxPoints_more))  # Calculate the cross-sectional water passage area
                L = zhouchang(beCalPoints, h_0, len(sameDmxPoints_more))  # Calculate the wet perimeter
                Q = Mnll(n=beCalPoints[0][3],
                         S=S, L=L,
                         I=beCalPoints[0][4])  # Calculate the peak flood flow of Manning
                Qm = beCalPoints[0][5]
                abs_error = abs(Q - Qm) / Qm
                if abs_error <= 0.01:
                    S_result = S
                    L_result = L
                    break
                else:
                    h_0 += step
                    count_step += 1
                    if len(sameDmxPoints_more) % 2 == 0:
                        beCalPoints = sameDmxPoints_more[mid - index:mid + index + 2]
                    else:
                        beCalPoints = sameDmxPoints_more[mid - index:mid + index + 1]
                    index += 1
                    temp_history.append(abs_error)

                if index >= mid + 1 or count_step >= max_count:
                    S_result = S
                    L_result = L
                    break
            h_result = h_0
            # Assign a value to the section line
            print(thisDmxPoint_dmxID)
            fc_line = layer_dmx.GetFeature(thisDmxPoint_dmxID)
            h_z = step * count_step
            if count_step == max_count:
                h_z = step * (temp_history.index(min(temp_history)) + 1)

            HC.UpdateField(layer_dmx, fc_line, field_z, h_z)
            HC.UpdateField(layer_dmx, fc_line, fields_S, S_result)
            HC.UpdateField(layer_dmx, fc_line, fields_L, L_result)
            thisDmxPoint_dmxID = fc_point.GetField(setting.dmxPoints_field['DmxID'])
            sameDmxPoints = []

    inLine = None
    dmxPointDs = None
    print("-----------------The operation of Calculating peak flood water level has been successful-----------------")
