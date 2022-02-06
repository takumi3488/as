#!/usr/bin/env python3
#
# == orderek ==
#   original code was written by H.Harima
#   This version is implementd version written in Python
#   from scratch by Hiroki Funashima to test spintexture suite
#
import re
import math
from py_module.ParseFort1 import *
from py_module.ParseFort2 import *
from py_module.PointGroup import *


class mod_orderek:
    def __init__(self, fort1, fort2):
        self.fort1_obj = ParseFort1(fort1)
        self.fort2_obj = ParseFort2(fort1, fort2)
        self.header()

    def header(self):
        print('============== welcome to orderek(ver.2.0.0) ===============')
        print('                    original orderek was written by H.Harima')
        print('                   python version was written by H.Funashima')
        print()
        self.display_pointgroup()

    def display_pointgroup(self):
        generator_list = self.fort1_obj.generator_list
        il = self.fort1_obj.il
        pg_obj = PointGroup(il, generator_list)
        self.fort1_obj.show_info()
        print('  point group: {}'.format(pg_obj.name))
        print()

    def standard_input(self):
        my_input = True
        while my_input:
            try:
                print('  ENTER K-POINT. KX,KY,KZ,IC')
                input_word = input()
                data = list(map(lambda x: int(x),
                            re.sub(',', ' ', input_word).split()))
                if len(data) == 4:
                    my_input = False
                else:
                    print('input is invalid. try again')
            except ValueError:
                print('input is invalid. try again')
        self.data = data
        return self

    def main(self):
        fort2_obj = self.fort2_obj
        kx, ky, kz, ic = self.data
        if ic == 0:
            print('bye')
        else:
            if ic < 0:
                kx = -kx
                ky = -ky
                kz = -kz
                ic = -ic
            iic = math.gcd(ic, fort2_obj.ic)
            if ic > fort2_obj.ic:
                print('no data')
            elif ic < fort2_obj.ic:
                if fort2_obj.ic % ic > 0:
                    print('no data')
                else:
                    iic = fort2_obj.ic // ic
                    kx = kx * iic
                    ky = ky * iic
                    kz = kz * iic
                    ic = fort2_obj.ic
        if fort2_obj.defined_k(kx, ky, kz):
            fort2_obj.orderek(kx, ky, kz)
        else:
            print('no data')
