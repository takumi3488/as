#!/usr/bin/env python3
#
#  parse fort.2(disperision and it's l-component)
#   written by Hiroki Funashima, Kobe Univeristy, 2017
#
#
from py_module.ParseFort1 import *
from py_module.Estate import *
import math
import re


class ParseFort2:
    def __init__(self, fort1, fort2):
        self.fort1_obj = ParseFort1(fort1)
        self.nka = self.fort1_obj.nka
        self.atom_name = self.fort1_obj.atom_name
        self.datafile = fort2
        self.search_size()
        self.init_kpoint_info()
        self.parse()
        self.nkp = self.ikp

    def parse(self):
        neig_lines = 0
        clpm_lines = self.set_clpm_lines()
        line_type = 0
        neig = None
        self.ikp = 0
        try:
            for line in open(self.datafile, 'r'):
                line_type += 1

                if re.search('[A-Z]', line.strip().upper()):

                    eeig = []
                    k_info = {}
                    clpm = []
                    line_type = 0
                    #
                    # read new data.
                    #
                    k_info = self.parse_k_info(line)
                    kx, ky, kz, ic = k_info['coordinate']
                    self.ic = ic
                    neig = k_info['neig']
                    neig_lines = self.set_neig_lines(neig)
                    for i in range(neig):
                        clpm.append([])

                else:
                    if line_type <= neig_lines:
                        for en in self.line2float(line):
                            eeig.append(en)
                        #
                        # self check for number of eigenenergy
                        #
                        if line_type == neig_lines:
                            if not len(eeig) == neig:
                                self.error_neig(neig, eeig)
                    else:
                        jline = line_type - neig_lines - 1
                        ieig = jline // clpm_lines
                        for x in self.parse_clpm_info(line):
                            clpm[ieig].append(x)
                if line_type == neig_lines + clpm_lines:
                    self.regist_info(k_info, eeig, clpm)

        except IOError:
            print("===== Error(ParseFort2) =====")
            print(" file:{} is not found.".format(self.datafile))
            exit()

    def error_neig(self, neig, eeig):
            print('=== Error(num of eeig) in ParseFort2 ===')
            print(' # of eigenenergy must be {0}'.format(neig))
            print('  but, accually {1}'.format(len(eeig)))
            exit()

    def parse_k_info(self, line):
        data = line.strip().upper().split()
        k_info = {}
        k_info['coordinate'] = [
                                int(data[1]), int(data[2]),
                                int(data[3]), int(data[4])
                               ]
        k_info['ispin'] = int(data[5])
        k_info['kname'] = data[6]
        k_info['ir'] = int(data[7])
        k_info['mwei'] = int(data[8])
        if len(data) == 10:
            temp_data = list(data[9])
            nstr = ''
            neig = ''
            for i, v in enumerate(temp_data):
                if i < len(temp_data) - 3:
                    nstr += v
                else:
                    neig += v
            k_info['nstr'] = int(nstr)
            k_info['neig'] = int(neig)
        else:
            k_info['nstr'] = int(data[9])
            k_info['neig'] = int(data[10])
        return k_info

    def parse_clpm_info(self, line):
        data = self.line2float(line)
        this_clpm = []
        if len(data) % 4 != 0:
            print("=== Error(clpm size) in ParseFort2 ===")
        for j in range(len(data) // 4):
            s = data[0 + 4 * j]
            p = data[1 + 4 * j]
            d = data[2 + 4 * j]
            f = data[3 + 4 * j]
            this_clpm.append([s, p, d, f])
        return this_clpm

    def set_neig_lines(self, neig):
        if neig % 8 == 0:
            neig_lines = neig // 8
        else:
            neig_lines = (neig // 8) + 1
        return neig_lines

    def set_clpm_lines(self):
        return (self.nka // 2) + 1

    def line2float(self, line):
        return list(map(lambda x: float(x), line.strip().split()))

    def search_size(self):
        self.nx = 0
        self.ny = 0
        self.nz = 0
        self.nspin = 0
        try:
            for line in open(self.datafile, 'r'):
                linedata = line.strip().upper()
                if re.search('[A-Z]', linedata):
                    data = line.split()
                    if self.nx < int(data[1]):
                        self.nx = int(data[1])
                    if self.ny < int(data[2]):
                        self.ny = int(data[2])
                    if self.nz < int(data[3]):
                        self.nz = int(data[3])
                    if self.nspin < int(data[5]):
                        self.nspin = int(data[5])

            maxnx = max([self.nx, self.ny, self.nz])
            self.nx = maxnx
            self.ny = self.nx
            self.nz = self.nx
        except IOError:
            print("===== Error(ParseFort2) =====")
            print(" file:{} is not found.".format(self.datafile))
            exit()

    def k_obj(self, kx, ky, kz, ic=None):
        #
        # interface function
        #   to access kpoint info interface
        #
        if ic is None:
            ic = self.ic
        else:
            if ic > self.ic:
                print("===== Error(k_obj) in ParseFort2 =====")
                print("ic must be less than equal {0}".format(self.ic))
                exit()
            else:
                if self.ic % ic == 0:
                    iic = self.ic // ic
                    kx *= iic
                    ky *= iic
                    kz *= iic
                else:
                    print("===== Error(k_obj) in ParseFort2 =====")
                    print("ic must be divisor of {}".format(self.ic))
                    exit()
        jx = self.math2spec(kx)
        jy = self.math2spec(ky)
        jz = self.math2spec(kz)
        if self.defined_k(kx, ky, kz):
            return self.kpoint_info[jx][jy][jz]
        else:
            print('====== Error(k_obj) in ParseFort2 =====')
            print('this kpoint is not defined.')
            exit()

    def math2spec(self, ix):
        #
        # function to convert from array index of mathmatical
        #                       to that of python's specification
        #
        #  -- mathmatical index --
        #    -nx <= ix <= nx
        #
        #  -- specification index --
        #    0 <= jx < 2*nx
        #
        # here, jx = ix + nx ( -nx <= ix <= nx )
        #
        if -self.nx <= ix <= self.nx:
            pass
        else:
            print('==== Error(math2spec) =====')
            print(' index is invalid.')
            print(' ix = {}'.format(ix))
            exit()
        return ix + self.nx

    def spec2math(self, jx):
        #
        # function to convert from array of python's specification
        #                       to that of mathematical index
        #
        #  -- specification index --
        #    0 <= jx < 2*nx
        #
        #  -- mathmatical index --
        #    -nx <= ix < nx
        #
        # here, ix = jx - nx ( -nx <= ix <= nx )
        #
        return jx - self.nx

    def init_kpoint_info(self):
        #
        # constract kpoint_info object
        #
        self.kpoint_info = []
        for kx in range(-self.nx, self.nx + 1):
            jx = self.math2spec(kx)
            if jx == len(self.kpoint_info):
                self.kpoint_info.append([])
            for ky in range(-self.ny, self.ny + 1):
                jy = self.math2spec(ky)
                if jy == len(self.kpoint_info[jx]):
                    self.kpoint_info[jx].append([])
                for kz in range(-self.nz, self.nz + 1):
                    jz = self.math2spec(kz)
                    if jz == len(self.kpoint_info[jx][jy]):
                        self.kpoint_info[jx][jy].append([])
                    self.kpoint_info[jx][jy][jz] = None
        #
        # check size
        #
        if len(self.kpoint_info) != 2 * self.nx + 1:
            print('===== Error(init_kpoint_info) =====')
            print('size is invalid for kx')
            exit()
        for kx in range(-self.nx, self.nx + 1):
            jx = self.math2spec(kx)
            if len(self.kpoint_info[jx]) != 2 * self.ny + 1:
                print('===== Error(init_kpoint_info) =====')
                print('size is invalid for ky')
                print('for kx = {}'.format(kx))
                print('array size must be {0}'.format(2 * self.ny + 1))
                print('but acctually size is {0}'.
                      format(len(self.kpoint_info[jx])))
                exit()
            for ky in range(-self.ny, self.ny + 1):
                jy = self.math2spec(ky)
                if len(self.kpoint_info[jx]) != 2 * self.ny + 1:
                    print('===== Error(init_kpoint_info) =====')
                    print('size is invalid for ky')
                    print('for kx = {0}, ky = {1}'.format(kx, ky))
                    print('array size must be {0}'.format(2 * self.nz + 1))
                    print('but acctually size is {0}'.
                          format(len(self.kpoint_info[jx][jy])))
                    exit()
                for kz in range(-self.nx, self.nz + 1):
                    jz = self.math2spec(kz)
                    if self.kpoint_info[jx][jy][jz] is not None:
                        print('===== Error(init_kpoint_info) =====')
                        print('size is invalid for ky')
                        print('for kx = {0}, ky = {1}, kz = {2}'.
                              format(kx, ky, kz))
                        print('kpoint_info has not been defined.')
                        exit()

    def defined_k(self, ix, iy, iz):
        jx = self.math2spec(ix)
        jy = self.math2spec(iy)
        jz = self.math2spec(iz)
        if self.kpoint_info[jx][jy][jz] is None:
            return False
        else:
            return True

    def constract_kpoint_info_obj(self, ix, iy, iz):
        jx = self.math2spec(ix)
        jy = self.math2spec(iy)
        jz = self.math2spec(iz)
        self.kpoint_info[jx][jy][jz] = Estate(self.fort1_obj,
                                              self.nspin,
                                              self.ikp)
        return self

    def regist_info(self, k_info, eeig, clpm):
        kx, ky, kz, ic = k_info['coordinate']
        if not self.defined_k(kx, ky, kz):
            self.ikp += 1
            self.constract_kpoint_info_obj(kx, ky, kz)
        jx = self.math2spec(kx)
        jy = self.math2spec(ky)
        jz = self.math2spec(kz)
        self.kpoint_info[jx][jy][jz].regist(k_info, eeig, clpm)
        return self

    def show_registered_data(self):
        print('  num of kpoint = {}'.format(self.nkp))
        print()
        for kz in range(-self.nz, self.nz+1):
            for ky in range(-self.ny, self.ny+1):
                for kx in range(-self.nx, self.nx+1):
                    if self.defined_k(kx, ky, kz):
                        self.show_k_info(kx, ky, kz)

    def show_k_info(self, kx, ky, kz):
        jx = self.math2spec(kx)
        jy = self.math2spec(ky)
        jz = self.math2spec(kz)
        self.kpoint_info[jx][jy][jz].show_info()

    def orderek(self, kx, ky, kz):
        jx = self.math2spec(kx)
        jy = self.math2spec(ky)
        jz = self.math2spec(kz)
        self.kpoint_info[jx][jy][jz].show_orderek_style()

    def show_k_info0(self, k_info, eeig, clpm):
        kx, ky, kz, ic = k_info['coordinate']
        iu = k_info['ispin']
        mr = k_info['kname']
        ir = k_info['ir']
        mwei = k_info['mwei']
        nstr = k_info['nstr']
        neig = k_info['neig']
        print('infomation about ({0}, {1}, {2} / {3})'.format(kx, ky, kz, ic))
        print('  name of kpoint: {}'.format(mr))
        print('  num of star = {}'.format(nstr))
        print('  spin state = ', end='')
        if iu == 1:
            print('up-spin')
        elif iu == 2:
            print('down-spin')
        print('  index of irreducible representation = {}'.format(ir))
        print('  degeneracy = {}'.format(mwei))
        print()
        print('  num of eigen state = {}'.format(neig))
        for i in range(neig):
            if i == 0:
                print('   E = ', end='')
            if i > 0 and i % 8 == 0:
                print('       ', end='')
            if eeig[i] >= 0:
                print(' {0:.5f} '.format(eeig[i]), end='')
            else:
                print('{0:.5f} '.format(eeig[i]), end='')

            if i % 8 == 7 or i == neig-1:
                print()
        print()
        print('  infomation for clpm:')
        for i in range(neig):
            if eeig[i] >= 0.0:
                print('       E =  {0:.5f}'.format(eeig[i]))
            else:
                print('       E = {0:.5f}'.format(eeig[i]))
            for ika in range(self.nka):
                #
                # for l-component
                #
                print('         ({0:2s}) '.format(self.atom_name[ika]), end='')
                print(' s = {0:.5f}'.format(clpm[i][ika][0]), end='')
                print(' p = {0:.5f}'.format(clpm[i][ika][1]), end='')
                print(' d = {0:.5f}'.format(clpm[i][ika][2]), end='')
                print(' f = {0:.5f}'.format(clpm[i][ika][3]))
        print()
