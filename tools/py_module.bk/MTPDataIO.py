#!/usr/bin/env python3

import re
class MTPDataIO:
    def __init__(self, inputObj):
        import os.path
        filename = inputObj.datafile
        if os.path.isfile(filename):
            self.filename = inputObj.datafile
        else:
            print("=============== Error ===============")
            print("  data file is not found.")
            print("=====================================")
            exit()
        self.get_mtinfo()
        self.parse_mtp()
        self.rmin = inputObj.rmin
        self.rmax = inputObj.rmax
        self.inputObj = inputObj

    def get_rmesh(self, strings):
        return int(re.sub('\)', '', re.sub('R\(', '', strings.split('=')[0])))

    def get_mtinfo(self):
        self.rmt = []
        self.rmeshdb = []
        for line in open(self.filename, 'r'):
            if re.search('R\(1\)', line) and re.search('SATOM', line):
                data_line = re.sub(',', ' ', re.sub('=\s+', '=', line.strip()))
                data_array = data_line.split()
                rmt = float(data_array[3].split('=')[1])
                rmesh = int(self.get_rmesh(data_array[2]))
                self.rmt.append(rmt)
                self.rmeshdb.append(rmesh)
        self.nkat = len(self.rmeshdb)

    def get_rj(self, iatom, j):
        #
        # convert R(ikat, i) -> real radius
        #
        # ref.) sub. struct in initc
        #
        import math
        xmin = -9.0
        jc = self.rmeshdb[iatom-1]
        dr = (self.rmt[iatom-1] - xmin)/(jc-1)
        xj = xmin + dr * (j-1)
        return math.exp(xj)

    def parse_mtp(self):
        import re
        import numpy as np
        search_line = False
        self.iter = 0
        for line in open(self.filename):
            if re.search('V\*R\s+IN\s+MT', line.upper()):
                self.iter += 1
                self.rmtp = []
                self.mtp = []
                self.r_i = []
                self.rmesh = []
                search_line = True
            else:
                if re.search('<<<', line):
                    search_line = False
                else:
                    if search_line:
                        if re.search('IKAT', line):
                            linedata = re.sub('\s*=\s*', '=', line).split()
                            ikat = int(linedata[0].split('=')[1])
                            isp = int(linedata[1].split('=')[1])
                            if len(self.rmtp) < ikat:
                                self.rmtp.append([])
                                self.mtp.append([])
                                self.r_i.append([])
                                self.rmesh.append([])
                            if len(self.rmtp[ikat-1]) < isp:
                                self.rmtp[ikat-1].append(np.array([]))
                                self.mtp[ikat-1].append(np.array([]))
                                self.r_i[ikat-1].append(np.array([]))
                                self.rmesh[ikat-1].append(np.array([]))
                        else:
                            linedata = line.strip().split()
                            r = 0.0
                            for index, value in enumerate(linedata):
                                if index % 2 == 0:
                                    ri = int(value)
                                    self.r_i[ikat-1][isp-1] =\
                                        np.append(self.r_i[ikat-1][isp-1], ri)
                                    r = self.get_rj(ikat, ri)
                                    self.rmesh[ikat-1][isp-1] =\
                                        np.append(self.rmesh[ikat-1][isp-1], r)
                                else:
                                    pot = float(re.sub('D', 'e', value))
                                    self.rmtp[ikat-1][isp-1] =\
                                        np.append(self.rmtp[ikat-1][isp-1], pot)
                                    self.mtp[ikat-1][isp-1] =\
                                        np.append(self.mtp[ikat-1][isp-1], pot / r)
        self.nsp = len(self.rmtp[0])
