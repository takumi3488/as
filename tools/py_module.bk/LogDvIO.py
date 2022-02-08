#!/usr/bin/env python3
#
# plot program for Logarithmic Derivative for KANSAI Pakage
#
#   written by Hiroki Funashima, 2001,2002,2003,2004,2005 Osaka University
#   modified by Hiroki Funashima, 2016,Kobe University
#    to rewrite for python3 style
#
#
class LogDvIO:
    def __init__(self, inputObj):
        import os.path
        import numpy as np
        if os.path.isfile(inputObj.datafile):
            self.filename = inputObj.datafile
        else:
            print("=============== Error ===============")
            print("  data file is not found.")
            print("=====================================")
            exit()

        self.maxiter = inputObj.iter
        self.nwindow = inputObj.nwindow
        if self.nwindow:
            print("multi-window mode is on")
            self.maxiter *= 2
        self.iter = 0
        #
        # nka: # of kind of atom
        # nsp: # of kind of spin
        # maxl: muximum of L
        #
        #
        self.nka = 0
        self.nsp = 0
        self.maxl = 0
        self.parse_kansai_format()

    def parse_kansai_format(self):
        import re
        import numpy as np
        self.lineno = 0
        #
        # array for self.logdevdata
        #  logdv[IA-1, ISP-1, L, attr] = [....] or value
        #   attr = ind0( # of node for phi0 )
        #        = ind1( # of node for phi1 )
        #        = energy( energy range)
        #        = phi0
        #        = phi1
        #
        self.logdv = []
        if self.nwindow:
            self.logdv2 = []
        search_line = False
        for line in open(self.filename, 'r'):
            if re.match("\s+IND0", line.upper()):
                data = re.sub("=\s+", "=", line.rstrip()).upper().split()
                self.data_dict = {}
                self.lineno = 0
                for comp in data:
                    key, value = comp.split('=')
                    self.data_dict[key.upper()] = int(value)
                ia = self.data_dict['IA']
                isp = self.data_dict['ISP']
                l = self.data_dict['L']

                #
                # check max value for nka, isp and l
                #
                if ia > self.nka:
                    self.nka = ia
                if isp > self.nsp:
                    self.nsp = isp
                if l > self.maxl:
                    self.maxl = l

                ind0 = self.data_dict['IND0']
                ind1 = self.data_dict['IND1']

                if self.nwindow:
                    if len(self.logdv) < ia:
                        self.logdv.append([])
                        self.logdv2.append([])
                else:
                    if len(self.logdv) < ia:
                        self.logdv.append([])

                if self.nwindow:
                    if len(self.logdv[ia-1]) < isp:
                        self.logdv[ia-1].append([])
                        self.logdv2[ia-1].append([])
                else:
                    if len(self.logdv[ia-1]) < isp:
                        self.logdv[ia-1].append([])

                if self.nwindow:
                    if len(self.logdv[ia-1][isp-1]) < l+1:
                        self.logdv[ia-1][isp-1].append({})
                        self.logdv2[ia-1][isp-1].append({})
                else:
                    if len(self.logdv[ia-1][isp-1]) < l+1:
                        self.logdv[ia-1][isp-1].append({})

                if self.nwindow:
                    if self.iter % 2 == 0:
                        self.logdv[ia-1][isp-1][l]['ind0'] = ind0
                        self.logdv[ia-1][isp-1][l]['ind1'] = ind1
                    else:
                        self.logdv2[ia-1][isp-1][l]['ind0'] = ind0
                        self.logdv2[ia-1][isp-1][l]['ind1'] = ind1
                else:
                    self.logdv[ia-1][isp-1][l]['ind0'] = ind0
                    self.logdv[ia-1][isp-1][l]['ind1'] = ind1
                if ia == 1 and isp == 1 and l == 0:
                    self.iter += 1
                if self.nwindow:
                    if self.iter <= self.maxiter+1:
                        search_line = True
                    else:
                        self.iter -= 2
                        break
                else:
                    if self.iter <= self.maxiter:
                        search_line = True
                    else:
                        self.iter -= 1
                        break

            else:
                if search_line:
                    #-----
                    # add by HF on 26 May 2017
                    #
                    if re.match("^\s+A\(KL\)\=", line.upper()):
                        continue
                    #-----
                    else:
                        self.lineno += 1
                        data = []
                        for x in line.rstrip().split():
                            data.append(float(x))
                        if self.lineno == 1:
                            keyword = 'energy'
                        elif self.lineno == 2:
                            keyword = 'phi0'
                        elif self.lineno == 3:
                            keyword = 'phi1'
                            search_line = False
                        if self.nwindow:
                            if self.iter % 2 == 0:
                                self.logdv[ia-1][isp-1][l][keyword] =\
                                        np.array(data)
                            else:
                                self.logdv2[ia-1][isp-1][l][keyword] =\
                                        np.array(data)
                        else:
                            self.logdv[ia-1][isp-1][l][keyword] = np.array(data)

    def get_logdv(self, ia, isp, l, keyword):
        return self.logdv[ia-1][isp-1][l][keyword.lower()]

    def get_logdv2(self, ia, isp, l, keyword):
        return self.logdv2[ia-1][isp-1][l][keyword.lower()]
