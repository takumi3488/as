#!/usr/bin/env python3
#
# plot for core level
#   written by Hiroki Funashima
#

class CoreDataIO:
    def __init__(self, inputObj):
        import os.path
        datafile = inputObj.datafile
        if os.path.isfile(datafile):
            self.datafile = datafile
        else:
            print("=============== Error ===============")
            print("  data file is not found.")
            print("=====================================")
            exit()
        self.maxiter = inputObj.iter
        self.parse()

    def parse(self):
        import re
        search_line = False
        self.iter = 1
        self.l_list = []
        self.energy_list = []
        self.energy_min = 0.0
        ikat = 0
        for line in open(self.datafile, 'r'):
            if re.match('\s+=+\s+EIGENENERGY\s+OF\s+CORE', line.upper()):
                search_line = True
                ikat = 0
                self.l_list = []
                self.energy_list = []
                self.energy_min = 0.0
            else:
                if search_line:
                    if re.match('\s+[0-9]', line):
                        data = line.strip().split()
                        ilevel = int(data[0])//10
                        energy = float(data[1])
                        if self.energy_min > energy:
                            self.energy_min = energy
                        if ilevel == 10:
                            ikat += 1
                            if len(self.l_list) < ikat:
                                self.l_list.append([])
                                self.energy_list.append([])
                        self.l_list[ikat-1].append(ilevel)
                        self.energy_list[ikat-1].append(energy)
                    else:
                        self.iter += 1
                        if self.iter > self.maxiter:
                            break
                        search_line = False
        self.nkat = ikat


class AtomName:
    def __init__(self, filename):
        self.filename = filename
        self.atom_name_list = []
        self.parse()

    def parse(self):
        import re
        search_line = False
        self.nkat = 0
        for line in open(self.filename, 'r'):
            if re.match('\s+INTERSTITIAL', line.upper()):
                search_line = True
                self.nkat = 0
            else:
                if search_line:
                    if re.match('\s+ATOM', line.upper()):
                        self.nkat += 1
                        atom = line.strip().split('*')[0].split('(')[1]
                        if len(self.atom_name_list) < self.nkat:
                            self.atom_name_list.append('')
                        self.atom_name_list[self.nkat-1] = atom
                    else:
                        search_line = False
