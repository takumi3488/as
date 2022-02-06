#!/usr/bin/env python3
#
# plot for core level
#   written by Hiroki Funashima
#

from py_module.CoreDataIO import *

class PlotCoreLevel:
    def __init__(self, inputObj):
        self.dataObj = CoreDataIO(inputObj)
        self.atom_name_list = AtomName(inputObj.datafile).atom_name_list
        self.inputObj = inputObj

    def lnum2lchar(self, lnum):
        if lnum == 0:
            return 's'
        elif lnum == 1:
            return 'p'
        elif lnum == 2:
            return 'd'
        else:
            return chr(lnum+99)

    def draw(self):
        import numpy as np
        import matplotlib.pyplot as plt

        nkat = self.dataObj.nkat
        l_list = self.dataObj.l_list
        energy_list = self.dataObj.energy_list
        #atom_name_list = self.atom_name_list
        atom_name_list = self.inputObj.fort1Obj.atom_name_list
        #
        # plotting parameter
        #
        line_length = 20.0
        char_space = 5.0
        space = 10.0
        xmax = nkat * (line_length + char_space + space)
        xx = xmax + space + space
        plt.xlim([-space, xmax + space])
        plt.ylim([-10000.0, -1.0])
        lcolor = ['r', '#ffa500', 'g', 'b']

        for ikat in range(nkat):
            ncomp = len(l_list[ikat])
            for i in range(ncomp):
                n = l_list[ikat][i]//10
                l = l_list[ikat][i] % 10
                l_name = str(n) + self.lnum2lchar(l)
                x1 = ikat * (line_length + char_space + space)
                x2 = x1 + line_length
                energy = energy_list[ikat][i]
                if i == 0:
                    plt.text(x1, -2.0, atom_name_list[ikat])
                plt.axhline(y=energy, xmin=x1/xx, xmax=x2/xx, color=lcolor[l])
                plt.text(x2-0.5*space, energy, l_name, ha='right', va='center')
        plt.yscale('symlog')
        plt.ylabel('(Ry.)')
        plt.xticks([], [])
        plt.title('Eigenenergy of core state')
        plt.show()
