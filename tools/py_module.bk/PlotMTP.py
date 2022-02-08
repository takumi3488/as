#!/usr/bin/env python3

class PlotMTP:
    def __init__(self, dataObj):
        self.dataObj = dataObj
        self.nkat = self.dataObj.nkat
        self.nsp = self.dataObj.nsp
        self.rmin = self.dataObj.rmin
        self.rmax = self.dataObj.rmax
        self.atom_name_list =\
            dataObj.inputObj.fort1Obj.atom_name_list

    def draw(self):
        import numpy as np
        import matplotlib.pyplot as plt
        from scipy.interpolate import interp1d
        #
        # auto radius range in axis
        #
        if self.rmin is None:
            rmin = 0.0
        if self.rmax is None:
            rmax = 0.0
            for ikat in range(self.nkat):
                isp = 1
                rmesh = self.dataObj.rmesh[ikat][isp-1]
                nrmesh = len(rmesh)
                if rmax < rmesh[nrmesh-1]:
                    rmax = rmesh[nrmesh-1]
            self.rmax = rmax + 5.0
        plt.xlim([self.rmin, self.rmax])
        colr = ['red', 'green', 'blue', 'k', 'orange', 'violet']
        for ikat in range(self.nkat):
            for isp in range(self.nsp):
                r = self.dataObj.rmesh[ikat][isp]
                nr = len(r)
                rpot = self.dataObj.rmtp[ikat][isp]
                pot = self.dataObj.mtp[ikat][isp]
                rf = interp1d(r, rpot, kind='slinear')
                rnew = np.linspace(r[0], r[nr-1], 20*nr+1)
                label = self.atom_name_list[ikat]
                label += ' (IA=' + str(ikat) + ')'
                if self.nsp == 2:
                    if isp == 1:
                        label += ' (up)'
                    else:
                        label += ' (down)'
                plt.plot(rnew, rf(rnew), '-',color=colr[ikat], label='')
                plt.plot(r, rpot, 'o',color=colr[ikat], label=label)
                plt.xlim([self.rmin, self.rmax])
                plt.xlabel('r[Bohr]')
                plt.ylabel('$rV(r)$')
                plt.title('Potential $r*V(r)$ in Muffin-Tin')
                plt.legend()
        plt.show()
