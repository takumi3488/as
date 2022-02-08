#!/usr/bin/env python3
#
# plot program for Logarithmic Derivative for KANSAI Pakage
#
#   written by Hiroki Funashima, 2001,2002,2003,2004,2005 Osaka University
#   modified by Hiroki Funashima, 2016,Kobe University
#    to rewrite for python3 style
#
#
from py_module.LogDvIO import *
class PlotLogDv:
    def __init__(self, inputObj):
        self.inputObj = inputObj
        self.dataObj = LogDvIO(inputObj)
        self.nwindow = inputObj.nwindow

    def nsp2num(self, keyword):
        if keyword.lower() == 'up':
            return 0
        elif keyword.lower() == 'down':
            return 1
        else:
            print('==== error(nsp2num) ====')
            print(' unknown spin type')
            print('========================')
            exit()

    def lchar2lnum(self, lchar):
        if lchar.lower() == 's':
            return 0
        elif lchar.lower() == 'p':
            return 1
        elif lchar.lower() == 'd':
            return 2
        else:
            return ord(lchar)-99

    def lnum2lchar(self, lnum):
        if lnum == 0:
            return 's'
        elif lnum == 1:
            return 'p'
        elif lnum == 2:
            return 'd'
        else:
            return chr(lnum+99)

    def lparity(self, lchar):
        if self.lchar2lnum(lchar) % 2 == 0:
            return 1
        else:
            return -1

    def phiparity(self, phi0, phi1):
        ne0 = len(phi0)
        ne1 = len(phi1)
        if (phi1[ne1-1] - phi1[0]) * (phi0[ne0-1] - phi0[0]) > 0:
            return 1
        else:
            return -1

    def even_or_odd(self, lchar, phi0, phi1):
        if self.lparity(lchar) * self.phiparity(phi0, phi1) == 1:
            return 'odd'
        else:
            return 'even'

    def draw(self):
        import numpy as np
        import matplotlib.pyplot as plt
        from scipy.interpolate import interp1d

        mark = ['o', '^', 'x', 'd'] # s, p, d, f
        colr = ['red', 'green', 'blue', 'k', 'orange', 'violet']

        energy_scale =\
            self.dataObj.get_logdv(ia=1, isp=1, l=0, keyword='energy')
        if self.inputObj.emin is None:
            emin = energy_scale[0]
        else:
            emin = self.inputObj.emin

        if self.inputObj.emax is None:
            emax = energy_scale[len(energy_scale) - 1]
        else:
            emax = self.inputObj.emax

        logdrv_min = self.inputObj.logdrv_min
        logdrv_max = self.inputObj.logdrv_max
        nsp = self.dataObj.nsp

        self.atom_name_list = self.inputObj.fort1Obj.atom_name_list
        self.valence = self.inputObj.fort1Obj.valence

        #
        # add by HF on 1st June, 2017
        #
        plt.figure(figsize=(16,12))

        plt.xlim([emin, emax])
        for obj in self.inputObj.dataset:
            ia = obj.ia
            if ia > self.dataObj.nka:
                continue
            #
            # implict nsp
            #
            if obj.nsp() == 0:
                if nsp == 1:
                    obj.spin_list = ['up']
                elif nsp == 2:
                    obj.spin_list = ['up', 'down']
            elif obj.nsp() == 2:
                if nsp == 1:
                    obj.spin_list = ['up']
            for spname in obj.spin_list:
                for l_name in obj.l_list:
                    isp = self.nsp2num(spname)
                    l = self.lchar2lnum(l_name)

                    if self.nwindow:
                        #
                        # valence
                        #
                        #
                        # escape non-valence case
                        #
                        if len(self.valence[ia - 1][l]) < 1:
                            continue
                        energy2 = self.dataObj.get_logdv(
                                ia=ia, isp=isp, l=l, keyword='energy')
                        phi02 = self.dataObj.get_logdv(
                                ia=ia, isp=isp, l=l, keyword='phi0')
                        phi12 = self.dataObj.get_logdv(
                                ia=ia, isp=isp, l=l, keyword='phi1')
                        f02 = interp1d(energy2, phi02, kind='cubic')
                        f12 = interp1d(energy2, phi12, kind='cubic')

                        ne2 = len(energy2)
                        enew2 =\
                            np.linspace(energy2[0], energy2[ne2-1], 20*ne2+1)
                        parity2 = self.even_or_odd(l_name, phi02, phi12)
                        parity_title2 = ' (n:' + parity2 + ')'
                        parity_title2 = ''

                        #
                        # semicore
                        #
                        energy = self.dataObj.get_logdv2(
                                ia=ia, isp=isp, l=l, keyword='energy')
                        phi0 = self.dataObj.get_logdv2(
                                ia=ia, isp=isp, l=l, keyword='phi0')
                        phi1 = self.dataObj.get_logdv2(
                                ia=ia, isp=isp, l=l, keyword='phi1')
                        f0 = interp1d(energy, phi0, kind='cubic')
                        f1 = interp1d(energy, phi1, kind='cubic')
   
                        ne = len(energy)
                        enew = np.linspace(energy[0], energy[ne-1], 20*ne+1)
                        parity = self.even_or_odd(l_name, phi0, phi1)
                        parity_title = ' (n:' + parity + ')'
                        parity_title = ''
                    else:
                        #
                        # escape non-valence case
                        #
                        if len(self.valence[ia - 1][l]) < 1:
                            continue
                        energy = self.dataObj.get_logdv(
                                ia=ia, isp=isp, l=l, keyword='energy')
                        phi0 = self.dataObj.get_logdv(
                                ia=ia, isp=isp, l=l, keyword='phi0')
                        phi1 = self.dataObj.get_logdv(
                                ia=ia, isp=isp, l=l, keyword='phi1')
                        f0 = interp1d(energy, phi0, kind='cubic')
                        f1 = interp1d(energy, phi1, kind='cubic')

                        ne = len(energy)
                        enew = np.linspace(energy[0], energy[ne-1], 20*ne+1)
                        parity = self.even_or_odd(l_name, phi0, phi1)
                        parity_title = ' (n:' + parity + ')'
                        parity_title = ''

                    try:
                        label = self.atom_name_list[ia-1] + '(IA=' + str(ia) + ')' +\
                        ' L=' + str(l) + \
                         parity_title
                    except IndexError:
                        label = '(IA=' + str(ia) + ')' +\
                        ' L=' + str(l) + \
                         parity_title

                    if nsp == 2:
                        if isp == 1:
                            label += ' (up)'
                        elif isp == 2:
                            label += ' (down)'

                    if self.nwindow:
                        #
                        # window 1
                        #
                        plt.subplot(3, 2, 1)
                        plt.plot(enew, f0(enew), '-',color=colr[ia-1], label='')
                        if len(self.valence[ia - 1][l]) > 1:
                            plt.plot(energy, phi0, marker=mark[l], markeredgecolor=colr[ia-1],\
                                    markerfacecolor="white", label=label)
                        else:
                            plt.plot(energy, phi0, marker=mark[l], color=colr[ia-1], label=label)
                        plt.xlim([emin, emax])
                        plt.xlabel('Energy(Ry.)')
                        plt.ylabel('$\phi(r;E)$')
                        plt.title('Window1(semicore)')
                        plt.legend()
                        plt.axhline(y=0, color='k')

                        plt.subplot(3, 2, 3)
                        plt.plot(enew, f1(enew), '-',color=colr[ia-1], label='')
                        if len(self.valence[ia - 1][l]) > 1:
                            plt.plot(energy, phi1, marker=mark[l], markeredgecolor=colr[ia-1],\
                                    markerfacecolor="white", label=label)
                        else:
                            plt.plot(energy, phi1, marker=mark[l], color=colr[ia-1], label=label)
                        plt.xlim([emin, emax])
                        plt.xlabel('Energy(Ry.)')
                        plt.ylabel('$\phi^\prime(r;E)$')
                        plt.legend()
                        plt.axhline(y=0, color='k')

                        plt.subplot(3, 2, 5)
                        plt.plot(enew, f1(enew)/f0(enew), '-',color=colr[ia-1], label='')
                        if len(self.valence[ia - 1][l]) > 1:
                            plt.plot(energy, phi1/phi0, marker=mark[l], markeredgecolor=colr[ia-1],\
                                    markerfacecolor="white", label=label)
                        else:
                            plt.plot(energy, phi1/phi0, marker=mark[l], color=colr[ia-1], label=label)
                        plt.xlim([emin, emax])
                        plt.ylim([logdrv_min, logdrv_max])
                        plt.xlabel('Energy(Ry.)')
                        plt.ylabel('$\phi^\prime(r;E)/\phi(r;E)$')
                        plt.legend()
                        plt.axhline(y=0, color='k')

                        #
                        # window 2
                        #
                        plt.subplot(3, 2, 2)
                        plt.plot(enew2, f02(enew2), '-',color=colr[ia-1], label='')
                        plt.plot(energy2, phi02, mark[l],color=colr[ia-1], label=label)
                        plt.xlim([emin, emax])
                        plt.xlabel('Energy(Ry.)')
                        plt.ylabel('$\phi(r;E)$')
                        plt.title('Window2(valence)')
                        plt.legend()
                        plt.axhline(y=0, color='k')

                        plt.subplot(3, 2, 4)
                        plt.plot(enew2, f12(enew2), '-',color=colr[ia-1], label='')
                        plt.plot(energy2, phi12, mark[l],color=colr[ia-1], label=label)
                        plt.xlim([emin, emax])
                        plt.xlabel('Energy(Ry.)')
                        plt.ylabel('$\phi^\prime(r;E)$')
                        plt.legend()
                        plt.axhline(y=0, color='k')

                        plt.subplot(3, 2, 6)
                        plt.plot(enew2, f12(enew2)/f02(enew2), '-',color=colr[ia-1], label='')
                        plt.plot(energy2, phi12/phi02, mark[l],color=colr[ia-1], label=label)
                        plt.xlim([emin, emax])
                        plt.ylim([logdrv_min, logdrv_max])
                        plt.xlabel('Energy(Ry.)')
                        plt.ylabel('$\phi^\prime(r;E)/\phi(r;E)$')
                        plt.legend()
                        plt.axhline(y=0, color='k')

                    else:
                        plt.subplot(3, 1, 1)
                        plt.plot(enew, f0(enew), '-', color=colr[ia-1], label='')
                        plt.plot(energy, phi0, mark[l],color=colr[ia-1], label=label)
                        plt.xlim([emin, emax])
                        plt.xlabel('Energy(Ry.)')
                        plt.ylabel('$\phi(r;E)$')
                        plt.legend()
                        plt.axhline(y=0, color='k')

                        plt.subplot(3, 1, 2)
                        plt.plot(enew, f1(enew), '-', color=colr[ia-1], label='')
                        plt.plot(energy, phi1, mark[l],color=colr[ia-1], label=label)
                        plt.xlim([emin, emax])
                        plt.xlabel('Energy(Ry.)')
                        plt.ylabel('$\phi^\prime(r;E)$')
                        plt.legend()
                        plt.axhline(y=0, color='k')

                        plt.subplot(3, 1, 3)
                        plt.plot(enew, f1(enew)/f0(enew), '-', color=colr[ia-1], label='')
                        plt.plot(energy, phi1/phi0, mark[l],color=colr[ia-1], label=label)
                        plt.xlim([emin, emax])
                        plt.ylim([logdrv_min, logdrv_max])
                        plt.xlabel('Energy(Ry.)')
                        plt.ylabel('$\phi^\prime(r;E)/\phi(r;E)$')
                        plt.legend()
                        plt.axhline(y=0, color='k')


        plt.suptitle(self.inputObj.datafile)
        #plt.show()
        plt.savefig('levsts.png')
