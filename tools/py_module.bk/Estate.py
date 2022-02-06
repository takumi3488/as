#!/usr/bin/env python


class Estate:
    def __init__(self, fort1_obj, nspin, ikp):
        self.ikp = ikp
        self.fort1_obj = fort1_obj
        self.nka = fort1_obj.nka
        self.nspin = nspin
        if nspin == 1 or nspin == 2:
            pass
        else:
            print('=== Error(Estate) ===')
            print('nspin is invalid.')
            print('nspin = {}'.format(nspin))
            exit()

        #
        # coordinate = [kx, ky, kz, ic]
        #
        max_ir = 12

        init_array = []
        for ir in range(max_ir):
            init_array.append(None)

        init_array2 = [init_array, init_array]

        #
        # initalize
        #
        self.ir = []
        self.nstar = []
        self.degeneracy = []
        self.neig = [[], []]
        self.eigen_energy = [[], []]
        self.clpm = [[], []]
        for ir in range(max_ir):
            self.ir.append(None)
            self.nstar.append(None)
            self.degeneracy.append(None)
            for isp in range(2):
                self.neig[isp].append(None)
                self.eigen_energy[isp].append(None)
                self.clpm[isp].append(None)

    def regist(self, k_info, eeig, clpm):
        max_ir = 12
        kx, ky, kz, ic = k_info['coordinate']
        iu = k_info['ispin']
        mr = k_info['kname']
        ir = k_info['ir']
        mwei = k_info['mwei']
        nstr = k_info['nstr']
        neig = k_info['neig']
        if 1 <= ir <= max_ir:
            pass
        else:
            print('=== Error(Estate) ===')
            print('ir is invalid.')
            print('ir = {}'.format(ir))
            exit()

        if 0 <= iu <= 1:
            pass
        else:
            print('=== Error(Estate) ===')
            print('ispin is invalid.')
            print('ispin = {}'.format(iu))
            exit()
        #
        # general info
        #
        self.coordinate = [kx, ky, kz, ic]
        self.name = mr
        self.nstar = nstr
        #
        # register infomation which depends on ispin and ir
        #
        self.ir[ir - 1] = True
        self.degeneracy[ir - 1] = mwei
        self.neig[iu - 1][ir - 1] = neig
        #
        # for 1st window
        #
        if self.eigen_energy[iu - 1][ir - 1] is None:
            self.eigen_energy[iu - 1][ir - 1] = eeig
            self.clpm[iu - 1][ir - 1] = clpm
        #
        # 2nd window
        #
        else:
            if eeig[0] > self.eigen_energy[iu - 1][ir - 1][-1]:
                for eig in eeig:
                    self.eigen_energy[iu - 1][ir - 1].append(eig)
                for xclpm in clpm:
                    self.clpm[iu - 1][ir - 1].append(xclpm)
            else:
                for eig in reversed(eeig):
                    self.eigen_energy[iu - 1][ir - 1].insert(0, eig)
                for xclpm in reversed(clpm):
                    self.clpm[iu - 1][ir - 1].insert(0, xclpm)
        return self

    def list_ir(self):
        ir_list = []
        for i, x in enumerate(self.ir):
            if x:
                ir_list.append(i+1)
        return ir_list

    def show_info(self):
        kx, ky, kz, ic = self.coordinate
        mr = self.name
        nstr = self.nstar
        print('infomation about ({0}, {1}, {2} / {3})'.format(kx, ky, kz, ic))
        print('  name of kpoint: {}'.format(mr))
        print('  num of star = {}'.format(nstr))
        print()
        if self.nspin == 1:
            print('  magnetic state: paramagnetic state')
        elif self.nspin == 2:
            print('  magnetic state: magnetic-ordered state')
        print()
        print('  eigen energies')
        self.show_energy_state(self.sort_energy())
        print('  irreducible representation = ', end='')
        ir_list = self.list_ir()
        if len(ir_list) == 1:
            print(ir_list[0])
        else:
            for i, ir in enumerate(ir_list):
                if i == len(ir_list) - 1:
                    print('{0}'.format(ir))
                elif i == len(ir_list) - 2:
                    print('{0} and '.format(ir), end='')
                else:
                    print('{0}, '.format(ir), end='')

        print()
        print('  infomation for each irreducible representation:')
        for ir in ir_list:
            for ispin in range(1, self.nspin+1):
                self.show_info_ispin_ir(ir, ispin)

    def show_info_ispin_ir(self, ir, ispin):
        atom_name = self.fort1_obj.atom_name
        print('  ** index of irreducible representation = {}'.format(ir))
        if self.nspin == 2:
            if ispin == 1:
                print('  spin state: up-spin')
            elif ispin == 2:
                print('  spin state: down-spin')
            else:
                print('==== Error(Estate) ====')
                print('ispin is invalid.')
                print(' ispin = {}'.format(ispin))
                exit()
        print('  degeneracy = {}'.
              format(self.degeneracy[ir - 1]))
        print('  num of eigen state = {}'.
              format(self.neig[ispin - 1][ir - 1]))
        neig = self.neig[ispin - 1][ir - 1]
        eeig = self.eigen_energy[ispin - 1][ir - 1]
        clpm = self.clpm[ispin - 1][ir - 1]

        for i in range(neig):
            if eeig[i] >= 0.0:
                print('       E =  {0:.5f}'.format(eeig[i]))
            else:
                print('       E = {0:.5f}'.format(eeig[i]))
            for ika in range(self.nka):
                #
                # for l-component
                #
                print('         ({0:2s}) '.format(atom_name[ika]), end='')
                print(' s = {0:.5f}'.format(clpm[i][ika][0]), end='')
                print(' p = {0:.5f}'.format(clpm[i][ika][1]), end='')
                print(' d = {0:.5f}'.format(clpm[i][ika][2]), end='')
                print(' f = {0:.5f}'.format(clpm[i][ika][3]))
        print()

    def show_energy_state(self, sorted_array):
        for ispin in range(self.nspin):
            ne = len(sorted_array[ispin])
            print('  ispin = {}'.format(ispin + 1))
            print('  num of total eigen state = {}'.format(ne))
            print('    E = ', end='')
            for i in range(ne):
                ir = sorted_array[ispin][i][0]
                e_index = sorted_array[ispin][i][1]
                energy = self.eigen_energy[ispin][ir][e_index]
                if i % 6 == 0 and i > 1:
                    print('        ', end='')
                if energy > 0.0:
                    print(' {0:.5f} '.format(energy), end='')
                else:
                    print('{0:.5f} '.format(energy), end='')
                if i % 6 == 5:
                    print()
                elif i == ne - 1 and i % 6 < 5:
                    print()
            print()

    def display_orderek_bar(self, symbol):
        for i in range(31):
            print('{}'.format(symbol), end='')
        for ika in range(self.nka):
            for i in range(13):
                print('{}'.format(symbol), end='')
        for i in range(6):
            print('{}'.format(symbol), end='')
        print()
        return self

    def show_orderek_style(self):
        kpname = self.name
        self.display_orderek_bar('=')
        print('nband  name IR nDeg  Energy(Ry)', end='')
        for ika in range(self.nka):
            print('   s  p  d  f', end='')
        print('   out')
        self.display_orderek_bar('-')
        for ispin in range(self.nspin):
            sorted_array = self.sort_energy()
            ne = len(sorted_array[ispin])
            total_band = 0
            for info in sorted_array[ispin]:
                total_band += self.degeneracy[info[0]]

            for info in reversed(sorted_array[ispin]):
                ir, e_index = info
                ndeg = self.degeneracy[ir]
                nband = total_band + 1 - ndeg
                total_band -= ndeg
                iir = ir + 1
                energy = self.eigen_energy[ispin][ir][e_index]
                clpm = self.clpm[ispin][ir][e_index]
                print('{0:5d}'.format(nband), end='')
                print('    {0:2s}'.format(kpname), end='')
                print(' {0:2d}'.format(ir + 1), end='')
                print('   {0:2d}    '.format(ndeg), end='')
                if energy >= 0.0:
                    print(' {0:.5f}'.format(energy), end='')
                else:
                    print('{0:.5f}'.format(energy), end='')
                out = 1.0
                for ika in range(self.nka):
                    print(' ', end='')
                    for i in range(4):
                        out -= clpm[ika][i]
                        print(' {0:2d}'.
                              format(int(clpm[ika][i] * 100)), end='')

                print('    {0:2d}'.format(int(out * 100)))
        self.display_orderek_bar('=')

    def check_sort(self, need_sort, control):
        local_check = []
        n = len(need_sort)
        for j in range(n):
            if control[j]:
                local_check.append(need_sort[j])
        return any(local_check)

    def sort_multi_array(self, ary, control):
        flag = []
        need_sort = []
        sorted_array = []

        #
        # initialize
        #
        n = len(ary)
        for j in range(n):
            need_sort.append(True)
            flag.append(0)

        while self.check_sort(need_sort, control):
            record_value = None
            for j in range(n):
                if not need_sort[j] or not control[j]:
                    continue
                else:
                    if record_value is None:
                        record_value = ary[j][flag[j]]
                        jj = j
                    else:
                        try_value = ary[j][flag[j]]
                        if record_value > try_value:
                            record_value = try_value
                            jj = j
            sorted_array.append([jj, flag[jj]])
            flag[jj] += 1
            if flag[jj] >= len(ary[jj]):
                need_sort[jj] = False
        return sorted_array

    def sort_energy(self):
        #
        #
        control = self.ir
        sorted_list = []
        for ispin in range(self.nspin):
            ary = self.eigen_energy[ispin]
            sorted_list.append(self.sort_multi_array(ary, control))
        return sorted_list

    def band_energy(self, band_index, ispin=None):
        #
        # interface function to n'th band energy
        #
        # ispin = 1 :upsin
        #       = 2 :down-spin
        #
        # band_index: band index
        #
        #
        if ispin is None:
            if self.nspin == 1:
                ispin = 1
            else:
                print('===== Error(band_energy) in Estate =====')
                print('ispin is not defined')
                exit()
        else:
            if self.nspin == 1:
                ispin = 1
            else:
                if 1 <= ispin <= 2:
                    pass
                else:
                    print('===== Error(band_energy) in Estate =====')
                    print('ispin must be 1 or 2')
                    exit()

        band_energy_list = []
        for info in self.sort_energy()[ispin - 1]:
            #
            #  this implement ir is [0...maxir - 1]
            # and band index is 0, 1, 2,...
            #
            ir = info[0]
            e_index = info[1]
            ndeg = self.degeneracy[ir]
            energy = self.eigen_energy[ispin - 1][ir][e_index]
            for i in range(ndeg):
                band_energy_list.append(energy)
        max_band_index = len(band_energy_list)
        if band_index > max_band_index:
            print('===== Error(band_energy) in Estate =====')
            print('band index must be less than equal {}'.
                  format(max_band_index))
            exit()
        else:
            return band_energy_list[band_index - 1]
