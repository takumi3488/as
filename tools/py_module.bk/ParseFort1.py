#!/usr/bin/env python3
#
# parse fort.1(crystal structure) for PointGroup
#   written by Hiroki Funashima, Kobe Univeristy, 2017
#
#
import math
import re


class ParseFort1:
    def __init__(self, filename):
        self.filename = filename
        self.symmorphic = False
        self.centrosymmetric = False
        self.generator_list = []
        self.translation_vector = []
        self.atom_name = []
        self.atomic_number = []
        self.orbital_info = []
        self.atomic_position = []
        self.parse()

    def parse(self):
        linecount = None
        self.centrosymmetric = False
        lineno = 0
        ika = 0
        try:
            zero_check = []
            for line in open(self.filename, 'r'):
                lineno += 1
                if lineno == 1:
                    self.title = line.strip()
                elif lineno == 2:
                    data = list(map(lambda x: int(x), line.strip().split()))
                    self.il = data[0]
                    self.ngen = data[1]
                    self.inv = data[2]
                elif 2 < lineno <= self.ngen + 2:
                    generator = int(line.strip().split()[0])
                    self.generator_list.append(generator)
                    data = []
                    for x in line.strip().split()[1:]:
                        data.append(int(x))
                    for i in range(3):
                        numerator = data[2*i]
                        denominator = data[2*i + 1]
                        gcd = self.gcd(numerator, denominator)
                        data[2*i] = data[2*i] // gcd
                        data[2*i+1] = data[2*i+1] // gcd

                    tvec = [[data[0], data[1]],
                            [data[2], data[3]],
                            [data[4], data[5]]
                            ]
                    if data[0] == 0 and data[2] == 0 and data[4] == 0:
                        zero_check.append(True)
                        if self.il > 0:
                            if generator == 25:
                                self.centrosymmetric = True
                        else:
                            if generator == 13:
                                self.centrosymmetric = True
                    else:
                        zero_check.append(False)
                    self.translation_vector.append(tvec)
                elif lineno == self.ngen + 2 + 1:
                    self.a, self.b, self.c = self.get_3value(line)
                elif lineno == self.ngen + 2 + 2:
                    self.ca, self.cb, self.cc = self.get_3value(line)
                elif lineno == self.ngen + 2 + 3:
                    data = line.split()
                    data = list(map(lambda x: int(x), line.split()))
                    self.nka = int(data[1])
                    for i in range(self.nka):
                        self.orbital_info.append([])
                    self.na = int(data[0])
                    self.ka = []
                    ka_part = True
                    coordinate_part = True
                    for i in range(2, len(data)):
                        self.ka.append(data[i])
                elif lineno > self.ngen + 2 + 3:
                    if re.match('^[A-Z,a-z]', line.strip()):
                        linecount = 0
                        self.atom_name.append(line.strip().split()[0])
                        z = float(line.strip().split()[1].replace('D', 'e'))
                        self.atomic_number.append(z)
                        ika += 1
                    else:
                        if linecount is None:
                            if not re.search('D', line.strip()):
                                if ka_part:
                                    for i in list(map(lambda x: int(x),
                                                      line.strip())):
                                            self.ka.append(i)
                                else:
                                    coordinate_part = False
                            else:
                                ka_part = False
                                if coordinate_part:
                                    self.atomic_position.\
                                        append(self.get_3value(line))

                        else:
                            linecount += 1
                            if linecount == 1:
                                self.orbital_info[ika - 1] = \
                                        self.line2orbital(line)
                            elif linecount == 2 and \
                                    not re.search('D', line.strip()):
                                for x in self.line2orbital(line):
                                    self.orbital_info[ika - 1].append(x)

                if all(zero_check):
                    self.symmorphic = True
                else:
                    self.symmorphic = False

        except FileNotFoundError:
            print("===== Error(ParseFort1) ======")
            print(" file:{} is not found.".format(self.filename))
            exit()

    def get_3value(self, line):
        value = []
        for x in line.split():
            for xx in self.devide_2floats(x):
                value.append(float(xx.replace('D', 'e')))
        return value

    def devide_2floats(self, data):
        array = []
        if re.search('0-0', data):
            data_list = data.split('0-0')
            for i, v in enumerate(data_list):
                element = v
                if i > 0:
                    element = '-0' + element
                if i < len(data_list) - 1:
                    element = element + '0'
                array.append(element)
        if len(array) == 0:
            return [data]
        else:
            return array

    def line2orbital(self, line):
        return list(map(lambda x: int(x), line.strip().split()))

    def show_info(self):
        print(' ** infomation for fort.1:')
        print('   lattice type: ', end='')
        if self.il == -1:
            lattice_type = 'rhombohedral'
        elif self.il == 0:
            lattice_type = 'hexagonal'
        elif self.il == 1:
            lattice_type = 'sigle'
        elif self.il == 2:
            lattice_type = 'face centered'
        elif self.il == 4:
            lattice_type = 'c-center base centered'
        else:
            print(' unknown lattice type(error)')
            exit()
        print("{0} lattice".format(lattice_type))
        print("   Num of generator = {0}".format(self.ngen))
        print("     -> rotational operation code in tspace: ", end='')
        order = len(self.generator_list)
        for i, generator in enumerate(self.generator_list):
            print('{}'.format(generator), end='')
            if order > 1:
                if i < order - 2:
                    print(', ', end='')
                elif i == order - 2:
                    print(' and ', end='')
                else:
                    pass
        print()
        print()
        print("  lattice constant")
        print("    a = {0:.5f},  b = {1:.5f},  c={2:.5f}".
              format(self.a, self.b, self.c))
        print("   ca = {0:.5f}, cb = {1:.5f}, cc={2:.5f}".
              format(self.ca, self.cb, self.cc))
        print()
        print("   num of kind of atoms = {0},   num of atoms = {1}".
              format(self.nka, self.na))
        for i, atom in enumerate(self.atom_name):
            print('     ({0:1d}):{1:2s}'.format(i+1, atom))
        print()

    def show_atom_info(self):
        for atom in self.atom_name:
            print(atom)

    def gcd(self, a, b):
        if a * b == 0:
            return 1
        else:
            za = abs(a)
            zb = abs(b)
            return math.gcd(za, zb)
