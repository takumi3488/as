#!/usr/bin/env python3
#
# Generate point group from generator
#    written by HF in Kobe on 28/07/2017
#
from py_module.PointGroupOperation import *
from py_module.PointGroupName import *
import re


class PointGroup:
    def __init__(self, il, generator_list):
        self.il = il
        self.generator_list = generator_list
        self.element_list = []
        for x in generator_list:
            self.element_list.append(x)
        self.add_indentity_operation()
        self.pg_obj = PointGroupOperation(self.il)
        self.generate_all_elements()
        self.pg_name_obj = PointGroupName(self)
        self.name = self.pg_name_obj.pg_name
        self.crystal_system = self.crystal_system_name()
        self.bravais_lattice = self.bravais_lattice_name()

    def crystal_system_name(self):
        if self.pg_name_obj.crystal_system == -1:
            return 'rhombohedral'
        elif self.pg_name_obj.crystal_system == 0:
            return 'hexagonal'
        elif self.pg_name_obj.crystal_system == 1:
            return 'cubic'
        elif self.pg_name_obj.crystal_system == 2:
            return 'tetragonal'
        elif self.pg_name_obj.crystal_system == 3:
            return 'orthorhombic'
        elif self.pg_name_obj.crystal_system == 4:
            return 'monoclinic'
        elif self.pg_name_obj.crystal_system == 5:
            return 'triclinic'
        else:
            print("====== Error(PointGroupName) ======")
            print(" unknown crystal system")
            print(self.pg_name_obj.crystal_system)
            exit()

    def bravais_lattice_name(self):
        if self.il <= 0:
            prefix = ''
        elif self.il == 1:
            prefix = 'simple'
        elif self.il == 2:
            prefix = 'face centered'
        elif self.il == 3:
            prefix = 'body centered'
        elif self.il == 4:
            prefix = 'bace centered'
        else:
            print('unknown lattice type')
            print('il = {0}'.format(self.il))
            exit()
        return prefix + ' ' + self.crystal_system + ' lattice'


    def show_generator_list_info(self):
        print('   ========= Generator List =========')
        self.show_list_info(self.generator_list)

    def show_element_list_info(self):
        print('   ========= Element List =========')
        self.show_list_info(self.element_list)

    def show_list_info(self, list_name):
        print("   Num of elements = {}".format(len(list_name)))
        print("   code   name       operation")
        print("   -------------------------------")
        for iop in list_name:
            print('   ',end='')
            self.element_info(iop)
        print()
        print("    -------------------------------")
        print()

    def element_info(self, iop):
        self.pg_obj.clear()
        self.pg_obj.operate(iop)
        icode = self.pg_obj.code_index()
        code_name = self.pg_obj.name()
        operator = self.pg_obj.index2operator(icode)
        print("  {0:2d}  ".format(icode), end='')
        if re.search('^I', code_name):
            print("{0:<5s}  ".format(code_name), end='')
        else:
            print(" {0:<4s}  ".format(code_name), end='')
        print("[", end='')
        for x in operator:
            print(" {0:>2s} ".format(x), end='')
        print("]")
        return self

    def add_indentity_operation(self):
        identity_operation = 1
        for list_name in [self.generator_list, self.element_list]:
            if identity_operation not in list_name:
                list_name.insert(0, identity_operation)

    def generate_all_elements(self):
        #
        #
        self.old_order = len(self.element_list)
        self.new_order = 0
        while self.new_order is not self.old_order:
            self.old_order = len(self.element_list)
            add_list = []
            for i in range(len(self.element_list)):
                for j in range(len(self.element_list)):
                    op1 = self.element_list[i]
                    op2 = self.element_list[j]
                    iop = self.multiply(op1, op2)
                    if iop not in add_list:
                        add_list.append(iop)

                for j in range(len(self.element_list)):
                    op1 = self.element_list[i]
                    op2 = self.element_list[j]
                    iop = self.multiply(op2, op1)
                    if iop not in add_list:
                        add_list.append(iop)

            for iop in add_list:
                self.add_new_element(iop)

            self.new_order = len(self.element_list)

    def multiply(self, op1, op2):
        self.pg_obj.clear()
        iop = self.pg_obj.multiply(op1, op2).code_index()
        self.pg_obj.clear()
        return iop

    def add_new_element(self, element):
        if element not in self.element_list:
            self.element_list.append(element)
            self.element_list.sort()
