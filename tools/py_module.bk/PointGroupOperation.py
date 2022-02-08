#!/usr/bin/env python3
#
# symmetry operation for Oh group based on TSPACE
#  original code was written in Ruby by HF, 2016
#  written by HF in Kobe 25/07/2017
#
#


class PointGroupOperation:

    def __init__(self, il):
        self.xyz = ['x', 'y', 'z']
        if -1 <= il <= 4:
            self.il = il
        else:
            print("===== Error(PointGroupOperation) =====")
            print(" latice type IL is invalid.")
            print(" IL = {0}".format(il))
            exit()

        if self.il > 0:
            self.define_rotation_operator_Oh()
            self.define_operator_name_Oh()
        else:
            self.define_rotation_operator_D6h()
            self.define_operator_name_D6h()

    def clear(self):
        self.xyz = ['x', 'y', 'z']
        return self

    def operate(self, operator):
        if isinstance(operator, int):
            operator = self.index2operator(operator)
        elif isinstance(operator, str):
            if operator.isdigit():
                operator = self.index2operator(operator)
            else:
                operator = self.operator_name2operator(operator)

        new_xyz = [0, 0, 0]
        for i in range(0, 3):
            new_xyz[i] = self.__substitute(operator[i], self.xyz)

        for i in range(0, 3):
            self.xyz[i] = new_xyz[i]

        return self

    def unknown_operator(self):
        print("===== error =====")
        print("unknown operator")
        exit()

    def name(self):
        return self.operator2operator_name(self.xyz)

    def code_index(self):
        return self.operator2index(self.xyz)

    def operator2operator_name(self, operator):
        index = self.operator2index(operator)
        if index is None:
            print("==== internal error ====")
            print("operator = {0}".format(operator))
            exit()
        return self.index2operator_name(index)

    def operator_name2operator(self, operator):
        index = self.operator2index(operator)
        return self.index2operator(index)

    def operator2index(self, operator):
        import re
        if isinstance(operator, list):
            try:
                index = self.operator_list.index(operator)
                index = index + 1
            except ValueError:
                index = None

        elif isinstance(operator, str):
            if re.match('^[E, C, I]', operator.upper()):
                index = self.operator_name_list.index(operator.upper())
                if index is None:
                    self.unknown_operator()
                else:
                    return index + 1
            elif operator.isdigit():
                return int(operator)
            else:
                self.unknown_operator()
        else:
            self.unknown_operator()
        return index

    def index2operator(self, i):
        return self.operator_list[i - 1]

    def index2operator_name(self, i):
        return self.operator_name_list[i - 1]

    def matrix(self):
        matrix = [0, 0, 0]
        for i in range(0, 3):
            matrix[i] = self.__matrix_representation(self.xyz[i])
        return matrix

    def multiply(self, op1, op2):
        self.clear()
        for op in [op1, op2]:
            self.operate(op)
        return self

    def include_inv(self):
        import re
        if re.match('^I', self.operator2operator_name(self.xyz)):
            return True
        else:
            return False

    #
    # private member
    #
    def __substitute(self, operator, xyz):
        import re
        if not len(xyz) == 3:
            print('invalid xyz')
            print('size xyz = {0}'.format(len(xyz)))
            exit()

        if operator.strip().lower() in ['x', '-x']:
            element = xyz[0]
        elif operator.strip().lower() in ['y', '-y']:
            element = xyz[1]
        elif operator.strip().lower() in ['z', '-z']:
            element = xyz[2]
        elif operator.strip().lower() in ['w', '-w']:
            if [xyz[0], xyz[1]] == ['x', 'y']:
                element = 'w'
            elif [xyz[0], xyz[1]] == ['x', 'w']:
                element = 'y'

            elif [xyz[0], xyz[1]] == ['y', '-w']:
                element = 'x'
            elif [xyz[0], xyz[1]] == ['y', 'x']:
                element = '-w'

            elif [xyz[0], xyz[1]] == ['w', 'x']:
                element = '-y'
            elif [xyz[0], xyz[1]] == ['w', '-y']:
                element = 'x'

            elif [xyz[0], xyz[1]] == ['-x', '-y']:
                element = '-w'
            elif [xyz[0], xyz[1]] == ['-x', '-w']:
                element = '-y'

            elif [xyz[0], xyz[1]] == ['-y', 'w']:
                element = '-x'
            elif [xyz[0], xyz[1]] == ['-y', '-x']:
                element = 'w'

            elif [xyz[0], xyz[1]] == ['-w', '-x']:
                element = 'y'
            elif [xyz[0], xyz[1]] == ['-w', 'y']:
                element = '-x'

            else:
                print("===== Error ======")
                print("unknown type D6h operation for w")
                exit()
        else:
            print("===== Error ======")
            print("unknown operator {0}".format(operator.strip()))
            exit()

        if re.match('^-', operator):
            return self.__minus_element(element)
        else:
            return element

    def __minus_element(self, char):
        import re
        if re.match('^-', char):
            return char.replace('-', '')
        else:
            return '-' + char

    def __matrix_representation(self, code):
        if code.strip().lower() == 'x':
            return [1, 0, 0]
        elif code.strip().lower() == 'y':
            return [0, 1, 0]
        elif code.strip().lower() == 'z':
            return [0, 0, 1]
        elif code.strip().lower() == '-x':
            return [-1, 0, 0]
        elif code.strip().lower() == '-y':
            return [0, -1, 0]
        elif code.strip().lower() == '-z':
            return [0, 0, -1]
        else:
            self.unknown_operator

    def define_rotation_operator_Oh(self):
        self.operator_list = [
                ['x', 'y', 'z'],
                ['x', '-y', '-z'],
                ['-x', 'y', '-z'],
                ['-x', '-y', 'z'],
                ['z', 'x', 'y'],
                ['-z', 'x', '-y'],
                ['-z', '-x', 'y'],
                ['z', '-x', '-y'],
                ['y', 'z', 'x'],
                ['y', '-z', '-x'],
                ['-y', 'z', '-x'],
                ['-y', '-z', 'x'],
                ['y', 'x', '-z'],
                ['-y', '-x', '-z'],
                ['z', '-y', 'x'],
                ['-x', 'z', 'y'],
                ['-z', '-y', '-x'],
                ['-x', '-z', '-y'],
                ['x', '-z', 'y'],
                ['z', 'y', '-x'],
                ['-y', 'x', 'z'],
                ['x', 'z', '-y'],
                ['-z', 'y', 'x'],
                ['y', '-x', 'z'],
                ['-x', '-y', '-z'],
                ['-x', 'y', 'z'],
                ['x', '-y', 'z'],
                ['x', 'y', '-z'],
                ['-z', '-x', '-y'],
                ['z', '-x', 'y'],
                ['z', 'x', '-y'],
                ['-z', 'x', 'y'],
                ['-y', '-z', '-x'],
                ['-y', 'z', 'x'],
                ['y', '-z', 'x'],
                ['y', 'z', '-x'],
                ['-y', '-x', 'z'],
                ['y', 'x', 'z'],
                ['-z', 'y', '-x'],
                ['x', '-z', '-y'],
                ['z', 'y', 'x'],
                ['x', 'z', 'y'],
                ['-x', 'z', '-y'],
                ['-z', '-y', 'x'],
                ['y', '-x', '-z'],
                ['-x', '-z', 'y'],
                ['z', '-y', '-x'],
                ['-y', 'x', '-z']
                ]

    def define_operator_name_Oh(self):
        self.operator_name_list = [
                'E',
                'C2X',  'C2Y',  'C2Z',
                'C31+', 'C32+', 'C33+', 'C34+',
                'C31-', 'C32-', 'C33-', 'C34-',
                'C2A',  'C2B',  'C2C',  'C2D',  'C2E',  'C2F',
                'C4X+', 'C4Y+', 'C4Z+',
                'C4X-', 'C4Y-', 'C4Z-',
                'IE',
                'IC2X',  'IC2Y',  'IC2Z',
                'IC31+', 'IC32+', 'IC33+', 'IC34+',
                'IC31-', 'IC32-', 'IC33-', 'IC34-',
                'IC2A',  'IC2B',  'IC2C',  'IC2D',  'IC2E',  'IC2F',
                'IC4X+', 'IC4Y+', 'IC4Z+',
                'IC4X-', 'IC4Y-', 'IC4Z-'
                ]

    def define_rotation_operator_D6h(self):
        self.operator_list = [
                ['x', 'y', 'z'],
                ['w', 'x', 'z'],
                ['-y', 'w', 'z'],
                ['-x', '-y', 'z'],
                ['-w', '-x', 'z'],
                ['y', '-w', 'z'],
                ['-w', 'y', '-z'],
                ['x', 'w', '-z'],
                ['-y', '-x', '-z'],
                ['w', '-y', '-z'],
                ['-x', '-w', '-z'],
                ['y', 'x', '-z'],
                ['-x', '-y', '-z'],
                ['-w', '-x', '-z'],
                ['y', '-w', '-z'],
                ['x', 'y', '-z'],
                ['w', 'x', '-z'],
                ['-y', 'w', '-z'],
                ['w', '-y', 'z'],
                ['-x', '-w', 'z'],
                ['y', 'x', 'z'],
                ['-w', 'y', 'z'],
                ['x', 'w', 'z'],
                ['-y', '-x', 'z']
                ]

    def define_operator_name_D6h(self):
        self.operator_name_list = [
                'E',
                'C6+', 'C3+',
                'C2',
                'C6-', 'C3-',
                'C211', 'C221', 'C231',
                'C212', 'C222', 'C232',
                'IE',
                'IC6+', 'IC3+',
                'IC2',
                'IC6-', 'IC3-',
                'IC211', 'IC221', 'IC231',
                'IC212', 'IC222', 'IC232',
                ]
