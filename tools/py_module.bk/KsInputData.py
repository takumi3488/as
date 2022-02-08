#!/usr/bin/env python3
#
# plot program for Logarithmic Derivative for KANSAI Pakage
#
#   written by Hiroki Funashima, 2001,2002,2003,2004,2005 Osaka University
#   modified by Hiroki Funashima, 2016,Kobe University
#    to rewrite for python3 style
#
#
from py_module.AtomicOrbital import *
from py_module.ParseFort1 import *
class KsInputData:
    def __init__(self, filename, argvs):
        import os.path
        import re


        #
        # get fort1 info
        #
        if os.path.isfile('../fort.1'):
            fort1 = '../fort.1'
        elif os.path.isfile('fort.1'):
            fort1 = 'fort.1'
        elif os.path.isfile('wk/fort.1'):
            fort1 = 'wk/fort.1'
        else:
            print("fort.1 is not found.")
            quit()

        self.fort1Obj = ParseFort1(fort1)


        self.emin = -3.0
        self.emax =  4.0
        self.logdrv_min = -10.0
        self.logdrv_max = 20.0
        self.datafile = ''
        self.iter = 1
        self.nwindow = False
        self.rmin = 0.0
        self.rmax = 40.0

        if os.path.isfile(filename):
            self.filename = filename
            self.parse()

        if argvs.filename is not None:
            self.datafile = argvs.filename
        else:
            if self.datafile == '' :
                self.datafile = self.last_outfile()

        if argvs.iteration is not None:
            if re.match('^-*[0-9]+$', argvs.iteration ):
                self.iter = int(argvs.iteration)

        try:
            if argvs.magnetic is not None:
                self.spin_list = ['up', 'down']
            else:
                self.spin_list = []
        except AttributeError:
            pass

        try:
            if argvs.nwindow is not None:
                self.nwindow = True
            else:
                self.nwindow = False
        except AttributeError:
            pass

        try:
            if argvs.rmin is not None:
                self.rmin = float(argvs.rmin)
        except AttributeError:
            pass

        try:
            if argvs.rmax is not None:
                self.rmax = float(argvs.rmax)
        except AttributeError:
            pass

        self.count_niter()
        if self.iter == 0:
            self.iter = 1
        elif self.iter < 0:
            self.iter = self.niter - self.iter + 1
            if self.iter < 0:
                self.iter = 1
        elif self.iter > self.niter:
            self.iter = self.niter

        self.construct_dataset()


    def parse(self):
        import re
        self.dataset = []
        for line in open(self.filename, 'r'):
            if not re.search('^#', line.strip()) and \
                    not re.search('^$', line.strip()):
                dataline = line.strip().split('#')[0].strip()
                if re.search('=', dataline):
                    data_ary = line.strip().split('=')
                    key = data_ary[0].strip()
                    value = data_ary[1].strip()
                    if key.lower() == 'emin':
                        self.emin = float(value)
                    elif key.lower() == 'emax':
                        self.emax = float(value)
                    elif key.lower() == 'datafile':
                        self.datafile = value
                    elif key.lower() == 'iter':
                        self.iter = int(value)
                    elif key.lower() == 'logdrv_min':
                        self.logdrv_min = float(value)
                    elif key.lower() == 'logdrv_max':
                        self.logdrv_max = float(value)
                    elif key.lower() == 'nwindow':
                        if re.match('^y', value.lower()) or\
                            re.match('^t', value.lower()):
                            self.nwindow = True
                    elif key.lower() == 'rmin':
                        self.rmin = float(value)
                    elif key.lower() == 'rmax':
                        self.rmax = float(value)

    def construct_dataset(self):
        self.dataset = []
        for ika in range(self.fort1Obj.nka):
            obj = AtomicOrbital()
            obj.ia = ika + 1
            obj.l_list = ['s', 'p', 'd', 'f']
            try:
                obj.spin_list = self.spin_list
            except AttributeError:
                obj.spin_list = []
            self.dataset.append(obj)

    def count_niter(self):
        import re
        self.niter = 0
        for line in open(self.datafile, 'r'):
            if re.search('0\s=+\sITER=',line.strip()):
                self.niter += 1

    def last_outfile(self):
        import subprocess
        return subprocess.check_output('ls -ltr ../*.out | tail -1', shell=True).decode('utf-8').strip().split()[-1]
