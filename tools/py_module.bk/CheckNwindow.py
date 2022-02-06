#!/usr/bin/env python3
#
#

from py_module.ParseFort1 import *
class CheckNwindow:
    def __init__(self):
        self.fort1_path = ['./fort.1', 'wk/fort.1', '../fort.1']
        self.find_fort1()
        self.check()

    def find_fort1(self):
        import os.path
        self.path_find = False
        for fort1_path in self.fort1_path:
            if os.path.exists(fort1_path):
                self.fort1 = fort1_path
                self.path_find = True
                break

    def check(self):
        if not self.path_find:
            return False

        self.max_comp = 0
        self.fort1_obj = ParseFort1(self.fort1)
        for ika in self.fort1_obj.valence:
            for l in ika:
                if len(l) > self.max_comp:
                    self.max_comp = len(l)

        if self.max_comp > 1:
            self.nwin =  True
        else:
            self.nwin = False

        return self.nwin
