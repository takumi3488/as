#!/usr/bin/env python3
#
# plot program for Logarithmic Derivative for KANSAI Pakage
#
#   written by Hiroki Funashima, 2001,2002,2003,2004,2005 Osaka University
#   modified by Hiroki Funashima, 2016,Kobe University
#    to rewrite for python3 style
#
#

class AtomicOrbital:
    def __init__(self):
        self.ia = 0
        self.l_list = []
        self.spin_list = []

    def nsp(self):
        return len(self.spin_list)

    def nl(self):
        return len(self.l_list)
