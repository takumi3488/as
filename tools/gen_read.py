#!/usr/bin/env python3
#
# gen_read.py
#   utility for `generator.table' for tspace 4th reviced version
#   written by Hiroki Funashima, 2000-2017 in Kobe
#
#   change log:
#     ver.1  written in perl
#     ver.2  rewritten in ruby
#     ver.3  include crystal database
#     ver.4  rewritten from scratch to refactored
#     ver.4.5  rewritten in python
#

from KansaiPackage.TSPACE.GenReadIO import GenReadIO

if __name__ == '__main__':
    generator_table = '/Volumes/exHD/funa/data/generator.table'
    crystal_database = '/Volumes/exHD/funa/data/gen_read_crystal.dbs'

    configure = {'generator_table': generator_table,
                 'crystal_database': crystal_database}
    gen_read = GenReadIO(configure)
    gen_read.header()
    gen_read.spg_set()
    if gen_read.spg is not None:
        gen_read.my_spg = gen_read.generator_table_obj.space_group(gen_read.spg)
        gen_read.show_space_group()
        gen_read.show_each_choice()
    gen_read.footer() 

