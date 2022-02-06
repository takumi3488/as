#!/usr/bin/env python3
import os
import shutil
import subprocess
import re
def get_tmpdir():
    command = 'find . -name tmp -print'
    proc = subprocess.Popen(command,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True)
    stdout, stderr = proc.communicate()
    dir_list = stdout.decode('utf-8').split('\n')
    dir_list.pop(-1)
    return dir_list


for dirname in get_tmpdir():
    if os.path.isdir(dirname):
        ans = input('remove dir:{}? yes/no[no] > '.format(dirname)).strip()
        if re.search('^y', ans.lower()):
            shutil.rmtree(dirname)
