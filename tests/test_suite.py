#!/usr/bin/env python
'''
Run tests
'''

import os
import sys
import glob
import shutil
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
from subprocess import call
from Bio import SeqIO

#$sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def load_data_loc():
    return os.path.join(str(os.getcwd()), \
        'test_data/')

def load_random_test_dir():
    loc = os.path.join(str(os.getcwd()), \
        'test_backend/testdir/')
    return loc

def get_script_loc(script):
    if script == 'rpxsuite':
        return os.path.join(str(os.getcwd()), \
            '../rpXsuite/RPxSuite.py')

class test_RPxSuite():
    def setUp(self):
        self.script = get_script_loc('rpxsuite')
        self.test_dir = load_random_test_dir()
        self.fastas = glob.glob(load_data_loc() + \
            '*.fasta')

        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def tearDown(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

    def run(self):
        self.setUp()
        self.test0()
        self.tearDown()

    def test0(self):
        '''
        Test that it runs at all
        '''
        # Run program
        base = self.test_dir + 'test'

        cmd = "{0} {1} -g Ribosomal_L6 -o {2}".format(self.script, ' '.join(self.fastas),
            self.test_dir)
        print(cmd)
        call(cmd, shell=True)

        # Make sure it produces files
        assert len(glob.glob(self.test_dir + '*')) == 11


if __name__ == '__main__':
    test_RPxSuite().run()
    print('everything is working swimmingly!')
