#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 11:32:57 2018

@author: ernestmordret
"""


import argparse
from params import *
import os

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--detect", action="store_true", 
                    help="run the detect script only.")

parser.add_argument("-q", "--quantify", action="store_true", 
                    help="run the quantify script only. Will not work\
                        unless the detect option was ran first")
   
args = parser.parse_args()

if os.path.isdir(output_dir):
    pass
else:
    print('creating output directory')
    os.mkdir(output_dir)
    
if args.detect:
    if os.path.isfile(output_dir+'/subs'):
        overwrite = raw_input("output directory already exists. Press y to overwrite")
        if overwrite == 'y':
            print('overwriting subs')
            try:
                os.remove(output_dir+'/subs')
                os.remove(output_dir+'/subs.csv')
                os.remove(output_dir+'/dp')
            except OSError:
                pass
            os.mkdir(output_dir)
        else:
            print('exit')
            exit()

    print "detecting substitutions. That step might take some time..."
    import detect
    
if args.quantify:
    if os.path.isfile(output_dir+'/subs'):
        if os.path.isfile(output_dir+'/qSubs'):
            overwrite = raw_input("output directory already exists. Press y to overwrite")
            if overwrite == 'y':
                print('overwriting qSsubs')
                try:
                    os.remove(output_dir+'/qSubs')
                except OSError:
                    exit()
        import quantify
    else:
        print "subs not found. Please run detect first"