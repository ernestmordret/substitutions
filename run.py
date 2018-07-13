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

parser.add_argument("-p", "--plot", action="store_true", 
                    help="plot the substitution heatmap, and\
                    output the unique substitutions count matrix.")
   
args = parser.parse_args()

if os.path.isdir(output_dir):
    pass
else:
    print('creating output directory')
    os.mkdir(output_dir)
    
if args.detect:
    if os.path.isfile(os.path.join(output_dir,'subs')):
        overwrite = raw_input("output directory already exists. Press y to overwrite")
        if overwrite == 'y':
            print('overwriting subs')
            try:
                os.remove(os.path.join(output_dir,'subs'))
                os.remove(os.path.join(output_dir,'subs.csv'))
            except OSError:
                pass
            os.mkdir(output_dir)
        else:
            print('exit')
            exit()

    print "detecting substitutions. That step might take some time..."
    try:
        import detect
    except:
        "this didn't go smoothly. Please check the parameters in the params.py file"
    
if args.quantify:
    if os.path.isfile(os.path.join(output_dir,'subs')):
        if os.path.isfile(os.path.join(output_dir,'qSubs')):
            overwrite = raw_input("output directory already exists. Press y to overwrite")
            if overwrite == 'y':
                print('overwriting qSsubs')
                try:
                    os.remove(output_dir+'/qSubs')
                except OSError:
                    exit()
        try:
            import quantify
        except:
            "this didn't go smoothly. Please check the parameters in the params.py file"
    else:
        print "subs not found. Please run detect first"
        
if args.plot:
    if os.path.isfile(os.path.join(output_dir,'subs')):
        try:
            import plot
        except:
            print "problems happened during the plotting..."
    else:
        print "subs not found. Please run detect first"