#!/usr/bin/env python3
# -*- coding: utf-8 -*-
######################################################################################################

import argparse
import struct
filename_default = "sol.bin"

parser = argparse.ArgumentParser(description='Information for Solution File.')
parser.add_argument('filename', nargs='?', default=filename_default, help='Solution File')
parser.add_argument('-c','--config', help='Print CONFIG INFO',action='store_true')
parser.add_argument('-g','--glob', help='use project global path',action='store_true')


args =  parser.parse_args()


if (args.config):
    print ( "PRinting CONFIG")
    
if (args.glob):
    if (args.filename == filename_default): # damit auch mit global Flag die Input Datei geändert werden kann
        args.filename = "output/"+filename_default    

print ("using file:",args.filename)

data = open(args.filename, "rb").read(27)

print(struct.unpack("6siiii3s",data))

#parser.print_help()