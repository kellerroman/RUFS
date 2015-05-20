#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#######################
import struct

from configparser import SafeConfigParser

parser = SafeConfigParser()
parser.read('config.cfg')

iter        = 1
solution    = 1
residual    = 1
equation    = 1
dimension   = 1
space_disc  = 1
space_order = 1
CFL         = 0.1
timestep    = 1E-4

for section_name in parser.sections():
    print ('Section:', section_name)
#     print ('  Options:', parser.options(section_name))
    for name, value in parser.items(section_name):
        print ('  %s = %s' % (name, value))
        if (section_name == "control" and name == "iterations"):
            iter = int(value)
        if (section_name == "output" and name == "solution"):
            solution = int(value)
        if (section_name == "output" and name == "residual"):
            residual = int(value)
        if (section_name == "disc" and name == "timestep"):
            timestep = float(value)
        if (section_name == "disc" and name == "CFL"):
            CFL = float(value)
#     print ()
    
# Write binary data to a file
with open('config.bin', 'wb') as f:
    f.write((1)             .to_bytes(4, byteorder='little', signed=True))
    f.write(iter            .to_bytes(4, byteorder='little', signed=True))
    f.write(solution        .to_bytes(4, byteorder='little', signed=True))
    f.write(residual        .to_bytes(4, byteorder='little', signed=True))
    f.write(equation        .to_bytes(4, byteorder='little', signed=True))
    f.write(dimension       .to_bytes(4, byteorder='little', signed=True))
    f.write(space_disc      .to_bytes(4, byteorder='little', signed=True))
    f.write(space_order     .to_bytes(4, byteorder='little', signed=True))
    f.write(struct.pack('d',CFL))
    f.write(struct.pack('d',timestep))