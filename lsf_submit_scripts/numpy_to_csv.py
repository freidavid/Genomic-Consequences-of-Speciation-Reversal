#!/bin/python
from optparse import OptionParser
...
parser = OptionParser()

parser.add_option('-i', '--infile', dest='infile', help='Input file name')
parser.add_option('-o', '--outfile', dest='outfile', help='Output file')

(options, args) = parser.parse_args()


import numpy as np
S = np.load(options.infile)
np.savetxt(options.outfile, S, delimiter=',')
