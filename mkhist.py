#!/usr/bin/python

#### mkhist ####
# (c) 1998-2018 Anton Feenstra k.a.feenstra@vu.nl
# 
# based on the template for the Umbrella Sampling MD practical
# of the Structural Bioformatics course 2011/2012
# which was in turn based on mkhist.awk, from the previous century.
#
# This may look like software, but is an actual research tool.
# It workes for the purpose it was written for. It may work for you.
# If not and/or it ends up destroying your universe and all you hold
# dear, I will feel sorry for you, but you cannot hold me responsible.
# Please do file a bug report :-)
# 
# newest version(s) available from:
# https://github.com/ibivu/ConsDYN.git 
# 
#    PLEASE CITE:
#    
#    Halima Mouhib, Akiko Higuchi, Sanne Abeln, 
#    Kei Yura, K. Anton Feenstra. 
#    "Showing the impact of pathogenic mutations of the glucose
#     transporter (GLUT1) on the channel dynamics using ConsDYN"
#    F1000, submitted (2019).
#    
#    :ETIC ESAELP
#    
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import math
from optparse import OptionParser, OptionGroup

# we will use a few global variables, because they're really global:
kB   = 8.3144621e-3  # boltzmann's constant in kJ/mol/K
binw = 0             # binwidth; will be set according to optparse default

######## COMMAND LINE / INPUT STUFF ##########

def parse_commandline():
    usage = "%prog <files> [options]"
    description = \
        "%prog processes data from umbrella sampling simulations." \
        "The energy file is expected to have time, umbrella energy, " \
        "and (any number of) potential energy (terms) white-space " \
        "separated and one line per time point." \
        "The distance file is expected to have time and distance, " \
        "also white-space separated and one line per time point."
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("-b", "--binw", dest="binw", type="float", 
                     help="Binwidth (%default)")
    parser.set_defaults(binw=10)
    parser.add_option("-d", "--data-col", dest="datacol", type="int",
                     help="column to take data from (%default)")
    parser.set_defaults(datacol=2)
    parser.add_option("-w", "--weight-col", dest="weightcol", type="int",
                     help="optional column to take weights from (%default)")
    parser.set_defaults(weightcol=None)
    parser.add_option("", "--bars", dest="drawbars", 
                      action="store_true",
                      help="Plot bar outlines explicitly (%default)")
    parser.set_defaults(drawbars=False)
    parser.add_option("-x", "--exit-on-error", dest="exitonerror", 
                      action="store_true",
                      help="Exit on error parseing input (%default)")
    parser.set_defaults(exitonerror=False)
    
    # get the options:
    (options, args) = parser.parse_args()
    
    if len(args)>1:
        print "Error: only one input file supported"
        sys.exit(-1)
        
    # clean up (recommended):
    del(parser)
    return options, args[0]

######## MISC FUNCTIONS ##########

def bin2data(b):
    return (b+0.5)*binw

def bin2datarange(b):
    return (b)*binw, (b+1)*binw

def data2bin(t):
  return int(t/binw)

######## MAIN ##########
# I use a 'main' function to avoid accidental (global) use 
# of variables defined in the main code block.
def main():
    global binw
    
    options, infile = parse_commandline()
        
    # set global bin width
    binw = options.binw
    
    # reset counting:
    minbin = 1e99
    maxbin = -1e99
    mind = 1e999
    maxd = -1e999
    sumd = 0 # sum of data (weighted)
    sum2 = 0 # sum of squares (weighted)
    sumw = 0 # sum of weights (nr of datapoints if unweighted)
    n    = 0 # nr of datapoints for reporting
    w    = 1 # default weight=1 if not reading from input
    n_bin= {}# hist data (dict, not list, so we can index)
    
    # we want one-pass weighted stdev (\sigma^2):
    #     \sigma^2 = 1/W \sum_i w_i x_i^2 - \mu^2
    # where:
    #     W = \sum_i w_i
    #     \mu = 1/W \sum_i w_i x_i
    # so we collect weighted sum (W) and weighted sum of squares
    
    # read energies from file, line by line:
    for line in open(infile):
        if ( line.startswith("#") or
             line.startswith("@") or
             line.startswith("&") ) : continue      # skip some lines
        words = line.split()
        try:
            d = float(words[options.datacol])
            if options.weightcol:
                w = float(words[options.weightcol])
        except ValueError:
            if options.exitonerror:
                print "Illegal input on line", "'"+line.rstrip()+"'"
                sys.exit(-1)
            else:
                continue # with next line
        
        # collect stats:
        mind = min(mind,d)
        maxd = max(maxd,d)
        sumd += w*d
        sum2 += w*d**2
        sumw += w
        n    += 1
        
        # calculate bins:
        bin = data2bin(d)
        minbin = min(minbin,bin)
        maxbin = max(maxbin,bin)
        
        # make sure each bin gets initialized:
        n_bin.setdefault(bin, 0)
        
        # now add up numbers for each bin:
        n_bin[bin] += w

    # calc stats (formula above for line loop):
    aver  = sumd/sumw
    stdev = math.sqrt(sum2/sumw - aver**2)
    print "# nr data points:", n
    print "# total weight:", sumw
    print "# min value:", mind
    print "# max value:", maxd
    print "# average value:", aver
    print "# standard deviation (sigma^2):", stdev
    
    # print hist data
    if options.drawbars:
        print "# bin", "N", "N/unit", "norm.dist."
        for i in sorted(n_bin.keys()):
            binstart,binend = bin2datarange(i)
            if not n_bin.has_key(i-1):
                print binstart, 0, 0, 0
            print binstart, n_bin[i], n_bin[i]/binw, n_bin[i]/sumw
            print binend, n_bin[i], n_bin[i]/binw, n_bin[i]/sumw
            if not n_bin.has_key(i+1):
                print binend, 0, 0, 0
    else:
        print "# bin", "N", "N/unit", "norm.dist."
        for i in sorted(n_bin.keys()):
            print bin2data(i), n_bin[i], n_bin[i]/binw, n_bin[i]/sumw

if __name__ == "__main__":
    main()

#last line
