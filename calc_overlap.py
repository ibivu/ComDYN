#!/usr/bin/python

# calc_overlap.py v0.1
# Copyright (c) 2018 K. Anton Feenstra (feenstra@few.vu.nl)
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

#import scipy.stats as stats
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#from matplotlib.backends.backend_pdf import PdfPages
from vararg import *

import numpy as np
from optparse import OptionParser

verbose=False

def parse_commandline():
    usage = "%prog <fasta> [<pdf>] [options]"
    description = \
        "%prog reads two xvg files and calculates the overlap between them " \
        "(assuming they are distributions based on the same bin sets)."
    epilog = \
        "Copyright (c) 2018 K. Anton Feenstra -- "\
        "feenstra@few.vu.nl -- www.few.vu.nl/~feenstra"
    parser = OptionParser(usage=usage, description=description,
                          epilog=epilog)
    # general options:
    parser.add_option("-v", "--verbose",   dest="verbose", action="store_true",
                      help="Verbose output (%default)")
    parser.set_defaults(verbose=False)
    # file options:
    parser.add_option("-r",  dest="ref", metavar="<file>",
                      help="reference input data file: xvg or x y text")
    parser.add_option("-f",  dest="xvgs", metavar="<file> [<file> ...]",
                      action="callback", callback=vararg_callback,
                      help="input data file(s) for which to calculate "\
                      "distance to reference: xvg or x y text")
    parser.set_defaults(xvgs=[])
    # parser.add_option("-o",  dest="out", metavar="<file>",
    #                   help="output plot file (pdf)")
    # parser.set_defaults(out=None)

    # get the options:
    (options, args) = parser.parse_args()

    print options.ref
    print options.xvgs
        
    if not options.ref:
        # check if we have an option left (to be used as input filename):
        try:
            options.ref=args.pop(0) # get first of remaining options
        except IndexError: # no more args to pop
            pass # ignore silently
        
    while len(options.xvgs)<1:
        # check if we have an option left (to be used as input filename):
        try:
            options.xvgs.append(args.pop(0)) # get first of remaining options
        except IndexError: # no more args to pop
            break # stop trying
    
    if not options.ref or not options.xvgs:
        print "ERROR: Need at least reference and one input file."
        print ""
        parser.print_help()
        print ""
        exit(-1)
    
    # if not options.out:
    #     # check if we have an option left (to be used as output filename):
    #     try:
    #         test = args.pop(-1) # get first of remaining options
    #         if test[-4:]=='.pdf':
    #             options.out=test #if pdf, keep as output name
    #         else:
    #             args.append(test) # else put back in args
    #     except IndexError:
    #         pass # not an error, this is optional!
    
    if len(args):
        print "WARNING: ignoring extra arguments:", args
    
    # clean up (recommended):
    del(parser)
    return options, args


def read_xvg(xvg, TYPE=None, maxcol=None):
    """ reads from xvg (either list, file or filename)
        and returns list of datasets.
        One dataset for each set in the input file, which is list of datapoints.
        Datapoints are list of (string) values
        type optional; it is either a single element or list, which will be
        applied to each of the data elements on each line. 
        Set maxcol to limit the number of columns being read. 
        """
    try:
        f_xvg=open(xvg,'r')
        if verbose: print "Opened", xvg
    except TypeError:
        # probably already file
        f_xvg=xvg

    list_of_sets=[]
    set_list=[]
    for line in f_xvg:
        if line[0] == '&':
            list_of_sets.append(set_list)
            set_list=[]
            continue
        
        if line[0] in "#@": 
            continue
        
        if not maxcol:
            words=line.split()
        else:
            words=line.split()[:maxcol]
        if not TYPE:
            vals=words
        else:
            try: # assume TYPE is list of type's:
                vals=[ t(w) for t,w in zip(TYPE, words) ]
            except TypeError: # else, it should be just one:
                vals=[ TYPE(w) for w in words ]
        set_list.append(vals)
    if set_list:
        if verbose: print "Appending last set", len(set_list)
        list_of_sets.append(set_list)
    
    f_xvg.close()
    
    return list_of_sets

def dict_find_keys(dic,val):
    r=[]
    for d in dic:
        if dic[d]==val:
            r.append(d)
    return r

if __name__ == "__main__":
    #Read data:
    options, args = parse_commandline()
    
    print "Reference", options.ref
    ds=read_xvg(options.ref, TYPE=float, maxcol=2)
    if options.verbose:
        print "Size of ref data:"
        print "sets:", len(ds)
        print "datapoints per set:", " ".join([str(len(d)) for d in ds])
        print "data depth per set:", " ".join([str(len(d[0])) for d in ds])

    # we may have multiple sets in the xvg: pack them together:
    d=[]
    map(d.extend, ds)
    # transpose set data into dict, with bins as keys:
    refd={b:c for b,c in d}
    refbins=set(refd.keys())
    reftot=sum(refd.values())
    refmax=max(refd.values())
    rmaxbins=dict_find_keys(refd, refmax)
    if options.verbose:
        print "Max", refmax, "at", rmaxbins
        print "bins:", refbins
        # print "values:", refd.values()
    
    for xvg in options.xvgs:
        print "Reading from", xvg
        ds=read_xvg(xvg, TYPE=float, maxcol=2)
        if options.verbose:
            print "Size of input data:"
            print "sets:", len(ds)
            print "datapoints per set:", " ".join([str(len(d)) for d in ds])
            print "data depth per set:", " ".join([str(len(d[0])) for d in ds])
        
        # we may have multiple sets in the xvg: pack them together:
        d=[]
        map(d.extend, ds)
        
        # transpose set data into dict, with bins as keys:
        d={b:c for b,c in d}
        bins=set(d.keys())
        tot=sum(d.values())
        dmax=max(d.values())
        dmaxbins=dict_find_keys(d, dmax)
        if options.verbose:
            print "Max", dmax, "at", dmaxbins
            print "bins:", bins
            #print "values:", refd.values()
        
        # get overlapping bins:
        binoverlap=list(refbins.intersection(bins))
        binoverlap.sort()
        
        if options.verbose:
            print binoverlap
        overlap=[ min(refd[b],d[b]) for b in binoverlap ] 
        if options.verbose:
            print overlap
        otot = sum(overlap)
        
        if options.verbose:
            print "ref", reftot, refmax, "set", tot, dmax, "overlap", otot
        print "Overlap:", otot/min(reftot,tot)
        print "Shift:", np.mean(dmaxbins)-np.mean(rmaxbins)

#last line
