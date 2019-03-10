#!/usr/bin/python

# plot_d1d2.py v0.1
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

import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
from vararg import *

import numpy as np
from optparse import OptionParser

def parse_commandline():
    usage = "%prog <fasta> [<pdf>] [options]"
    description = \
        "%prog reads xvg files and plots contour densities based on "\
        "data point density from the input. Densities are estimated "\
        "using KDE. Multiple datasets in the xvg input are separated "\
        "by a line containing a single '&'; these will be concatenated "\
        "Datasets of one length will not be included in the densities "\
        "but instead be plotted as single larger dot."
    epilog = \
        "Copyright (c) 2018 K. Anton Feenstra -- "\
        "feenstra@few.vu.nl -- www.few.vu.nl/~feenstra"
    parser = OptionParser(usage=usage, description=description,
                          epilog=epilog)
    # file options:
    parser.add_option("-v", "--verbose",   dest="verbose", action="store_true",
                      help="Verbose output (%default)")
    parser.set_defaults(verbose=False)
    parser.add_option("-f",  dest="xvg", metavar="<file> [<file> ...]",
                      action="callback", callback=vararg_callback,
                      help="input data file(s): xvg or x y text")
    parser.set_defaults(xvg=None)
    parser.add_option("-o",  dest="out", metavar="<file>",
                      help="output plot file (pdf)")
    parser.set_defaults(out=None)
    parser.add_option("-t",  dest="title", metavar="<title>",
                      help="title for plot")
    parser.set_defaults(xlabel=None)
    parser.add_option("-x",  dest="xlabel", metavar="<label>",
                      help="label for plot x axis")
    parser.set_defaults(xlabel=None)
    parser.add_option("-y",  dest="ylabel", metavar="<label>",
                      help="label for plot y axis")
    parser.set_defaults(ylabel=None)
    parser.add_option("", "--xrange",  dest="xrange", metavar="<start,end>",
                      help="x range for KDE grid (default based on data)")
    parser.set_defaults(xrange=None)
    parser.add_option("", "--yrange",  dest="yrange", metavar="<start,end>",
                      help="y range for KDE grid (default based on data)")
    parser.set_defaults(yrange=None)
    parser.add_option("-g", "--grid",  dest="gridspace", metavar="<g>", type=float,
                      help="spacing for KDE grid (default 100 grid points)")
    parser.set_defaults(gridspace=None)

    # get the options:
    (options, args) = parser.parse_args()
    
    print args
    if not options.xvg:
        # check if we have an option left (to be used as input filename):
        try:
            options.xvg = [args.pop(0)] # get first of remaining options
        except IndexError:
            print "Need at least an input file (xvg or x y text)"
            print ""
            parser.print_help()
            print ""
            exit(-1)
    
    if not options.out:
        # check if we have an option left (to be used as output filename):
        try:
            test = args.pop(-1) # get first of remaining options
            if test[-4:]=='.pdf':
                options.out=test #if pdf, keep as output name
            else:
                args.append(test) # else put back in args
        except IndexError:
            pass # not an error, this is optional!
    
    if len(args): # assume they're all input files:
        if not options.xvg: options.xvg=[]
        options.xvg.extend(args)
    
    # clean up (recommended):
    del(parser)
    return options, args


def read_xvg(xvg, TYPE=None):
    """ reads from xvg (either list, file or filename)
        and returns list of datasets.
        One dataset for each set in the input file, which is list of datapoints.
        Datapoints are list of (string) values
        type optional; it is either a single element or list, which will be
        applied to each of the data elements on each line. 
        """
    try:
        f_xvg=open(xvg,'r')
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
        
        words=line.split()
        if not TYPE:
            vals=words
        else:
            try: # assume TYPE is list of type's:
                vals=[ t(w) for t,w in zip(TYPE, words) ]
            except TypeError: # else, it should be just one:
                vals=[ TYPE(w) for w in words ]
        set_list.append(vals)

    f_xvg.close()
    
    return list_of_sets

    
if __name__ == "__main__":
    #Read data:
    options, args = parse_commandline()
    
    # set Plot ready for the results:
    fig, ax = plt.subplots()

    if options.title: ax.set_title(options.title)
    if options.xlabel: ax.set_xlabel(options.xlabel)
    if options.ylabel: ax.set_ylabel(options.ylabel)
    
    print "Reading from", options.xvg
    colors=cm.get_cmap("Set1").colors[2:]
    if options.verbose:
        print colors[0]
        print colors[1]
        print colors[2]
    linewidths=[1,1,3]
    # build list of contour levels, max should be at 200, logarithmic scale:
    contour_levels = [ 200*(2**l) for l in range(-9,1) ]
    
    for set, xvg in enumerate(options.xvg):
        ds=read_xvg(xvg, TYPE=float)
        if options.verbose:
            print "Size of input data:"
            print "sets:", len(ds)
            print "datapoints per set:", " ".join([str(len(d)) for d in ds])
            print "data depth per set:", " ".join([str(len(d[0])) for d in ds])
        # we may have multiple sets in the xvg: pack them together:
        d=[]
        map(d.extend, ds)
        if len(ds[-1])==1: # if last set is only one, remove it:
            s=d.pop(-1) # and keep it for plotting
        if options.verbose:
            print "concatenated sets:"
            print len(d)
            print len(d[0])
        # transpose set data into separate x and y arrays (m1, m2):
        m1,m2 = map(list, zip(*d))
        if options.verbose:
            print "transposed data:", len(m1), len(m2)

        # find data ranges; add 0.05 (nm) on all sides:        
        if not options.xrange:
            xmin = min(m1)-0.05
            xmax = max(m1)+0.05
        else:
            xmin, xmax = [ float(w) for w in options.xrange.split(",") ]
        if not options.yrange:
            ymin = min(m2)-0.05
            ymax = max(m2)+0.05
        else:
            ymin, ymax = [ float(w) for w in options.yrange.split(",") ]
        # # expand ranges by 10% on all sides:
        # xr=xmax-xmin; xmin-=0.1*xr; xmax+=0.1*xr
        # yr=ymax-ymin; ymin-=0.1*yr; ymax+=0.1*yr
        if options.verbose:
            print "x from", xmin, "to", xmax
            print "y from", ymin, "to", ymax
        
        #Perform a kernel density estimate on the data:
        
        if not options.gridspace:
            space=100j
        else:
            space=options.gridspace
        X, Y = np.mgrid[xmin:xmax:space, ymin:ymax:space]
        positions = np.vstack([X.ravel(), Y.ravel()])
        values = np.vstack([m1, m2])
        kernel = stats.gaussian_kde(values)
        Z = np.reshape(kernel(positions).T, X.shape)
        if options.verbose:
            print Z.shape
            print "z from", Z.min(), Z.max()
        
        # now plot:
        #    ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
        #              extent=[xmin, xmax, ymin, ymax])
        ax.plot(m1, m2, 'k.', markersize=2, color=colors[set])
        ax.plot(s[0], s[1], 'k.', markersize=10, color=colors[set])
        ax.contour(X, Y, Z, contour_levels, 
                   colors=[colors[set]],
                   linewidths=linewidths[min(set,len(linewidths)-1)])
        # only enlarge plot when new sets are added:
        xmin=min(ax.get_xlim()[0],xmin)
        xmax=max(ax.get_xlim()[1],xmax)
        ymin=min(ax.get_ylim()[0],ymin)
        ymax=max(ax.get_ylim()[1],ymax)
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
        ax.set_aspect('equal')
    if not options.out:
        plt.show()
    else:
        pp=PdfPages(options.out)
        pp.savefig(plt.gcf())
        pp.close()

#last line
