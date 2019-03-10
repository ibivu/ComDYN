#!/usr/bin/awk -f

# filter_conserved_rubbers.awk v0.1
# Copyright (c) 2017 K. Anton Feenstra (feenstra@few.vu.nl)
# 
# Syntax:
# filter_conserved_rubbers.awk list_of_residue_pairs input.top > output.top
# 
# Reads a list of residue pairs, and a gromacs topology file (.itp)
# with MARTINI moleculetype definition, and keeps only those
# 'rubber-band' interactions for the listed residue pairs.
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

BEGIN {
    npairs=0;
    nrubber=0;
    nkept=0;
    debug=0;
}

# keep track of which file we're reading:
FNR==1 {
    filenr+=1;
}

# find the lines that define martini rubber bands to remember for the second file:
# we also remmeber the distance here, to compare to those in the the second file.
filenr==1 {
    anr1=$1;
    anr2=$2;
    respairs[npairs]=anr1"-"anr2;
    npairs++;
    next;
}

# in the second file, we use the [atoms] section to find the bb bead atom number for each residue number

# beginning of topology block ('[ <name> ]') line:
/^[ \t]*\[/ && /\][ \t]*$/ {
  block = $0; 
  sub("^[ \t]*\[[ \t]*","",block);
  sub("[ \t]*\][ \t]*$","",block);
  tolower(block);
}

# atoms block: read them in, renumber, and print:
block=="atoms" {
    anr=$1;
    atp=$2;
    rnr=$3;
    rnm=$4;
    anm=$5;
    if ( anm=="BB") {
	restobbb[rnr]=anr;
	anrtores[anr]=rnr;
	if(debug) print ";DEBUG", anm, anr, rnr, restobbb[rnr];
    }
}

# at start of bonds block, make list of BB bead pairs:
/^[ \t]*\[[ \t]*bonds[ \t]*\][ \t]*$/ {
    for(i=0; i<npairs; i++) {
	n=split(respairs[i], respair, "-");
	if (n!=2) { print "DEATH HORROR:", i, respairs[i], n; exit(-1); }
	anr1=restobbb[respair[1]];
	anr2=restobbb[respair[2]];
	bbbpairs[anr1,anr2]=i;
	if(debug) print ";DEBUG", i, respairs[i], n, respair[1], respair[2], anr1, anr2
    }
}

# now for the 'rubber band' bonds, we do the filtering based on the bbb pairs:
block=="bonds" && /^[^#].*RUBBER_FC/ && NF>=5 {
    nrubber++;
    anr1=$1;
    anr2=$2;
    if(debug) print ";DEBUG", nrubber, anr1, anr2, bbpairs[anr1,anr2], nkept
    if ( (anr1,anr2) in bbbpairs)
	nkept++;
    else {
	r1=anrtores[anr1];
	r2=anrtores[anr2];
	printf "; SKIPPING dist %d (%d-%d)", r2-r1, r1, r2
    }
}

# now print everything; if something didn't need printing, we skipped it already.
# (the rubbers to be filtered out, already have a comment token on the output line.)
{ print }

#filenr==2 { print }

END {
    print ";read", npairs, "residue pairs"
    print ";found", nrubber, "rubbers in topology"
    print ";kept", nkept, "rubbers"
}

# last line
