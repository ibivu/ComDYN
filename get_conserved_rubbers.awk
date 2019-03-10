#!/usr/bin/awk -f

# get_conserved_rubbers.awk v0.1
# Copyright (c) 2017 K. Anton Feenstra (feenstra@few.vu.nl)
# 
# Syntax:
# get_conserved_rubbers.awk one.top two.top > list_of_rubber_bands.txt
# 
# Reads two gromacs topology file (.itp) with MARTINI moleculetype
# definition, extracts the 'rubber-band' interactions and prints those
# that are common between the two, including both reference distance
# and the difference between them.
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

# keep track of which file we're reading:
FNR==1 {
    filenr+=1;
}

# find the lines that define martini rubber bands to remember for the second file:
# we also remmeber the distance here, to compare to those in the the second file.
filenr==1 && /^[^#].*RUBBER_FC/ {
    d1nr++;
    anr1=$1
    anr2=$2
    d=$4
    dist[anr1,anr2]=d;
}

# in the second file, we'll process the rubber bands to find conserved ones:
filenr==2 && /^[^#].*RUBBER_FC/ {
    d2nr++;
    anr1=$1;
    anr2=$2;
    d2=$4;
    # find atom pair in dist list from first file:
    d1=dist[anr1,anr2];
    if (d1) {
	dboth++;
	print "DIST BOTH:", anr1, anr2, d2, d1, d2-d1;
    } else
	print "DIST SECOND:", anr1, anr2, d2;
}

#filenr==2 { print }

END {
    print "read", d1nr, "from first file"
    print "read", d2nr, "from second file"
    print "found", dboth, "in both files"
}

# last line
