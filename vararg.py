# run_mindist v0.1
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

def vararg_callback(option, opt_str, value, parser):
    """ Callback function for OptionParser to allow multiple arguments
    to a single option.
    """
    
    assert value is None
    done = 0
    value = []
    #print parser.option_list
    rargs = parser.rargs
    while rargs:
        arg = rargs[0]

        # Stop if we hit a known option arg, according to parser.has_option
        # we must also handle things like "-s0.2" and "--shcut=0.2"

        is_option = parser.has_option(arg)
        if not is_option and len(arg)>=2 and arg[0] == "-" and arg[1] != "-":
            # this may be something like "-s0.2",
            # check if first part ("-s") is option:
            is_option = parser.has_option(arg[:2])
        if not is_option and arg[:2] == "--":
            # this may be something like "--shcut=0.2",
            # check if part before '=' is option:
            is_option = parser.has_option(arg.split("=")[0])

        if is_option:
            break
        else:
            value.append(arg)
            del rargs[0]
    
    setattr(parser.values, option.dest, value)


def check_vararg_value(values):
    """ check if values is list; if only one string element,
    check if comma separated
    """

    if values and len(values)==1:
        try:
            words = values[0].split(',')
            if len(words)>1:
                values=words
        except AttributeError:
            pass
    
    return values
