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
