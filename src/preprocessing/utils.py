import pandas as pd


def vartype(x,
            letters = ['A', 'T', 'C', 'G'],
            len_SV_lim = 100
            ):
    """
    Define the TYPE of a variant
    """
    if ">" in (x["REF"] + x["ALT"]) or "<" in (x["REF"] + x["ALT"]):
        return "SV"
    
    elif len(x["REF"]) > (len_SV_lim+1) or len(x["ALT"]) > (len_SV_lim+1) :
        return "SV"
    
    elif x["REF"] in letters and x["ALT"] in letters:
        return "SNV"
    
    elif len(x["REF"]) == len(x["ALT"]):
        return "MNV"
    
    elif x["REF"] == "-" or ( len(x["REF"]) == 1 and x["ALT"].startswith(x["REF"]) ):
        return "INSERTION"
    
    elif x["ALT"] == "-" or ( len(x["ALT"]) == 1 and x["REF"].startswith(x["ALT"]) ):
        return "DELETION"
    
    return "COMPLEX"
