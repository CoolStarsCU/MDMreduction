#!/usr/bin/env python
#!/user/covey/iraf/mypython

###############################################################################
#split_strings.py
#
# Split up a string into a list of substrings separated by arbitrary numbers
# of spaces or tabs.
#
#Kevin Covey
#Version 1.0 (7-23-03)
###############################################################################

import string

def split_strings(inputstring):

# cut off the endline character

    noendline = inputstring[:-1]

    #print noendline

# split the string up along any spaces we can find
    
    nospaces = string.split(noendline, ' ')

    #print nospaces
    
# now cycle through each item in each list and split on any tabs.

    newitems = len(nospaces)
    k = 0
    tabseplist = ['crap']
    while k < newitems:
        gottabs = nospaces[k]
        splitup = string.split(gottabs,'\t')
        newentries = len(splitup)
        l = 0
        while l < newentries:
            tabseplist.append(splitup[l])
            l = l + 1
        k = k + 1

    #print tabseplist

    tabseplist = tabseplist[1:]

# get rid of any empty items in the list

    items = len(tabseplist)
    i = 0
    j = 0
    while i < items:
        if tabseplist[j] == '':
            del tabseplist[j]
            i = i + 1
            j = j
        else:
            i = i + 1
            j = j + 1

    #print tabseplist
        
    return tabseplist
    
