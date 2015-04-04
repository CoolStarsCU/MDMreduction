#! /usr/bin/env python
#! /user/covey/iraf/mypython

############################################################################
#
# sky_checker.py
#
# perform sky subtraction on a spectrum, then check to make sure that the
# user thinks a good job has occured by displaying the spectra and asking
# for confirmation of goodness.  If the spectra isn't good enough, allow
# user to define a number of regions to find the best scaling for optimal
# sky subtraction.
#
#Kevin R. Covey
#Version 1.0 (9-9-03)
############################################################################

from pyraf import iraf
from pyraf.iraf import images, imutil, tv
import split_strings, os
from split_strings import split_strings

def sky_checker(readprefix,target,skysubtractor,writeprefix,core):

# tweak the display parameters

    iraf.unlearn(iraf.images.tv.display)
    iraf.images.tv.display.frame = '1'
    iraf.images.tv.display.fill = 'yes'

# perform the initial sky subtraction (check for signs of a scaled version
# of this image, which saves us the trouble of perfecting the subtraction)

    step = .25
    scale = 1
    scaleup = 3

    previous = 0
    previous = os.access('skysub/scaled.' + target, os.F_OK)

    if previous == 1:
        iraf.imarith(operand1 = 'skysub/scaled.' + target, op = '-', operand2 = readprefix + core + skysubtractor + '.fits', result = writeprefix + target)
    elif previous ==0:
        iraf.imarith(operand1 = readprefix + target, op = '-', operand2 = readprefix + core + skysubtractor + '.fits', result = writeprefix + target)
        iraf.imarith(operand1 = readprefix + core + skysubtractor + '.fits', op = '*', operand2 = scale, result = 'skysub/scaled.' + core + skysubtractor + '.fits')

# display to make sure user accepts the result

    iraf.images.tv.display(image = writeprefix + target, frame = '1')
    done = raw_input('Are skylines removed acceptably? (y or n) ')

# enter a loop to maximize subtraction

    definedregions = 'n'
    while done == 'n':

# delete result of the last attempt at subtraction
        os.remove(writeprefix + target)
        os.remove('skysub/scaled.' + core + skysubtractor + '.fits')
        
# do subtraction with current scale value
        iraf.imarith(operand1 = readprefix + core + skysubtractor + '.fits', op = '*', operand2 = scale, result = 'skysub/scaled.' + core + skysubtractor + '.fits')
        iraf.imarith(operand1 = readprefix + target, op = '-', operand2 = 'skysub/scaled.' + core + skysubtractor + '.fits', result = writeprefix + target)

# display results, if only for testing purposes
        iraf.images.tv.display(image = writeprefix + target, frame = '1')

# if it hasn't already happened, define regions for statistical tests
        if definedregions == 'n':
             definedregions = 'y'
             nregions = int(raw_input('how many lines shall we optimize subtraction for? ' ))

             skyxs = [0]
             skyys = [0]
             linexs = [0]
             lineys = [0]

             i = 0
             while i <= nregions-1:
                 print 'define lower left corner of blank sky region: (place cursor, press x, then q)'
                 llsky = iraf.images.tv.imexamine(input = writeprefix + target, Stdout = 1)
                 print 'define upper right corner of blank sky region: (place cursor, press x, then q)'
                 ursky = iraf.images.tv.imexamine(input = writeprefix + target, Stdout = 1)
                 llskybreak = split_strings(llsky[0])
                 startx = llskybreak[0]
                 starty = llskybreak[1]
                 urskybreak = split_strings(ursky[0])
                 endx = urskybreak[0]
                 endy = urskybreak[1]
                 skyxs = skyxs+[startx, endx]
                 skyys = skyys+[starty, endy]

                 print 'define lower left corner of sky line region: (place cursor, press x, then q)'
                 llline = iraf.images.tv.imexamine(input = writeprefix + target, Stdout = 1)
                 print 'define upper right corner of sky line region: (place cursor, press x, then q)'
                 urline = iraf.images.tv.imexamine(input = writeprefix + target, Stdout = 1)
                 lllinebreak = split_strings(llline[0])
                 startx = lllinebreak[0]
                 starty = lllinebreak[1]
                 urlinebreak = split_strings(urline[0])
                 endx = urlinebreak[0]
                 endy = urlinebreak[1]
                 linexs = linexs+[startx, endx]
                 lineys = lineys+[starty, endy]

                 i = i + 1

# calculate sky levels (and sigmas) for each sky line region 

        skylevels = [0]
        skysigmas = [0]
        linelevels = [0]
        linesigmas = [0]
        
        j = 0
        while j <= nregions - 1:
            number1 = str(skyxs[2*j+1:2*j+2])
            number1 = number1[2:-5]
            number2 = str(skyxs[2*j+2:2*j+3])
            number2 = number2[2:-5]
            number3 = str(skyys[2*j+1:2*j+2])
            number3 = number3[2:-5]
            number4 = str(skyys[2*j+2:2*j+3])
            number4 = number4[2:-5]
            
            skytestregion = '['+number1+':'+number2+','+number3+':'+number4+']'

#	    print skytestregion

            skystats = iraf.imstat(images = writeprefix + target + skytestregion, fields = 'mean,stddev', Stdout=1) 

#	    print skystats

            splitup = split_strings(skystats[1])
            skylevels = skylevels+[splitup[0]]
            skysigmas = skysigmas+[splitup[1]]
            
# also find the means and sigmas for the sky line regions

            number5 = str(linexs[2*j+1:2*j+2])
            number5 = number5[2:-5]
            number6 = str(linexs[2*j+2:2*j+3])
            number6 = number6[2:-5]
            number7 = str(lineys[2*j+1:2*j+2])
            number7 = number7[2:-5]
            number8 = str(lineys[2*j+2:2*j+3])
            number8 = number8[2:-5]
            
            linetestregion = '['+number5+':'+number6+','+number7+':'+number8+']'

#	    print linetestregion

            linestats = iraf.imstat(images = writeprefix + target + linetestregion, fields = 'mean,stddev', Stdout=1) 

#	    print linestats

            breakup = split_strings(linestats[1])
            linelevels = linelevels+[breakup[0]]
            linesigmas = linesigmas+[breakup[1]]

            j = j + 1

# evaluate how well the selected sky lines have been subtracted

        k = 1
        totallinequality = 0
        while k <= nregions:
            thislinequality = ( float(linelevels[k]) - float(skylevels[k]) ) / float(skysigmas[k])
            totallinequality = totallinequality + thislinequality
            k = k + 1

# adjust the scale and step accordingly
        if totallinequality > 0:
            if scaleup == 0:
                step = step / 2.0
            scaleup = 1
            scale = scale + step
        elif totallinequality <= 0:
            if scaleup == 1:
                step = step / 2.0
            scaleup = 0
            scale = scale - step
            
        
        done = raw_input('Are we done yet? (y or n) ') 

    donestring = 'done with '+target        
    print donestring

    iraf.flprcache()
