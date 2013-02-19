# rlc color_b.py version 6.0

import colorsys,sys
from pymol import cmd
import math

"""
	
AUTHOR 

	Robert L. Campbell
        http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/

        with modifications for 'bw' and 'wr' gradients by Christopher L. McClendon

USAGE

	color_b(selection='sel',
	        gradient='bgr' or 'rgb' or 'bwr' or 'rwb' or 'bmr' or 'rmb' or 'gray' 
					or 'reversegray'
	        mode='hist' or 'ramp', 
					nbins=11, sat=1.0, value=1.0)

    This function allows coloring of a selection as a function of
    B-value, following a gradient of colours.  The gradients can be:

		'bgr': blue -> green   -> red
		'rgb': red  -> green   -> blue
		'bwr': blue -> white   -> red
		'rwb': red  -> white   -> blue
		'bmr': blue -> magenta -> red
		'rmb': red  -> magenta -> blue

		('rainbow' and 'reverserainbow' can be used as synonyms for 
		'bgr' and 'rgb' respectively and 'grey' can be used as a synonym for 'gray').

    The division of B-value ranges can in either of two modes: 'hist' or
    'ramp'. 'hist' is like a histogram (equal-sized B-value increments
    leading to unequal numbers of atoms in each bin). 'ramp' as a ramp
    of B-value ranges with the ranges chosen to provide an equal number
    of atoms in each group.

    You can also specify the saturation and value (i.e. the "s" and "v"
    in the "HSV" color scheme) to be used for the gradient. The defaults
    are 1.0 for both "sat" and "value".

		In the case of the gray scale gradients, "sat" sets the minimum intensity 
		(normally black) and "value" sets the maximum (normally white)

	usage:
	  from within PyMOL do "run color_b.py" to load the function definition.  
		Then you can use for example:

		    color_b (c;a or c;b),mode=ramp,gradient=bwr,nbins=30,sat=.5, value=1.

		to color chains A and B with the Blue-White-Red gradient in 30 colors of equal 
		numbers of atoms in each color.
"""

def whitebluegradient(n,nbins):
    nbins *= 1.0
    if float(n/(nbins)) <   0: n = 0
    if float(n/(nbins)) > 1.0:  n = nbins
    r = 1.0-n/(nbins)
    g = r
    b = 1.0 
    return (r,g,b)

def bluewhitegradient(n,nbins):
    nbins *= 1.0
    if float(n/(nbins)) <   0: n = 0
    if float(n/(nbins)) > 1.0:  n = nbins
    r = 1.0-n/nbins
    g = r
    b = 1.0
    return (1.0 - r, 1.0 - g, b)

def whiteredgradient(n,nbins):
    nbins *= 1.0
    if float(n/(nbins)) <   0: n = 0
    if float(n/(nbins)) > 1.0:  n = nbins
    r = 1.0
    g = 1.0-n/nbins
    b = g
    
    return (r,g,b)

# function for creating the gradient
def make_gradient(sel,gradient,nbins,sat,value):
	if gradient == 'bgr' or gradient == 'rainbow':
		col=[]
		coldesc=[]
		for j in range(nbins):
			# must append the str(sel[j]) to the color name so that it is unique 
			# for the selection
			coldesc.append('col' + str(sel[j]) + str(j))

			# create colors using hsv scale (fractional) starting at blue(.6666667) 
			# through red(0.00000) in intervals of .6666667/(nbins -1) (the "nbins-1" 
			# ensures that the last color is, in fact, red (0)
			# rewrote this to use the colorsys module to convert hsv to rgb
			hsv = (colorsys.TWO_THIRD - colorsys.TWO_THIRD * float(j) / (nbins-1), sat, value)
			#convert to rgb and append to color list
			rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

			col.append(rgb)
			cmd.set_color("col" + str(sel[j]) + str(j),col[j])

	elif gradient == 'rgb' or gradient == 'reverserainbow':
		col=[]
		coldesc=[]
		for j in range(nbins):
			# must append the str(sel[j]) to the color name so that it is unique 
			# for the selection
			coldesc.append('col' + str(sel[j]) + str(j))

			# create colors using hsv scale (fractional) starting at red(.00000) 
			# through blue(0.66667) in intervals of .6666667/(nbins -1) (the "nbins-1" 
			# ensures that the last color is, in fact, red (0)
			# rewrote this to use the colorsys module to convert hsv to rgb
			hsv = (colorsys.TWO_THIRD * float(j) / (nbins-1), sat, value)
			#convert to rgb and append to color list
			rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

			col.append(rgb)
			cmd.set_color("col" + str(sel[j]) + str(j),col[j])

	elif gradient == 'bmr':
		col=[]
		coldesc=[]
		for j in range(nbins):
			coldesc.append('col' + str(sel[j]) + str(j))
			# create colors in a gradient from blue through magenta to red
			rgb = [min(1.0, float(j)*2/(nbins-1)), 0.0, min(1.0, float(nbins-j-1)*2/(nbins-1))]

			# convert rgb to hsv,  modify saturation and value and convert back to rgb
			hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
			hsv[1] = hsv[1]*sat
			hsv[2] = hsv[2]*value
			rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

			col.append(rgb)
			cmd.set_color("col" + str(sel[j]) + str(j),col[j])

	elif gradient == 'rmb':
		col=[]
		coldesc=[]
		for j in range(nbins):
			coldesc.append('col' + str(sel[j]) + str(j))
			# create colors in a gradient from blue through magenta to red
			rgb = [min(1.0, float(nbins-j-1)*2/(nbins-1)), 0.0, min(1.0, float(j)*2/(nbins-1))]

			# convert rgb to hsv,  modify saturation and value and convert back to rgb
			hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
			hsv[1] = hsv[1]*sat
			hsv[2] = hsv[2]*value
			rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

			col.append(rgb)
			cmd.set_color("col" + str(sel[j]) + str(j),col[j])

	elif gradient == 'bwr':
		col=[]
		coldesc=[]
		for j in range(nbins/2):
			coldesc.append('col' + str(sel[j]) + str(j))
			# create colors in a gradient from blue to white 
			rgb = [min(1.0, float(j)*2/(nbins-1)), min(1.0,float(j)*2/(nbins-1)), min(1.0, float(nbins-j-1)*2/(nbins-1))]

			# convert rgb to hsv,  modify saturation and value and convert back to rgb
			hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
			hsv[1] = hsv[1]*sat
			hsv[2] = hsv[2]*value
			rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

			col.append(rgb)
			cmd.set_color("col" + str(sel[j]) + str(j),col[j])

		for j in range(nbins/2,nbins):
			coldesc.append('col' + str(sel[j]) + str(j))
			# create colors in a gradient from white to red
			rgb = [min(1.0, float(j)*2/(nbins-1)), min(1.0,float(nbins-j-1)*2/(nbins-1)), min(1.0, float(nbins-j-1)*2/(nbins-1))]

			# convert rgb to hsv,  modify saturation and value and convert back to rgb
			hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
			hsv[1] = hsv[1]*sat
			hsv[2] = hsv[2]*value
			rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

			col.append(rgb)
			cmd.set_color("col" + str(sel[j]) + str(j),col[j])

	elif gradient == 'bw':
		col=[]
		coldesc=[]
		for j in range(nbins):
			coldesc.append('col' + str(sel[j]) + str(j))
			# create colors in a gradient from blue to white 
			rgb = [min(1.0, float(j)/(nbins-1)), min(1.0,float(j)/(nbins-1)), min(1.0, float(nbins - j -1)/(nbins-1))]
			rgb = bluewhitegradient(j,nbins)
			print rgb
			# convert rgb to hsv,  modify saturation and value and convert back to rgb
			hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
			hsv[1] = hsv[1]*sat
			hsv[2] = hsv[2]*value
			rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

			col.append(rgb)
			cmd.set_color("col" + str(sel[j]) + str(j),col[j])

	elif gradient == 'wr':
		col=[]
		coldesc=[]
		for j in range(nbins):
			coldesc.append('col' + str(sel[j]+str(nbins)) + str(j))
			# create colors in a gradient from white to red
			#value = nbins/2 + float(j/2.0) #nbins/2 + float(j/2.0)
			#value = j
			#rgb = [min(1.0, float(value)*2/(nbins-1)), min(1.0, float(nbins -value-1)*2/(nbins-1)), min(1.0, float(nbins -value -1)*2/(nbins-1))]
			rgb = whiteredgradient(j,nbins) #[float(j)/(nbins-1), 1.0, 1.0] # float(nbins-j+1)/(nbins-1), 1.0]#float(nbins-j+1)/(nbins-1)]
			#rgb = [min(1.0, float(j)/(nbins-1)), min(1.0,float(nbins-j-1)/(nbins-1)), min(1.0, float(nbins-j-1)/(nbins-1))]

			# convert rgb to hsv,  modify saturation and value and convert back to rgb
			hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
			hsv[1] = hsv[1]*sat
			hsv[2] = hsv[2]*value
			rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

			col.append(rgb)
			cmd.set_color("col" + str(sel[j]+str(nbins)) + str(j),col[j])

	elif gradient == 'rwb':
		col=[]
		coldesc=[]
		for j in range(nbins/2):
			coldesc.append('col' + str(sel[j]) + str(j))
			# create colors in a gradient from red to white
			rgb = [min(1.0, float(nbins-j-1)*2/(nbins-1)), min(1.0,float(j)*2/(nbins-1)), min(1.0, float(j)*2/(nbins-1))]

			# convert rgb to hsv,  modify saturation and value and convert back to rgb
			hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
			hsv[1] = hsv[1]*sat
			hsv[2] = hsv[2]*value
			rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

			col.append(rgb)
			cmd.set_color("col" + str(sel[j]) + str(j),col[j])

		for j in range(nbins/2,nbins):
			coldesc.append('col' + str(sel[j]) + str(j))
			# create colors in a gradient from white to blue
			rgb = [min(1.0, float(nbins-j-1)*2/(nbins-1)), min(1.0,float(nbins-j-1)*2/(nbins-1)), min(1.0, float(j)*2/(nbins-1))]

			# convert rgb to hsv,  modify saturation and value and convert back to rgb
			hsv = list(colorsys.rgb_to_hsv(rgb[0],rgb[1],rgb[2]))
			hsv[1] = hsv[1]*sat
			hsv[2] = hsv[2]*value
			rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

			col.append(rgb)
			cmd.set_color("col" + str(sel[j]) + str(j),col[j])

	elif gradient == 'gray' or gradient == 'grey':
		col=[]
		coldesc=[]
		for j in range(nbins):
			coldesc.append('col' + str(sel[j]) + str(j))
			# create colors in a gradient of grays from "sat" to "value"

			hsv = [0, 0, sat + (value-sat)*float(j)/(nbins-1)]
#			hsv[1] = hsv[1]*sat
#			hsv[2] = hsv[2]*value
			rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

			col.append(rgb)
			cmd.set_color("col" + str(sel[j]) + str(j),col[j])

	elif gradient == 'reversegray' or gradient == 'reversegrey':
		col=[]
		coldesc=[]
		for j in range(nbins):
			coldesc.append('col' + str(sel[j]) + str(j))
			# create colors in a gradient of grays from "sat" to "value"

			hsv = [0, 0, value - (value-sat)*float(j)/(nbins-1)]
#			hsv[1] = hsv[1]*sat
#			hsv[2] = hsv[2]*value
			rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])

			col.append(rgb)
			cmd.set_color("col" + str(sel[j]) + str(j),col[j])


# return the gradient as a color description that include the selection as 
# part of its name
	return coldesc

# main function called from within PyMOL
def color_b(selection="all",mode="hist",gradient="bgr",nbins=11,sat=1.,value=1.):

	nbins=int(nbins)
	sat=float(sat)
	value=float(value)

# make sure lowercase
	gradient.lower()
	mode.lower()

# Sanity checking
	if nbins == 1:
		print "\n     WARNING: You specified nbins=1, which doesn't make sense...resetting nbins=11\n"
		nbins=11

	if mode not in ('hist','ramp'):
		print "\n     WARNING: Unknown mode ",mode, "    ----->   Nothing done.\n"
		return
	elif gradient not in ('bgr','rgb','rainbow','reverserainbow','bwr','bw','wr','rwb',
	                      'bmr','rmb','gray','grey','reversegray','reversegrey'):
		print "\n     WARNING: Unknown gradient: ",gradient, "    ----->   Nothing done.\n"
		return

	print "MODE, GRADIENT, NBINS:", mode,gradient, nbins

# get list of B-factors from selection
	m = cmd.get_model(selection)
	sel = []
	b_list = []

	if len(m.atom) == 0:
		print "Sorry, no atoms selected"

	else:
		for i in range(len(m.atom)):
			b_list.append(m.atom[i].b)

			# create list of b-limits (bottom limit of bin)
			b_lim=[]

		max_b = max(b_list)
		min_b = min(b_list)
		print "Minimum and Maximum B-values: ", min_b, max_b

		if mode == 'ramp':
			# color in bins of equal numbers of atoms
			b_list.sort()

			# subtract 0.1 from the lowest B in order to ensure that the single
			# atom with the lowest B value doesn't get omitted
			b_list[0] = b_list[0] - 0.1

			bin_num = int(len(b_list)/nbins)
			for j in range(nbins):
				sel.append(selection + " and b > " + str(b_list[j*bin_num]))
				#print "Color select: ",sel[j]

		elif mode == 'hist':
			# histogram:
			# color in bins of equal B-value ranges
			# subtract 0.1 from the lowest B in order to ensure that the single
			# atom with the lowest B value doesn't get omitted
			min_b = min_b - 0.1
			bin_width = (max_b - min_b)/nbins
			for j in range(nbins):
				sel.append(selection + " and b > " + str(min_b + j*bin_width))
				#print "Color select: ",sel[j]

# call the function to create the gradient which returns a list of colours
		colours = make_gradient(sel,gradient,nbins,sat,value)

# do the colouring now
		for j in range(nbins):
			#print "Color select: ",sel[j]
			cmd.color(colours[j],sel[j])

# allow calling without parentheses: color_hist_b [selection=], [mode= ],[gradient= ],[nbins= ]
cmd.extend("color_b",color_b)
