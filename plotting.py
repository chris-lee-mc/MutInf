#!/usr/bin/python

# Functions for plotting
# import this module before pylab to enable non-interactive figure creation.

import matplotlib
matplotlib.use('Agg') # so plots aren't made interactively
import pylab, optparse
import numpy as num
import PDBlite

# make a plot of S2s over sequence
# S2s_list and label_list are lists of S2s & labels to be plotted on the same figure
# each S2s element is an numarray of floats
def plot_S2s_over_sequence(S2s_list, label_list, plot_fn, ss_info=None, legend=False, errors_list=None):
    if errors_list != None:
        for S2s, errors, label in zip(S2s_list, errors_list, label_list):
            print range(1, S2s.shape[0]), S2s[1:]
            pylab.errorbar(range(S2s.shape[0]), S2s, fmt="b-", label=label, yerr=errors)
    else:
        for S2s, label in zip(S2s_list, label_list):
            pylab.plot(S2s, "-", label=label)

    # add faded background in SS regions
    if ss_info != None:
        for res_num in ss_info.get_res_nums():
            if ss_info.is_structured(res_num):
                #print "Found SS:", res_num
                x=res_num
                pylab.fill([x-.5,x-.5,x+.5,x+.5], [0,1,1,0], alpha=.3, edgecolor='w')

    #pylab.title("NH order parameters")
    #pylab.ylabel("Order parameter")
    pylab.xlabel("Residue number")
    pylab.ylim(ymax=1)
    pylab.grid()
    if legend: pylab.legend(label_list, prop=matplotlib.font_manager.FontProperties(size='6'), loc='lower right')
    print "Writing ", plot_fn
    pylab.savefig(plot_fn)
    

# reads in s2s_fn and returns an array of s2s over sequence
def load_s2s_file(s2s_fn):
    s2s = {}

    # read in s2s file
    for line in open(s2s_fn):
        if line.strip() == "": continue
        res_num, s2 = line.split()
        s2s[int(res_num)] = float(s2)

    # convert to an array
    s2s_array = num.zeros((1+max(s2s.keys()))) + num.nan
    for res_num in sorted(s2s.keys()):
        s2s_array[res_num] = s2s[res_num]

    return s2s_array




# plot a matrix with a colorbar
#   countour_levels: optional list of numbers to use as filled contour levels
#   range: optional range ([min, max]) of z values to enforce
#   split_symmetric: plot positive values in one triangle and negative values in the other
#   ss_info: optional PDBlite.SS_info object; causes boxes to be created around secondary structure regions
def plot_matrix(m, out_fn, contour_levels=None, range=None, split_symmetric=False, ss_info=None):
    if split_symmetric: m2 = split_symmetric_matrix(m)
    else: m2 = m

    pylab.clf()
    if contour_levels == None:
        if range != None:
            vmin, vmax = range
            pylab.imshow(m2, origin='lower', vmin=vmin, vmax=vmax, interpolation=None)
        else:
            pylab.imshow(m2, origin='lower', interpolation=None)
            
        pylab.colorbar()
        pylab.title(out_fn)
    else:
        pylab.contourf(m2, contour_levels, origin='lower', interpolation=None)
        pylab.colorbar()
        pylab.title(out_fn + "; contours=%s" % " ".join(map(str, contour_levels)), fontsize=8)
    print "Writing matrix to '%s'" % out_fn

    # add boxes around SS regions
    if ss_info != None:
        # create a list of ss_segments; intent is to use these to create less fill objects for speed or to make outlined boxes
        #ss_segments = []
        #curr_ss_segment = None
        #for res_num in ss_info.get_res_nums():
        #    is_structured = ss_info.is_structured(res_num)
        #    if curr_ss_segment == None and is_structured : # initiate segment
        #        curr_ss_segment = [res_num]
        #        ss_segments.append(curr_ss_segment)
        #    elif curr_ss_segment != None and not is_structured: # end segment
        #        curr_ss_segment = None
        #    elif curr_ss_segment != None and is_structured: # add to segment
        #        curr_ss_segment += [res_num]
        #print ss_segments

            
        for res1 in ss_info.get_res_nums():
            for res2 in ss_info.get_res_nums():
                if ss_info.is_structured(res1) and ss_info.is_structured(res2):
                    pylab.fill([res1-.5,res1-.5,res1+.5,res1+.5], [res2-.5,res2+.5,res2+.5,res2-.5], alpha=1, linewidth=.5, fill=False)
    pylab.savefig(out_fn)


if __name__ == '__main__':
    # Parse the input arguments
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(usage)
    parser.add_option("-s", "--s2s_fns", default=None, type="string", help="filenames to load S2s from; colon separated (e.g. fn1:fn2); each file has a 'res_num S2val' on each line")
    parser.add_option("-d", "--dssp_fn", type="string", default=None, help="filename containing DSSP output")
    (opts, args) = parser.parse_args()
    if opts.s2s_fns == None: parser.error("ERROR --s2s_fns option required")

    s2s_fns = filter(lambda fn: fn.strip()!="", opts.s2s_fns.split(":"))

    # Load SS info
    ss_info = None
    if opts.dssp_fn != None: ss_info = PDBlite.SS_info(opts.dssp_fn)

    # Load S2 data
    if opts.s2s_fns != None:
        s2s_list, names = [], []
        for s2s_fn in s2s_fns:
            s2s = load_s2s_file(s2s_fn)
            s2s_list.append(s2s)
            names.append(s2s_fn)

        plot_S2s_over_sequence(s2s_list, names, "S2s.png", ss_info, False)
