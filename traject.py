#!/usr/bin/python

# Analyze simulation trajectories
# Reads in files in PDB format (either as list of files or all in 1 file)
# Calculated the cross-correlation matrix, the covariance matrix, etc.
# Makes plots

from plotting import plot_matrix
import pylab, sys, os.path
from optparse import OptionParser
import PDBlite
import numpy as num

# split a symmetric matrix: positive values in the upper right triangle, negative values
# in the lower left triangle
def split_symmetric_matrix(m):
    m2 = m.copy()
    nrow, ncol = m2.shape

    for i in range(nrow):
        for j in range(ncol):
            if i < j: m2[i,j] = min(m2[i,j], 0) # upper right triangle
            elif i > j: m2[i,j] = max(m2[i,j], 0) # lower left triangle
    m2 = abs(m2)
    return m2


class Trajectory:
    def __init__(self, traj_fn):
        self.traj_fn = traj_fn
        self.name = os.path.basename(traj_fn).replace(".lst","").replace(".pdb","")

    # private: get the x,y,z coords for a residue based on the atom_selection method
    def _get_res_coords(self, res, atom_selection):
        xyz = num.zeros((3,))
        count = 0
        for atom in res.iter_atoms():
            include = False
            if atom_selection == "heavy" and atom.elem != "H": include = True
            elif atom_selection == "CA" and atom.atomName == "CA": include = True

            if include:
                xyz += [atom.x, atom.y, atom.z]
                count += 1
        return xyz/count

    # private: get the x,y,z coords for all residues in a pdb (dimensions: nres x 3) based on the atom_selection method
    def _get_coord_vect(self, pdb, atom_selection, res_range=None):
        nres = pdb.len()
        if res_range != None: nres = res_range[1] - res_range[0] + 1
        
        pdb_coords = num.zeros((nres, 3))
        res_count = 0
        for res in pdb.iter_residues():
            res_num = int(res.res_num)
            if res_range and (res_num < res_range[0] or res_num > res_range[1]): continue
            
            xyz = self._get_res_coords(res, atom_selection)
            #print coord_mat[res_num*3:res_num*3+3, pdb_num].shape, pos.shape
            pdb_coords[res_count, :] = xyz
            res_count += 1
        return pdb_coords

    # returns the nres x npdb x 3 array of coordinates
    # res_range: [start_res, end_res] to process
    def load_coords(self, atom_selection, res_range=None):
        pdbtraj = PDBlite.PDBTrajectory(self.traj_fn)
        self.npdb = pdbtraj.parse_len()

        print "Expecting '%d' pdbs in trajectory" % self.npdb
        pdb_num = 0
        for pdb in pdbtraj.get_next_pdb():
            #print "Loading coords for: ", pdb
            if pdb_num == 0:
                self.nres = pdb.len()
                self.all_coords = num.zeros((self.nres, self.npdb, 3))

            if self.nres != pdb.len():
                print "ERROR: structures in the trajectory don't have the same number of residues (expecting '%d')" % nres
                sys.exit(1)

            pdb_coords = self._get_coord_vect(pdb, atom_selection, res_range)
            self.all_coords[:, pdb_num, :] = pdb_coords
            pdb_num += 1
        print "Loaded coords for '%d' pdbs in trajectory" % pdb_num

        print "Coordinate matrix shape: " + str(self.all_coords.shape)
        print "Loaded %d pdbs with %d residues from trajectory '%s'" % (self.npdb, self.nres, self.traj_fn)

    # Returns the covariance & cross correlation matrices.
    # reference_pdb_num is the index (0-based) of the pdb to use as the reference,
    # or None to use the mean residue positions.
    def get_cov_matrix(self, reference_pdb_num):
        nres, npdb = self.nres, self.npdb

        if reference_pdb_num != None and reference_pdb_num > 0:
            print "%s: Using model '%d' as the reference" % (self.name, reference_pdb_num)
            ref_pdb_coords = self.all_coords[:,reference_pdb_num,:]
        else:
            print "%s: Using average model coordinates as the reference" % self.name
            ref_pdb_coords = self.all_coords.sum(axis=1)/npdb # sum along the pdb axis to create an (nres, 3) shape array

        all_coords2 = self.all_coords.copy()
        for pdb_num in range(npdb):
            all_coords2[:,pdb_num,:] -= ref_pdb_coords

        # calculate the MSD for each residue
        dist2 = num.zeros((nres,))
        for pdb_num in range(npdb):
            for i in range(nres):
                dist2[i] += num.dot(all_coords2[i,pdb_num, :], all_coords2[i,pdb_num, :])
        dist2 /= npdb

        #print "Res 5: x: ",  self.all_coords[4, :, 0]
        #print "Res 5: x: ",  all_coords2[4, :, 0]
            
        # calculate the covariance:
        #    Cij = sum over structs(Ri . Rj)/num_structures
        #       (where Ri is the coordinate vector for residue i)
        mcov = num.zeros((nres, nres))
        mcc = num.zeros((nres, nres))
        for i in range(nres):
            coords_i = all_coords2[i,:,:]
            for j in range(nres):
                coords_j = all_coords2[j,:,:]
                mcov[i,j] = num.sum(coords_i * coords_j)
                #mcov[i,j] += num.dot(self.all_coords[i, pdb_num, :], self.all_coords[j, pdb_num, :])
        mcov /= npdb
        mrmsd = num.sqrt(num.outer(dist2, dist2))
        mcc = mcov / mrmsd
        
        return mcov, mcc, mrmsd
        
# Parse the input arguments
usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-n", "--name", default=None, type="string", help="name to be used for output files")
parser.add_option("-m", "--matrix_fn", default=None, type="string", help="filename to load CC matrix from")
parser.add_option("-t", "--traj_fns", default=None, type="string", help="filenames to load PDB trajectories from in colon separated format (e.g. traj1:traj2)")
parser.add_option("-r", "--res_range", default=None, type="string", help="range of residues to use (e.g. 2:72)")
parser.add_option("-a", "--atom_selection", default="CA", type="string", help="which atoms to use in calculating the residue centroid: CA|heavy")
parser.add_option("-s", "--split_plots", action="store_true", default=False, help="split plots into + and - triangles")
parser.add_option("-f", "--ref_pdb_num", default=None, type="int", help="the pdb num in the trajectory to use as the reference; if 'None' then the averaged residue positions are used")
parser.add_option("-c", "--calc_cc", default=False, action="store_true", help="calculate the cross correlation and covariance matrices for the trajector(y|ies)")
parser.add_option("-d", "--dssp_fn", type="string", default=None, help="filename containing DSSP output")
(opts, args) = parser.parse_args()

if len(args) != 0: parser.error("Incorrect number of arguments")
if opts.matrix_fn == None and opts.traj_fns == None:
    parser.error("Must specify a trajectory or a matrix file to load")

# Init the arg traj_fns
traj_fns = None
if opts.traj_fns != None: traj_fns = filter(lambda fn: fn!="", opts.traj_fns.split(":"))

# Verify the arg res_range
res_range = None
if opts.res_range != None:
    try:
        res_range = map(int, opts.res_range.split(":"))
        if len(res_range) != 2: raise Exception()
    except: parser.error("Invalid syntax for --res_range")

# Verify the atom_selection arg
if opts.atom_selection not in ("CA", "heavy"):
     parser.error("Invalid value for --atom_selection")
print "Using '%s' atoms to specify residue coordinates" % opts.atom_selection

# init the output filename prefix
if opts.name != None: out_fn_prefix = opts.name
elif traj_fns != None:
    if len(traj_fns) == 1: out_fn_prefix = os.path.basename(traj_fns[0])
    else: out_fn_prefix = os.path.basename(traj_fns[0]) + "_MULTTRAJ"
else: out_fn_prefix = "tmp"


### DO STUFF 

#if opts.matrix_fn:
#    mat = pylab.load(opts.matrix_fn)
#    plot_matrix(mat, out_fn_prefix + ".png", opts.split_plots)

# Load trajectories and their coordinates
if traj_fns != None:
    trajs = [Trajectory(traj_fn) for traj_fn in traj_fns]
    num_trajs = len(trajs)
    for traj_num, traj in zip(range(num_trajs), trajs):
        traj.load_coords(opts.atom_selection, res_range)

if opts.calc_cc:
    mcovs, mccs, mrmsds = {}, {}, {}
    for traj_num, traj in zip(range(num_trajs), trajs):
        mcovs[traj_num], mccs[traj_num], mrmsds[traj_num] = traj.get_cov_matrix(None)

    # if only 1 trajectory given, use it; otherwise take the average of the matrices
    if len(trajs) == 1:
        mcov, mcc, mrmsd = mcovs[0], mccs[0], mrmsds[0]
    else:
        mcov, mcc, mrmsd = mcovs[0].copy(), mccs[0].copy(), mrmsds[0].copy()
        for traj_num in range(1,num_trajs):
            mcov, mcc, mrmsd = mcov+mcovs[traj_num], mcc+mccs[traj_num], mrmsd+mrmsds[traj_num]
        mcov, mcc, mrmsd = mcov/num_trajs, mcc/num_trajs, mrmsd/num_trajs

    mcc_threshold = mcc.copy()
    mcc_threshold[mrmsd<.2] = 0

    # get info about SS regions
    ss_info = None
    if opts.dssp_fn != None: ss_info = PDBlite.SS_info(opts.dssp_fn)
        
    plot_matrix(mcov, out_fn_prefix + "_COV.png", split_symmetric=opts.split_plots)
    plot_matrix(mcc, out_fn_prefix + "_CC.png", range=[-1,1], split_symmetric=opts.split_plots)
    plot_matrix(mcc_threshold, out_fn_prefix + "_CC_THRESH.2.png", range=[-1,1], split_symmetric=opts.split_plots, ss_info=ss_info)
    plot_matrix(mcc, out_fn_prefix + "_CC_contour.png", contour_levels=[-.8, -.65, -.5, -.35, -.2, .2, .35, .5, .65, .8],
                split_symmetric=opts.split_plots, ss_info=ss_info)
    plot_matrix(mrmsd, out_fn_prefix + "_RMSD.png", split_symmetric=opts.split_plots)

