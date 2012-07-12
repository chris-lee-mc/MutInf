#!/usr/bin/python
# Kullback-Leibler Divergence Program, Copyright 2011 Christopher McClendon
# Released under Lesser GNU Public License

from dihedral_mutent import *
from optparse import OptionParser
from utils import arr2str2, fmt_floats
from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser
#from scipy import special

global adaptive_partitioning

def load_two_resfiles(run_params1, run_params2, load_angles=True, all_angle_info=None):
        rp1 = run_params1
        rp2 = run_params2
        if rp1.num_structs == None: rp1.num_structs = 1500000
        if rp2.num_structs == None: rp2.num_structs = 1500000
        sequential_num = 0
        resfile1=open(rp1.resfile_fn,'r')
        reslines1=resfile1.readlines()
        resfile1.close()
        resfile2=open(rp2.resfile_fn,'r')
        reslines2=resfile2.readlines()
        resfile2.close()
        reslist1 = []
        reslist2 = []
        for resline1, resline2 in zip(reslines1, reslines2):
            if (len(resline1.strip()) == 0) or  ((len(resline1.strip()) == 0)): continue
            for iteration in range(2):
                xvg_resnum, res_name, res_numchain = resline1.split()
                myexpr = re.compile(r"([0-9]+)([A-Z]*)")
                matches = myexpr.match(res_numchain)
                res_num = matches.group(1)
                if matches.group(2) != None:
                    res_chain = matches.group(2)
                else:
                    res_chain = " "
                if load_angles: 
                    if(iteration == 0):
                        minmax1 = ResidueChis(res_name,res_num, xvg_resnum, rp1.xvg_basedir, rp1.num_sims, rp1.num_structs, rp1.xvgorpdb, rp1.binwidth, rp1.sigalpha, rp1.permutations, rp1.phipsi, rp1.backbone_only, rp1.adaptive_partitioning, rp1.which_runs, rp1.pair_runs, bootstrap_choose = rp1.bootstrap_choose, calc_variance=rp1.calc_variance, all_angle_info=all_angle_info, xvg_chidir=rp1.xvg_chidir, skip=rp1.skip,skip_over_steps=rp1.skip_over_steps, calc_mutinf_between_sims=rp1.calc_mutinf_between_sims,max_num_chis=rp1.max_num_chis, sequential_res_num = sequential_num, pdbfile=rp1.pdbfile, xtcfile=rp1.xtcfile, output_timeseries=rp1.output_timeseries, bailout_early=True )
                        print "min value: "+str(minmax1.minmax[0])
                    else: #this time with min/max supplied
                        reslist1.append(ResidueChis(res_name,res_num, xvg_resnum, rp1.xvg_basedir, rp1.num_sims, rp1.num_structs, rp1.xvgorpdb, rp1.binwidth, rp1.sigalpha, rp1.permutations, rp1.phipsi, rp1.backbone_only, rp1.adaptive_partitioning, rp1.which_runs, rp1.pair_runs, bootstrap_choose = rp1.bootstrap_choose, calc_variance=rp1.calc_variance, all_angle_info=all_angle_info, xvg_chidir=rp1.xvg_chidir, skip=rp1.skip,skip_over_steps=rp1.skip_over_steps, calc_mutinf_between_sims=rp1.calc_mutinf_between_sims,max_num_chis=rp1.max_num_chis, sequential_res_num = sequential_num, pdbfile=rp1.pdbfile, xtcfile=rp1.xtcfile, output_timeseries=rp1.output_timeseries, minmax=minmax1.minmax ))
                        del minmax1
                if (load_angles != True): 
                    if (iteration == 1):
                        reslist.append(ResListEntry(res_name,res_num,res_chain))
                
                xvg_resnum, res_name, res_numchain = resline2.split()
                myexpr = re.compile(r"([0-9]+)([A-Z]*)")
                matches = myexpr.match(res_numchain)
                res_num = matches.group(1)
                if matches.group(2) != None:
                    res_chain = matches.group(2)
                else:
                    res_chain = " "
                if load_angles: 
                    if(iteration == 0):
                        minmax2 = ResidueChis(res_name,res_num, xvg_resnum, rp2.xvg_basedir, rp2.num_sims, rp2.num_structs, rp2.xvgorpdb, rp2.binwidth, rp2.sigalpha, rp2.permutations, rp2.phipsi, rp2.backbone_only, rp2.adaptive_partitioning, rp2.which_runs, rp2.pair_runs, bootstrap_choose = rp2.bootstrap_choose, calc_variance=rp2.calc_variance, all_angle_info=all_angle_info, xvg_chidir=rp2.xvg_chidir, skip=rp2.skip,skip_over_steps=rp2.skip_over_steps, calc_mutinf_between_sims=rp2.calc_mutinf_between_sims,max_num_chis=rp2.max_num_chis, sequential_res_num = sequential_num, pdbfile=rp2.pdbfile, xtcfile=rp2.xtcfile, output_timeseries=rp2.output_timeseries, bailout_early = True)
                        print "min value: "+str(minmax2.minmax[0])
                    else: #this time with min/max supplied
                        reslist2.append(ResidueChis(res_name,res_num, xvg_resnum, rp2.xvg_basedir, rp2.num_sims, rp2.num_structs, rp2.xvgorpdb, rp2.binwidth, rp2.sigalpha, rp2.permutations, rp2.phipsi, rp2.backbone_only, rp2.adaptive_partitioning, rp2.which_runs, rp2.pair_runs, bootstrap_choose = rp2.bootstrap_choose, calc_variance=rp2.calc_variance, all_angle_info=all_angle_info, xvg_chidir=rp2.xvg_chidir, skip=rp2.skip,skip_over_steps=rp2.skip_over_steps, calc_mutinf_between_sims=rp2.calc_mutinf_between_sims,max_num_chis=rp2.max_num_chis, sequential_res_num = sequential_num, pdbfile=rp2.pdbfile, xtcfile=rp2.xtcfile, output_timeseries=rp2.output_timeseries, minmax=minmax2.minmax ))
                        del minmax2
                if (load_angles != True): 
                    if (iteration == 1):
                        reslist.append(ResListEntry(res_name,res_num,res_chain))
            sequential_num += 1 
        return reslist1, reslist2


def Grassberger_KLdiv(nj, ni, numangles1, numangles2):
    ## Here, nj are the reference counts
    # This is based on the Renyi generalized divergence using Grassberger's (1988 Phys Lett A) estimate of pi and pj
    # the following is written in LaTeX notation:
    # D_{\alpha} (P || Q) = 1/(\alpha-1) * log (sum_{i=1}^{n} \frac{p_i^{\alpha}}{q_i^{\alpha-1}}
    # Kullback-Leibler Divergence (\alpha=1): lim(\alpha->1) D_\alpha (P || Q)
    
    #To obtain a finite-sample size correction to the Kullback-Leibler Divergence, we will follow the approach of (Grassberger 1988).
    #We will consider the Kullback-Leibler Divergence as a limit of the Renyi Generalized Divegence,

    # \begin{equation}
    #  KL_{n} = lim(\alpha->1) D_{\alpha} (p_{i} || p_{i^*})
    # \end{equation}

    # where
    # \begin{equation}
    # D_{\alpha}((p_{i} || p_{i^*}) = 1/(\alpha-1) * \ln (\sum_{i=1}^{n} \frac{p_{i}^{\alpha}}{p_{i^*}^{\alpha-1}})
    # \end{equation}

    # For finite sample sizes there will be some uncertainty in the $p_i$. Considering the actual histogram counts, we write:
    # \begin{equation}
    # p_{i}^{\alpha} = (<n_{i}>/N)^{\alpha},\, p_{i^*}^{\alpha} = (<n_{i^*}>/N)^{\alpha}
    # \end{equation}

    # To obtain $<n>^{\alpha}$, we assume a Poisson distribution for $n_{i}$ in successive realizations (i.e. assuming we are using a fine enough discretization such that $p_{i}$<<1).
    # For integer $\alpha$, we would then have
    # \begin{equation}
    # <n>^{\alpha} = \Big<\frac{n!}{(n-\alpha)!}\Big>, \alpha>=1
    # \end{equation}

    # However, as we consider the limit as $\alpha$ approaches 1, we need a continuous analog using $\Gamma$ functions.
    # Grassberger found an asymptotic expansion for $<n_{i}>^{\alpha}$ and showed that two terms gave numerically robust results for Shannon etropies. 
    # \begin{equation}
    # <n>^{\alpha} = {\frac {\Gamma   ( n+1  ) }{\Gamma   ( n-a+1  ) }}-{\frac {  ( -1  ) ^{n}\Gamma   ( a+1  ) \sin  ( \pi \,a  ) }{\pi   ( n+1  ) }}
    # \end{equation}

    # This same approximation is used in our previously-published MutInf method (ref).
    # Then, using this expression for $<n>$ Evaluating the Renyi Generalized Divergence in the $\alpha \rightarrow 1$ limit to give us the Kullback-Leibler Divergence,
    # we invoke L'Hopital's Rule to obtain: 
    # \begin{equation}
    # \lim_{\alpha \rightarrow 1} D_{\alpha} ((p_{i} || p_{i^*}) =  
    # \lim _{\alpha\rightarrow 1} \Big( \sum _{i=1}^{{\it nbins}}{\frac {Nf \Big( {\it n_{i^*}},\alpha-1 \Big) }{f \Big( {\it n_i},\alpha \Big) }}
    # \mbox{}{\frac {\partial }{\partial \alpha}}\sum _{i=1}^{{\it nbins}}{\frac {f \Big( {\it n_i},\alpha \Big) }{Nf \Big( {\it n_{i^*}},\alpha-1 \Big) }} \Big) 
    # \end{equation}

    # \begin{equation}
    # \lim_{\alpha \rightarrow 1} D_{\alpha} ((p_{i} || p_{i^*}) = \sum_{i}^{nbins} \frac{\Psi ( {\it n_i} ) {\it n_i}\,{\it n_{i^*}}+ 
    # ( -1 ) ^{{\it n_i}}{\it n_{i^*}}+ ( -1 ) ^{{\it n_{i^*}}}{\it n_i}\,{\it n_{i^*}}-\Psi ( {\it n_{i^*}} ) {\it n_i}\,{\it n_{i^*}}-{\it n_i}}{{\it n_{i^*}}}
    # \end{equation}

    # However, this expression is not numerically robust, so we truncate the expression for $<n>^{alpha}$ at the first term:
    # \begin{equation}
    # <n>^{\alpha} = {\frac {\Gamma   ( n+1  ) }{\Gamma   ( n-a+1  ) }}
    # \end{equation}
    # which then provides a more robust estimate for $D_{\alpha} (P || Q)$:
    # \begin{equation}
    # KL_{n} = \lim_{\alpha \rightarrow 1} D_{\alpha} ((p_{i} || p_{i^*}) = \sum_{i}^{nbins} \frac{n_i}{N}\Big(\Psi ( {\it n_i} ) - \Psi ( {\it n_{i^*}} )-\frac{1}{\it n_{i^*}}\Big)
    # \end{equation}
    # Using a series approximation of the digamma function, $\Psi(x) \approx \ln(x) - \frac{1}{2x}$,
    # it can be readily seen that the regular Kullback-Leibler Divergence is recovered along with a correction term that decreases in size as histogram counts increase.

    psi_ni = (special.digamma(ni + SMALL))[nj > 0]
    psi_nj = (special.digamma(nj + SMALL))[nj > 0]
    n_i = ni[nj > 0]
    n_j = nj[nj > 0]
    sign1 = ((-1) ** int16(n_i % 2))
    sign2 = ((-1) ** int16(n_j % 2))
    return (1 / float64(numangles1)) * sum( n_i * (psi_ni - psi_nj - 1.0 / float64(n_j)))

        ## a reordered version of the following, to work when ni=0 by setting to zero that term above
        ##(psi_ni * n_i * n_j + n_j * sign1  \
        ##  + n_i * n_j * sign2 - psi_nj * ni * nj - n_i) \
        ##    / (n_i * n_j))



def filter_div(div_ref, div_ref2, div1):
    #  If the "target" ensemble is the same as the equilibrium ensemble, this quantity will be zero. However, this is not often the case due to sample variability. Furthermore, if applied naively, it might be difficult to extract meaningful population shifts due to different simulation conditions versus artefactual popualtion shifts due to sample variability. In order to improve the signal-to-noise ratio in our calculation the Kullback-Leibler Divergence and thereby capture meaningful differences between conformational ensembles, we will calculate the K-L Divergence expected from sample variability in the "reference" ensemble and use it for a significance test and to correct the calculated values.
    
    
    #  To generate a realistic measure of sample variability, we use a statistical bootstrapping approach. We split the full reference ensembles into $nsims$ blocks (usually corresponding to clones of the same system with different random number seeds), and take half of the blocks at a time as a surrogate target ensemble and the complementing half as a surrogate reference ensemble. We aggregate the counts for the torsions to construct probability distributions and calculate the K-L divergence between all combinations of surrogate distributions. Any non-zero average K-L divergence between these distributions is a measure of average bias that we can later subtract from the total K-L divergence between the full "reference" ensemble and the full "target" ensemble, when it is significant.  The K-L Divergence under the null hypothesis that the average K-L divergence is no greater than that expected from sample variability in the reference ensemble is then given by:
    
    #  \begin{equation}
    #  KL_{i}^{H_0} = {nsims \choose nsims/2}^{-1} \sum_{blocks}^{\frac{nsims}{2}} \sum_{i}^{nbins}p_{i}\ln \frac{p_{i}^{S}}{p_{i^{S^C}}}
    #  \end{equation}

    #  where $S$ denotes subsamples and $S^C$ are their complements.
    #  To test for statistical significance of the obseved K-L divergence, we use the distribution of these surrogate K-L divergence values to obtain a p-value for the null hypothesis that the average K-L divergence is no greater than that expected from sample variability in the reference ensemble. If this p-value for a particular torsion is less than the significance level (in this case, set at a permissive $\alpha=0.1$), then the K-L divergence is set to zero; if not, then the average K-L divergence between the surrogate distributions described above is subtracted from the total:

    #  \begin{equation}
    #  KL_{i} =  {nsims \choose \frac{nsims}{2}}^{-1} \sum_{blocks}^{\frac{nsims}{2}} \sum_{i}^{nbins}p_{i}\ln \frac{p_i^S}{p_{i^{S^C}}}
    #  \end{equation}

    counter = 0 #residue list counter -- increments along residue list
    for res1, res2, res3 in zip(div1.reslist1, div_ref.reslist1, div_ref2.reslist1):
       print "Residue (referece): "+str(res1.name) + " "+str(res1.num)
       for mychi in range(6):
            print "Dihedral: "+ str(mychi + 1)
            kldiv_target_avg = average(div1.kldiv_res[counter,:,mychi])
            kldiv_refs = div_ref.kldiv_res[counter,:,mychi]
            kldiv_targs = div1.kldiv_res[counter,:,mychi]
            kldiv_null_hyp_ref = average(div_ref.kldiv_res[counter,:,mychi]) #average over bootstraps
            kldiv_num_greater_than_ref = len(kldiv_refs[kldiv_refs > kldiv_target_avg]) #how many of ref sub-ensembles have a KLdiv wrt full ref ensemble greater than the KLdiv of the target
            bootstrap_sets = div1.bootstrap_sets
            jsdiv_target_avg = average(div1.jsdiv_res[counter,:,mychi])
            jsdiv_refs  = div_ref.jsdiv_res[counter,:,mychi]
            jsdiv_refs2 = div_ref2.jsdiv_res[counter,:,mychi]
            jsdiv_targs = div1.jsdiv_res[counter,:,mychi]
            jsdiv_null_hyp_ref = (average(div_ref.jsdiv_res[counter,:,mychi]) + average(div_ref2.jsdiv_res[counter,:,mychi])) / 2.0 #average over bootstraps
            jsdiv_num_greater_than_ref = len(jsdiv_refs[jsdiv_refs > jsdiv_target_avg]) + len(jsdiv_refs[jsdiv_refs2 > jsdiv_target_avg]) #how many of sub-ensembles of reference or sub-ensembles of target  have a JSdiv greater than the avg JSdiv between ref and target
            
            #filter
            print "KLdiv avg total: "+str(kldiv_target_avg)
            print "KLdiv null hyp:  "+str(kldiv_null_hyp_ref)
            print "kldiv null hyp samples > target: "+str(kldiv_num_greater_than_ref)
            print "KLdiv corrected: "+str(div1.kldiv_res[counter,:,mychi]) + " bias: "+str( kldiv_null_hyp_ref)
            if(kldiv_num_greater_than_ref/bootstrap_sets > div1.sigalpha):
                div1.kldiv_res[counter,:,mychi] = 0 #zero all bootstrap_sets if it is filtered out. Hope variance is zero-safe
                print "this torsion's KL divergence was not significant."
            else:
                #remove bias due to kldiv of reference with respect to itself
                print "significant KL divergence: correcting for bias."
                kldiv_targs -= kldiv_null_hyp_ref  
                kldiv_targs[kldiv_targs < 0] = 0
                div1.kldiv_res[counter,:,mychi] = kldiv_targs 
            print 
                
            #filter
            print "JSdiv avg total: "+str(jsdiv_target_avg)
            print "JSdiv null hyp:  "+str(jsdiv_null_hyp_ref)
            print "jsdiv null hyp samples > target: "+str(jsdiv_num_greater_than_ref)
            print "JSdiv corrected: "+str(div1.jsdiv_res[counter,:,mychi]) + " bias: "+str( jsdiv_null_hyp_ref)
            if(jsdiv_num_greater_than_ref/(bootstrap_sets * 2) > div1.sigalpha): # divide by two times bootstrap_sets for since added two tail distributions above for JSDiv within target and within reference sub-ensembles 
                div1.jsdiv_res[counter,:,mychi] = 0 #zero all bootstrap_sets if it is filtered out. Hope variance is zero-safe
                print "this torsion's JS divergence was not significant. and will be set to zero"
            else:
                #remove bias due to kldiv of reference with respect to itself
                print "significant JS divergence."
                jsdiv_targs -= jsdiv_null_hyp_ref  
                jsdiv_targs[jsdiv_targs < 0] = 0
                div1.jsdiv_res[counter,:,mychi] = jsdiv_targs 
            print
       counter += 1 



 #########################################################################################################################################
 ###  Class KLdiv: computes Kullback-Leibler Divergence Expansion first and second-order terms for all residues given two residue lists
 ##########################################################################################################################################

class KLdiv:
    reslist1 = []
    reslist2 = []
    maxfile_kl = []
    sumfile_kl = []
    lastchi_kl = []
    maxfile_chisq = []
    sumfile_chisq = []
    lastchi_kl = []
    maxfile_chisq = []
    sumfil_chisq = []
    kldiv = []
    jsdiv = []
    chisq = []
    bootstrap_sets = 0
    nbins = 0
    nsims = 0
    numangles1 = 0
    numangles2 = 0
    min_numanlges = 0
    nchi1 = 0
    nchi2 = 0
    nchi_to_use = 0
    kldiv_done = 0 #flag, not used for anything yet
    kldiv_cov = []
    chisq_cov = []
    kldiv_corr = []
    chisq_corr = []
    sigalpha = 0
    backbone_only = 0
    blocks_ref = 6
    
    def cov(self):
        #self.kldiv_cov = zeros((len(self.reslist1), len(self.reslist2)), float64)
        #self.chisq_cov = zeros((len(self.reslist1), len(self.reslist2)), float64)
        #self.kldiv_cov = cov(self.kldiv, rowvar=0) #variables in columns (fastest-varying), observations in rows (slowest-varying)
        #self.chisq_cov = cov(self.chisq, rowvar=0) #variables in columns (fastest-varying), observations in rows (slowest-varying)
        #sum over torsions in a res before calculating covariance matrix
        print "kldiv_res sum over chis:"
        print sum(self.kldiv_res,axis=-1)
        print "cov:"
        print array(cov(sum(self.kldiv_res,axis=-1)))
        return [array(cov(sum(self.kldiv_res,axis=-1), rowvar=1)), array(cov(sum(self.chisq_res,axis=-1), rowvar=1)), array(cov(sum(self.jsdiv_res,axis=-1), rowvar=1))] # as observations in bootstrap_sets, which are faster-varying than te res index
    
    def corrcoef(self):
        #self.kldiv_cov = zeros((len(self.reslist1), len(self.reslist2)), float64)
        #self.chisq_cov = zeros((len(self.reslist1), len(self.reslist2)), float64)
        #self.kldiv_cov = cov(self.kldiv, rowvar=0) #variables in columns (fastest-varying), observations in rows (slowest-varying)
        #self.chisq_cov = cov(self.chisq, rowvar=0) #variables in columns (fastest-varying), observations in rows (slowest-varying)
        #sum over torsions in a res before calculating correlation matrix
        return (corrcoef(sum(self.kldiv_res,axis=-1), rowvar=1), cov(sum(self.chisq_res,axis=-1), rowvar=1),  cov(sum(self.jsdiv_res,axis=-1), rowvar=1)) # as observations in bootstrap_sets, which are faster-varying than the res index

    def replace_bfac_res(self, replacements, structure, outpdb):        
        model = structure[0]
        counter1 = 0
        numchain_re = re.compile(r'\s*([0-9]+)([A-Z]*).*')
        numchain_re2 = re.compile(r'\s*([0-9]+).*')
        for res1 in self.reslist1:
            test = None
            
            print "numchain: "+str(res1.num)+ " resname: "+str(res1.name)
            mymatch = numchain_re.match(str(res1.num))
            if mymatch != None:
                (resnum, reschain) = mymatch.groups()
                print "match: "+str(res1.num)+" "+str(resnum)
                if(reschain == ''):
                    reschain = ' '
                test = 1
            else:  
                try:
                    #model must have only one chain or resfile didn't specify chain or can't find match
                    for chain in model:
                        reschain = chain.id
                        print "reschain: "+str(reschain)
                        #test = model[reschain]
                        mymatch = numchain_re2.match(str(res1.num))
                        if mymatch != None:
                            (resnum) = mymatch.groups()
                            print "match: "+str(res1.num)+" "+str(resnum)
                            test = 1
                        else:
                            print "did not match"
                except:
                    print "did not match"
            try:
                if(test != None):
                    print "chain: "+str(reschain)
                    for atom in model[reschain][int(resnum)]:
                        print "atom: "+str(atom.id)
                        atom.bfactor = sum(average(replacements[counter1,:,:],axis=0))
            except:
                for chain in model:
                    print "empty chain name: "+str(chain.id)
                    for atom in chain[int(resnum)]:
                        print "atom: "+str(atom.id)
                        atom.bfactor = sum(average(replacements[counter1,:,:],axis=0))
            counter1 +=1 
            
        w = PDBIO()
        w.set_structure(structure)
        w.save(outpdb)
            
    
    def write_pdbs(self, postfix = None):
        if(postfix == None):
            postfix = ""
        PDB_input = self.run_params1.pdbfile
        print "Reading PDB file: "+str(PDB_input)
        parser = PDBParser()
        prefix=str(self.run_params1.resfile_fn)+str(self.run_params2.resfile_fn)
        kl_pdbfile = prefix + "_kldiv" + postfix + ".pdb"
        js_pdbfile = prefix + "_jsdiv" + postfix + ".pdb"
        try:
            structure_kl = parser.get_structure('self1', PDB_input)
            self.replace_bfac_res(self.kldiv_res, structure_kl, kl_pdbfile)
            structure_js = parser.get_structure('self2', PDB_input)
            self.replace_bfac_res(self.jsdiv_res, structure_js, js_pdbfile)
        except:
            print "Sorry, BioPython didn't like your pdb most likely, you will have to use replace_res_bfac_mut.pl to replace the b-factors in the pdb for this case\n"

    
    def output_cov(self):
        prefix=str(self.run_params1.resfile_fn)+str(self.run_params2.resfile_fn)
        try:
            if(self.kldiv_cov == [] or self.chisq_cov == []):
                (self.kldiv_cov, self.chisq_cov, self.jsdiv_cov) = self.cov()
            if(self.kldiv_corr == [] or self.chisq_corr == []):
                (self.kldiv_corr, self.chisq_corr, self.jsdiv_corr) = self.corrcoef()
        except:
            self.kldiv_cov = zeros((len(self.reslist1), len(self.reslist2)), float64)
            self.jsdiv_cov = zeros((len(self.reslist1), len(self.reslist2)), float64)
            self.chisq_cov = zeros((len(self.reslist1), len(self.reslist2)), float64)
            self.kldiv_corr = zeros((len(self.reslist1), len(self.reslist2)), float64)
            self.jsdiv_corr = zeros((len(self.reslist1), len(self.reslist2)), float64)
            self.chisq_corr = zeros((len(self.reslist1), len(self.reslist2)), float64)
        
        ##WRITE COVARIANCE MATRIX USING ROUTINE FROM dihedral_mutent.py
        name_num_list = make_name_num_list(self.reslist1)
        print self.kldiv_cov
        try:
            output_matrix(prefix+"_bootstrap_kldiv_covar_res.txt",            self.kldiv_cov ,name_num_list,name_num_list)
            output_matrix(prefix+"_bootstrap_kldiv_covar_res_0diag.txt",      self.kldiv_cov ,name_num_list,name_num_list, zero_diag=True)
            output_matrix(prefix+"_bootstrap_kldiv_corr_res.txt",             self.kldiv_corr,name_num_list,name_num_list)
            output_matrix(prefix+"_bootstrap_kldiv_corr_res_0diag.txt",       self.kldiv_corr,name_num_list,name_num_list, zero_diag=True)
            
            output_matrix(prefix+"_bootstrap_jsdiv_covar_res.txt",            self.jsdiv_cov ,name_num_list,name_num_list)
            output_matrix(prefix+"_bootstrap_jsdiv_covar_res_0diag.txt",      self.jsdiv_cov ,name_num_list,name_num_list, zero_diag=True)

            output_matrix(prefix+"_bootstrap_chisq_covar_res.txt",            self.chisq_cov ,name_num_list,name_num_list)
            output_matrix(prefix+"_bootstrap_chisq_covar_res_0diag.txt",      self.chisq_cov ,name_num_list,name_num_list, zero_diag=True)
        except: pass 


    def output(self):
        prefix=str(self.run_params1.resfile_fn)+str(self.run_params2.resfile_fn)
        maxfile_kl = open(prefix+"_max_kldiv.txt",'w')
        sumfile_kl = open(prefix+"_sum_kldiv.txt",'w')
        terfile_kl = open(prefix+"_ter_kldiv.txt",'w')
        allfile_kl = open(prefix+"_all_kldiv.txt",'w')
        maxfile_chisq = open(prefix+"_max_chisq.txt",'w')
        sumfile_chisq = open(prefix+"_sum_chisq.txt",'w')
        terfile_chisq = open(prefix+"_ter_chisq.txt",'w')
        allfile_chisq = open(prefix+"_all_chisq.txt",'w')
        maxfile_js = open(prefix+"_max_jsdiv.txt",'w')
        sumfile_js = open(prefix+"_sum_jsdiv.txt",'w')
        terfile_js = open(prefix+"_ter_jsdiv.txt",'w')
        allfile_js = open(prefix+"_all_jsdiv.txt",'w')
        pymolfile_kl = open(prefix+"_kldiv.pml",'w')
        pymolfile_js = open(prefix+"_jsdiv.pml",'w')
        PDB_input = self.run_params1.pdbfile
        # first, output data for each bootstrap separately
        bootstraps = self.bootstrap_sets
        for mybootstrap in range(bootstraps):
            sumfile_kl_boot = open(prefix+"_sum_kldiv_bootstrap"+str(mybootstrap)+".txt",'w')
            allfile_kl_boot = open(prefix+"_all_kldiv_bootstrap"+str(mybootstrap)+".txt",'w')
            sumfile_js_boot = open(prefix+"_sum_jsdiv_bootstrap"+str(mybootstrap)+".txt",'w')
            allfile_js_boot = open(prefix+"_all_jsdiv_bootstrap"+str(mybootstrap)+".txt",'w')
            counter1 = 0 #residue list counter
            for res1, res2 in zip(self.reslist1, self.reslist2):
                sumfile_kl_boot.write(str(res1.name)+str(res1.num)+" "+str(sum(self.kldiv_res[counter1,mybootstrap,:]))+"\n")
                sumfile_js_boot.write(str(res1.name)+str(res1.num)+" "+str(sum(self.jsdiv_res[counter1,mybootstrap,:]))+"\n")
                allfile_kl_boot.write(str(res1.name)+str(res1.num)+fmt_floats(list((self.kldiv_res[counter1,mybootstrap,:]).flatten()), digits=6, length=9)+"\n")
                allfile_js_boot.write(str(res1.name)+str(res1.num)+fmt_floats(list((self.jsdiv_res[counter1,mybootstrap,:]).flatten()), digits=6, length=9)+"\n")
            sumfile_kl_boot.close()
            sumfile_js_boot.close()
            allfile_kl_boot.close()
            allfile_js_boot.close()
        
        # now, output averages over boostraps #
        counter1 = 0 #residue list counter
        for res1, res2 in zip(self.reslist1, self.reslist2):
            
            print "output "+str(res1.name)+" "+str(res1.num)
            print res1.name, res1.num
            #print arr2str2(chi_pop_hist1[:,:], precision=3, length=6)
            #print
            #print arr2str2(chi_pop_hist2[:,:], precision=3, length=6)
            #print
            
            maxfile_kl.write(str(res1.name)+str(res1.num)+" "+str(amax(average(self.kldiv_res[counter1,:,:],axis=0)))+"\n") #average over boostraps
            sumfile_kl.write(str(res1.name)+str(res1.num)+" "+str(sum(average(self.kldiv_res[counter1,:,:],axis=0)))+"\n")
            terfile_kl.write(str(res1.name)+str(res1.num)+" "+str((average(self.kldiv_res[counter1,:,-1],axis=0)))+"\n")
            allfile_kl.write(str(res1.name)+str(res1.num)+" "+fmt_floats(list((average(self.kldiv_res[counter1,:,:],axis=0)).flatten()), digits=6, length=9)+"\n")
            
            maxfile_js.write(str(res1.name)+str(res1.num)+" "+str(amax(average(self.jsdiv_res[counter1,:,:],axis=0)))+"\n") #average over boostraps
            sumfile_js.write(str(res1.name)+str(res1.num)+" "+str(sum(average(self.jsdiv_res[counter1,:,:],axis=0)))+"\n")
            terfile_js.write(str(res1.name)+str(res1.num)+" "+str((average(self.jsdiv_res[counter1,:,-1],axis=0)))+"\n")
            allfile_js.write(str(res1.name)+str(res1.num)+" "+fmt_floats(list((average(self.jsdiv_res[counter1,:,:],axis=0)).flatten()), digits=6, length=9)+"\n")
            
            
            maxfile_chisq.write(str(res1.name)+str(res1.num)+" "+str(amax(average(self.chisq_res[counter1,:,:],axis=0)))+"\n") #average over boostraps
            sumfile_chisq.write(str(res1.name)+str(res1.num)+" "+str(sum(average(self.chisq_res[counter1,:,:],axis=0)))+"\n")
            terfile_chisq.write(str(res1.name)+str(res1.num)+" "+str(average(self.chisq_res[counter1,:,-1],axis=0))+"\n")
            allfile_chisq.write(str(res1.name)+str(res1.num)+" "+fmt_floats(list((average(self.chisq_res[counter1,:,:],axis=0)).flatten()), digits=6, length=9)+"\n")
            counter1 += 1

        ### Write Pymol Session Files ####
        print "Writing Pymol Session Files"
	
	postfix = ""
        prefix=str(self.run_params1.resfile_fn)+str(self.run_params2.resfile_fn)
        kl_pdbfile = prefix + "_kldiv" + postfix + ".pdb"
        js_pdbfile = prefix + "_jsdiv" + postfix + ".pdb"
        
        pymolfile_kl.write("from pymol import cmd"+"\n")
        pymolfile_kl.write("load "+str(kl_pdbfile)+", system \n")
        pymolfile_kl.write("preset.b_factor_putty('system')"+"\n")  #,_self=cmd"+"\n")
        pymolfile_kl.write("sele b0, b < 0.00001"+"\n")
        pymolfile_kl.write("cmd.color(5278, 'b0')"+"\n")
        pymolfile_kl.write("cmd.disable('b0')"+"\n")
        pymolfile_kl.write("cmd.bg_color('white')"+"\n")
        
        pymolfile_js.write("from pymol import cmd"+"\n")
        pymolfile_js.write("load "+str(js_pdbfile)+", system \n")
        pymolfile_js.write("preset.b_factor_putty('system'),_self=cmd"+"\n")
        pymolfile_js.write("sele b0, b < 0.00001"+"\n")
        pymolfile_js.write("cmd.color(5278, 'b0')"+"\n")
        pymolfile_js.write("cmd.disable('b0')"+"\n")
        pymolfile_js.write("cmd.bg_color('white')"+"\n")

        maxfile_kl.close()
        sumfile_kl.close()
        terfile_kl.close()
        allfile_kl.close()
        maxfile_js.close()
        sumfile_js.close()
        terfile_js.close()
        allfile_js.close()
        maxfile_chisq.close()
        sumfile_chisq.close()
        terfile_chisq.close()
        allfile_chisq.close()
	pymolfile_kl.close()
	pymolfile_js.close()
    
    #########################################################################################################################################
    ### _KL_innerloop_calc_ : For two vectors of pdf's (maximum 6 dimensions), calculates 1st-order KLdiv, JSdiv, Chi-Squared  
    #########################################################################################################################################
    
    def _KL_innerloop_calc_(self, chi_pop_hist1, chi_pop_hist2, chi_counts1, chi_counts2, nchi_to_use):
                kldiv1 = zeros((6),float64)
                kldiv2 = zeros((6),float64)
                jsdiv2 = zeros((6),float64)
                kldiv4 = zeros((6),float64)
                jsdiv4 = zeros((6),float64)
                if(self.run_params1.options.grassberger == "no" and self.run_params1.options.abs == "no"):
                    for mychi in range(nchi_to_use):
                        try:
                            pi = chi_pop_hist1[mychi,:]
                            pj = chi_pop_hist2[mychi,:]
                        except IndexError:
                            print "Cannot find torsion angle histogram for this residue, please check residue list and input torsion angles."
                            exit 
                        #kldiv1[mychi] = sum(pi[pj > SMALL] * log (pi[pj > 0 + SMALL]/pj[pj > 0 + SMALL]), axis=-1)
                        #Kullback-Leibler Divergence
                        #print "kldiv test"+str( sum(pj[pi > SMALL] * log (pj[pi > 0 + SMALL]/pi[pi > 0 + SMALL]), axis=-1))
                        kldiv2[mychi] = sum(pj[pi > SMALL] * log (pj[pi > 0 + SMALL]/pi[pi > 0 + SMALL]), axis=-1)
                        #mypi = pi[(pi > SMALL or pj > SMALL)]
                        #mypj = pj[(pi > SMALL or pj > SMALL)]
                        #print "pi or pj > 0, shape:"
                        #print mypi.shape
                        #print mypj.shape
                        #Jensen-Shannon Divergence
                        #print "jsdiv test"+str(0.5 * sum(pj * log (pj/(0.5*pi+0.5*pj+SMALL)), axis=-1) + 0.5 * sum(pi * log (pi/(0.5*pi + 0.5*pj+SMALL)), axis=-1))
                        jsdiv2[mychi] = 0.5 * sum(pj * log (pj/(0.5*pi+0.5*pj+SMALL)), axis=-1) + 0.5 * sum(pi * log (pi/(0.5*pi + 0.5*pj+SMALL)), axis=-1)

                if(self.run_params1.options.grassberger == "no" and self.run_params1.options.abs == "yes"):
                    for mychi in range(nchi_to_use):
                        try:
                            pi = chi_pop_hist1[mychi,:]
                            pj = chi_pop_hist2[mychi,:]
                        except IndexError:
                            print "Cannot find torsion angle histogram for this residue, please check residue list and input torsion angles."
                            exit 
                        #kldiv1[mychi] = sum(abs(pi[pj > SMALL] * log (pi[pj > 0 + SMALL]/pj[pj > 0 + SMALL])), axis=-1)
                        kldiv2[mychi] = sum(abs(pj[pi > SMALL] * log (pj[pi > 0 + SMALL]/pi[pi > 0 + SMALL])), axis=-1)
                if(self.run_params1.options.grassberger == "yes"):
                    for mychi in range(nchi_to_use):
                        #kldiv1 = sum(Grassberger_KLdiv(chi_counts1, chi_counts2, self.numangles1, self.numangles2), axis=-1)
                        #check to make sure the counts are there before running calculation
                        try:  
                            counts_i = chi_counts1[mychi,:]
                            counts_j = chi_counts2[mychi,:]
                        except IndexError:
                            print "Cannot find torsion angle histogram for this residue, please check residue list and input torsion angles."
                            exit 
                        kldiv2[mychi] = Grassberger_KLdiv(chi_counts1[mychi,:], chi_counts2[mychi,:], self.numangles1, self.numangles2)
                        jsdiv2[mychi] = 0.5 * Grassberger_KLdiv(chi_counts1[mychi,:], int16(0.5 * chi_counts2[mychi,:] + 0.5 * chi_counts1[mychi,:]), self.numangles1,self.numangles2) + \
                                        0.5 *  Grassberger_KLdiv(chi_counts2[mychi,:], int16(0.5 * chi_counts2[mychi,:] + 0.5 * chi_counts1[mychi,:]),self.numangles1,self.numangles2) 
                #if(self.run_params1.options.reference == 0):
                #    kldiv4 = (kldiv1 + kldiv2) / 2
                #if(self.run_params1.options.reference == 1):
                kldiv4 = kldiv2
                jsdiv4 = jsdiv2
                #else:
                #    kldiv4 = kldiv1
                #print "kldiv shape:"+str(kldiv4.shape)
                chisq = (chi_pop_hist1-chi_pop_hist2)**2
                chisq = sqrt(chisq[:,:].sum(axis=-1))
                chisq4 = zeros((6),float64)
                chisq4[:nchi_to_use] = chisq
                return (kldiv4, jsdiv4, chisq4)

    #########################################################################################################################################
    ### __KL_resi_resj__ : For a single residue, calculates KLdiv, JSdiv, Chi-Squared
    #########################################################################################################################################
    def __KL_resi_resj__(self,res1,res2):
        #We are interested in local population shifts caused by perturbations that reflect subtle changes in structure and/or dynamics. We can visualize these most readily using the first-order terms from our expansion. Consider the terms in the Kullback-Leibler Divergence arising from a particular degree of freedom. These we will denote the "local" Kullback-Leibler Divergence and provide an information-theoretic, quantitative measure of the extent to which the p.d.f. for a given degree of freedom deviates from the equilibrium p.d.f, i.e. a much less biased measure than the chi-squared statistic. 

        #\begin{equation}
        # KL_n = \sum_{i}^{nbins}p_{i}\ln \frac{p_i}{p_{i^*}}
        # \end{equation}

        # Note that for the local K-L Div over a set of d.o.f, we might want to include second-order terms to obtain more accurate (though not necessarily as stable) results. The biggest impact would be calculating these second-order terms within a single residue's torsions to improve the estimate of the Kullback-Leibler divergence for that residue. Currently, however, we focus on first-order terms, as these are most readily calculated.
        self.kldiv[:,:] = 0
        self.jsdiv[:,:] = 0
        self.chisq[:,:] = 0
        #nchi1, nchi2 = res1.chi_pop_hist.shape[1], res2.chi_pop_hist.shape[1]
        nchi1 = res1.get_num_chis(res1.name)* (1 - self.backbone_only) + self.run_params1.phipsi 
        nchi2 = res2.get_num_chis(res2.name)* (1 - self.backbone_only) + self.run_params1.phipsi 
        if self.run_params1.backbone == "coarse_phipsi":
            nchi1 = 1
        if self.run_params2.backbone == "coarse_phipsi":
            nchi2 = 1
        nchi_to_use = min(nchi1, nchi2)                
        if(self.which_runs_ref == None): #this indicates we aren't doing bootstrap resampling of reference during this function call
            for mybootstrap in range(self.bootstrap_sets): #bootstraps of target
                print "bootstrap: "+str(mybootstrap)
                #print res1.pop_hist, res1.chi_counts # bootstraps, nchi, nbins
                # reference now always has bootstrap_sets = 1
                ## TARGET DISTRIBUTION
                chi_counts2 = res2.chi_counts[mybootstrap,:nchi_to_use,:]  #grab data from a bootstrap sample of the target data
                chi_pop_hist2 =  chi_counts2 / (1.0 * resize(sum(chi_counts2, axis = -1), chi_counts2.shape))
                #REFERENCE DISTRIBUTION
                chi_counts1 = sum(res1.chi_counts[:,:nchi_to_use,:],axis=0)  #average over bootstrap dimension for reference in this case
                chi_pop_hist1 = chi_counts1 / (1.0 * resize(sum(chi_counts1, axis = -1),chi_counts1.shape))
                chi_pop_hist1[chi_pop_hist1==0] = SMALL * 0.5 # to avoid NaNs
                chi_pop_hist2[chi_pop_hist2==0] = SMALL * 0.5 # to avoid NaNs
                #RUN CALCULATION, AVOIDING NaNs
                print "counts i:"+str(int32(chi_counts1))+"\n"
                print "counts j:"+str(int32(chi_counts2))+"\n"
                print "pi      :"+str(chi_pop_hist1)+"\n"
                print "pj      :"+str(chi_pop_hist2)+"\n"
                (kldiv4, jsdiv4, chisq4) = self._KL_innerloop_calc_(chi_pop_hist1, chi_pop_hist2, chi_counts1, chi_counts2, nchi_to_use)
                self.kldiv[mybootstrap,:nchi_to_use] = kldiv4[:nchi_to_use]
                self.jsdiv[mybootstrap,:nchi_to_use] = jsdiv4[:nchi_to_use]
                self.chisq[mybootstrap,:nchi_to_use] = chisq4[:nchi_to_use]
                
        else:  #usually here, res1 and res2 are from the same residue list, just different "runs" will be used through which_runs_ref and its complement
            print "res1 chi counts shape:"+str(res1.chi_counts.shape)
            minangles = min(sum(res1.chi_counts[0,0]), sum(res2.chi_counts[0,0])) #minimum of number of datapoints in each 
            subsample_block_size = int(minangles / (self.run_params1.blocks_ref * 1.0))
            subsample_choose_ref = int(self.run_params1.blocks_ref / 2.0) #subsample_choose_ref should be half of blocks_ref
            #bootstrap_sets =  int(misc.comb(self.run_params1.blocks_ref,subsample_choose_ref)) #like bootstrap_sets, just for the reference, so it gets a 
            print "bootstrap sets reference: "+str(self.bootstrap_sets)
            #bootstrap_sets = self.bootstrap_sets
            which_runs_ref = []
            
            #NOTE THIS NEW WAY OF DOING BOOTSTRAPS IS PURE PYTHON AND MEMORY EFFICIENT
            for mybootstrap in range(self.bootstrap_sets): 
                kldiv1 = zeros((6),float64)
                kldiv2 = zeros((6),float64)
                jsdiv2 = zeros((6),float64)
                kldiv4 = zeros((6),float64)
                jsdiv4 = zeros((6),float64)
                chi_counts1 = zeros((nchi_to_use,self.nbins),float64)
                chi_counts2 = zeros((nchi_to_use,self.nbins),float64)
                for myblock in range(self.subsample_choose_ref):
                    #we're only selecting the histogram of bin number = nbins, which is first variable entry "1"; a more optimal bin size could be chosen by computing kldiv from multiple
                    #print "which runs ref:"+str(self.which_runs_ref[mybootstrap,myblock])
                    #print "which runs ref comp:"+str(self.which_runs_ref_complement[mybootstrap,myblock])
                    #print  "shape of source1:" + str(shape(chi_counts1[:nchi_to_use,:]))
                    #print  "shape of target1:" + str(shape(res1.chi_counts_sequential_varying_bin_size[1,self.which_runs_ref[mybootstrap,myblock],:nchi_to_use,:self.nbins]))
                    #print  "shape of source2:" + str(shape(chi_counts2[:nchi_to_use,:]))
                    #print  "shape of target2:" + str(shape(res2.chi_counts_sequential_varying_bin_size[1,self.which_runs_ref_complement[mybootstrap,myblock],:nchi_to_use,:self.nbins]))
                    chi_counts1 += res1.chi_counts_sequential_varying_bin_size[1,self.which_runs_ref[mybootstrap,myblock],:nchi_to_use,:self.nbins]
                    chi_counts2 += res2.chi_counts_sequential_varying_bin_size[1,self.which_runs_ref_complement[mybootstrap,myblock],:nchi_to_use,:self.nbins]
                    
                    #...if we wanted to instead grab from the original bins data .... 
                    #leftstop1  = subsample_block_size * (self.which_runs_ref_complement[mybootstrap,myblock])
                    #rightstop1 = subsample_block_size * (self.which_runs_ref_complement[mybootstrap,myblock] + 1)
                    #leftstop_new1 = subsample_block_size * (myblock)
                    #rightstop_new1 = subsample_block_size * (myblock + 1)
                    #leftstop2  = subsample_block_size * (self.which_runs_ref[mybootstrap,myblock])
                    #rightstop2 = subsample_block_size * (self.which_runs_ref[mybootstrap,myblock] + 1)
                    #leftstop_new2 = subsample_block_size * (myblock)
                    #rightstop_new2 = subsample_block_size * (myblock + 1)
                    #print "minangles:"+str(minangles)
                    #print "subsample block size:"+str(subsample_block_size)
                    #print "leftstop reference runs: "+str(leftstop1)
                    #print "rightstop reference runs: "+str(rightstop1)
                    #print "leftstop new: "+str(leftstop1)
                    #print "rightstop new: "+str(rightstop1)
                    #print "counts full chi 0:   "+str(res1.chi_counts[:,0,:])
                    #print "counts example:"+str(res1.chi_counts[:,:nchi_to_use,leftstop1:rightstop1])
                
                #print "chi_counts1:"+str(int16(chi_counts1))
                #print "chi_counts2:"+str(int16(chi_counts2))
                chi_pop_hist1 = chi_counts1 / (1.0 * resize(sum(chi_counts1, axis = -1),chi_counts1.shape)) #normalization
                chi_pop_hist2 = chi_counts2 / (1.0 * resize(sum(chi_counts2, axis = -1),chi_counts2.shape)) #normalization
                chi_pop_hist1[chi_pop_hist1==0] = SMALL * 0.5 # to avoid NaNs
                chi_pop_hist2[chi_pop_hist2==0] = SMALL * 0.5 # to avoid NaNs
                (kldiv, jsdiv, chisq) = self._KL_innerloop_calc_(chi_pop_hist1, chi_pop_hist2, chi_counts1, chi_counts2, nchi_to_use)
                self.kldiv[mybootstrap,:nchi_to_use] = kldiv[:nchi_to_use]
                self.jsdiv[mybootstrap,:nchi_to_use] = jsdiv[:nchi_to_use]
                self.chisq[mybootstrap,:nchi_to_use] = chisq[:nchi_to_use]

        print "kldiv:\n"+str(self.kldiv)
        print "KLDIV", res1.name, res1.num, fmt_floats(list((average(self.kldiv,axis=0)).flatten()), digits=6, length=9), "JSDIV", fmt_floats(list((average(self.jsdiv,axis=0)).flatten()), digits=6, length=9)
        return self.kldiv, self.chisq, self.jsdiv



    #########################################################################################################################################
    ##### do second order and __calc_pair_mutdiv__: Mutual Divergence For All Pairs of Torsions: 2nd term in Kullback-Leibler Divergence Expansion ################## 
    #########################################################################################################################################


    # Calculate the mutual information divergence between all pairs of residues.
    # The MI between a pair of residues is the sum of the MI between all combinations of res1-chi? and res2-chi?.
    # Returns mut_info_res_matrix, mut_info_uncert_matrix
    def do_second_order(self):
        name_num_list = make_name_num_list(self.reslist1)
        mut_div_res_matrix, mut_div_uncert_matrix, mut_div_res_matrix_different_sims, dKLtot_dresi_dresj_matrix = self.__calc_pair_mutdiv__()
        mut_div_res_matrix_avg = zeros(mut_div_res_matrix.shape[1:3],float32)
        mut_div_res_sumoverchis_matrix_avg = zeros(mut_div_res_matrix.shape[2:3],float32)
        for i in range(mut_div_res_matrix.shape[1]):
         for j in range(i, mut_div_res_matrix.shape[2]):
            for k in range(mut_div_res_matrix.shape[3]):
                for m in range(mut_div_res_matrix.shape[4]):
                    if(sum(mut_div_res_matrix[:,i,j,k,m]) > 0):
                        mutdiv_boots = mut_div_res_matrix[:,i,j,k,m].copy()
                        mutdiv_boots[mutdiv_boots < 0] = 0 #use negative and zero values for significance testing but not in the average
                        #zero values of mutdiv_boots will include those zeroed out in the permutation test
                        mutdiv = mut_div_res_matrix_avg[i,j,k,m] = average(mutdiv_boots[mutdiv_boots > 0 ],axis=0) 
                        #uncert = mut_div_uncert_matrix_avg[i,j,k,m] = sqrt(cov(mut_div_res_matrix[:,i,j,k,m]) / (run_params.num_sims))
                        #use vfast_cov_for1D_boot as a fast way to calculate a stdev instead of the function std()
                        uncert = pval = None
                        if(mut_div_res_matrix.shape[0] >= 10):
                            uncert = mut_div_uncert_matrix_avg[i,j,k,m] = sqrt((vfast_cov_for1D_boot(reshape(mut_div_res_matrix[:,i,j,k,m],(mut_div_res_matrix.shape[0],1)))[0,0]) / (run_params.num_sims))
                            #pval = mut_div_pval_matrix[i,j,k,m] =(stats.wilcoxon(mut_div_res_matrix[:,i,j,k,m]+SMALL))[1] / 2.0 #offset by SMALL to ensure no values of exactly zero will be removed in the test
                            #if pval == 0.00 and (i != j): pval = 1.0
                            #if pval <= 0.05:
                            #    mut_div_res_matrix_sig_05[:,i,j,k,m] = mutdiv_boots.copy()
                            #if pval <= 0.01:
                            #    mut_div_res_matrix_sig_01[:,i,j,k,m] = mutdiv_boots.copy()

                        else:
                            uncert = mut_div_uncert_matrix_avg[i,j,k,m] = sqrt((vfast_cov_for1D_boot(reshape(mut_div_res_matrix[:,i,j,k,m],(mut_div_res_matrix.shape[0],1)))[0,0]) / (run_params.num_sims))
                        #    pval = 0
                            
                        if(True):
                        #if pval <= 0.05 and (i != j or k == m):
                            mut_div_res_sumoverchis_matrix_avg[i,j] += mutdiv
                        #    mut_div_res_sumoverchis_matrix_sig_[:,i,j] += mutdiv_boots.copy()

                        #    if(i == j):
                        #        tot_ent_sig05 += mutdiv_boots.copy()
                        #    else:
                        #        tot_ent_sig05 -= mutdiv_boots.copy()


         prefix=str(self.run_params1.resfile_fn)+str(self.run_params2.resfile_fn)
         output_matrix(prefix+"_bootstrap_avg_mutdiv_res_sum.txt",            mut_div_res_sumoverchis_matrix_avg ,name_num_list,name_num_list)
         output_matrix(prefix+"_bootstrap_avg_mutdiv_res_sum_0diag.txt",      mut_div_res_sumoverchis_matrix_avg ,name_num_list,name_num_list, zero_diag=True)


    def __calc_pair_mutdiv__(self):  
        rp1 = self.run_params1 #this one is the reference
        rp2 = self.run_params2
        #which_runs2 = rp2.which_runs
        #print rp2.which_runs
        #bootstrap_sets = len(which_runs)
        #initialize the mut info matrix
        check_for_free_mem()
        mut_div_res_matrix = zeros((self.bootstrap_sets, len(self.reslist1),len(self.reslist1),6,6),float32)
        mut_div_res_matrix_different_sims = zeros((self.bootstrap_sets, len(self.reslist1),len(self.reslist1),6,6),float32)
        mut_div_uncert_matrix = zeros((self.bootstrap_sets, len(self.reslist1),len(self.reslist1),6,6),float32)
        dKLtot_dresi_dresj_matrix = zeros((self.bootstrap_sets, len(self.reslist1),len(self.reslist1)),float32)
        bootstrap_choose_star = rp2.bootstrap_choose
        bootstrap_choose = rp1.bootstrap_choose
        numangles_bootstrap = zeros(self.bootstrap_sets)
        numangles_bootstrap_star = zeros(self.bootstrap_sets)
        print "picking n sims at a time from reference, n = "+str(bootstrap_choose)
        print "picking n sims at a time from target, n = "   +str(bootstrap_choose_star)
        print "bootstrap sets = " + str(self.bootstrap_sets)
        #fill bootstraps of ref if reslist1 (the reference) only has one bootstrap, but the target has multiple
        for res_indk_ref, myresk_ref, myres2 in zip(range(len(self.reslist1)), self.reslist1, self.reslist2):
            
            #temp_shape = myresk_ref.bins.shape
            #new_shape = array(temp_shape,int64)
            #print new_shape
            #new_shape[2] = myres2.bins.shape[2]
            #myresk.bins = zeros(new_shape, int8) #re-initialize
        
            if(myresk_ref.bins.shape[2] == myres2.bins.shape[2]): #reference has the same number of bootstraps in this case ...
                numangles_bootstrap = repeat(myresk_ref.numangles_bootstrap, self.bootstrap_sets, axis=0) #replicate this over all bootstrap sets for reference
                numangles_bootstrap_star = myres2.numangles_bootstrap #this one is already is replicated over bootstrap sets
                #myresk_ref.bins = temp_bins[:,:,:,:] 
            else:
                assert(myresk_ref.bins.shape[2] == 1) # ... otherwise, reference must have only one bootstrap 
                temp_bins = myresk_ref.bins[:,:,:,:]
                myresk_ref.bins = repeat(temp_bins, self.bootstrap_sets, axis=2)
                myresk_ref.chi_counts = repeat(myresk_ref.chi_counts[:,:,:], self.bootstrap_sets, axis=0) #replicate this over all boostrap sets
                numangles_bootstrap = repeat(myresk_ref.numangles_bootstrap, self.bootstrap_sets, axis=0)  #replicate this over all bootstrap sets for reference
                numangles_bootstrap_star = myres2.numangles_bootstrap #this one is already is replicated over bootstrap sets
                #for mychi in temp_shape[0]:
                #    for permut in temp_shape[1]:
                #        myresk.bins[mychi,permut,:,:] = repeat(temp_bins[mychi,permut,0,:], self.bootstrap_sets) # to make sure this repeat works properly
                    
        print "numangles_bootstrap reference: "+str(numangles_bootstrap)
        print "numangles_bootstrap target: "+str(numangles_bootstrap_star)
        print "bins shape :" + str(self.reslist1[0].bins)
        print "bins_star shape: "+ str(self.reslist2[0].bins)
        #Loop over the residue list
        
        for res_ind1, myres1, res_ind1_ref, myres1_ref in zip(range(len(self.reslist2)), self.reslist2, range(len(self.reslist1)), self.reslist1):
           for res_ind2, myres2, res_ind2_ref, myres2_ref in zip(range(res_ind1, len(self.reslist1)), self.reslist2[res_ind1:], range(res_ind1, len(reslist1)), self.reslist1[res_ind1:]):
            print "\n#### Working on residues %s and %s (%s and %s):" % (myres1.num, myres2.num, myres1.name, myres2.name) , utils.flush()
            max_S = 0.
            if(OFF_DIAG == 1):
              for mychi1 in range(min(myres1_ref.nchi, myres1.nchi)):  #whichever has fewer chis
                 for mychi2 in range(min(myres2_ref.nchi, myres2.nchi)):   #whichever has fewer chis
                     print 
                     print "%s %s , %s %s chi1/chi2: %d/%d" % (myres1.name,myres1.num, myres2.name,myres2.num, mychi1+1,mychi2+1)
                     check_for_free_mem()
                     mutinf_thisdof = var_mi_thisdof = mutinf_thisdof_different_sims = dKLtot_dKL1_dKL2 = 0 #initialize
                     cross_mutinf_thisdof = cross_var_mi_thisdof = cross_mutinf_thisdof_different_sims = cross_dKLtot_dKL1_dKL2 = pvalue_cross = 0 #initialize
                     angle_str = ("%s_chi%d-%s_chi%d"%(myres1, mychi1+1, myres2, mychi2+1)).replace(" ","_")
                     
                     mutdiv_thisdof = zeros((self.bootstrap_sets), float64)  #no permutations
                     if((res_ind1 != res_ind2) or (res_ind1 == res_ind2 and mychi1 > mychi2)):
                         mutinf_thisdof, var_mi_thisdof, mutinf_thisdof_different_sims, dKLtot_dKL1_dKL2 = \
                                     calc_excess_mutinf(myres1_ref.chi_counts[:,mychi1,:],myres2_ref.chi_counts[:,mychi2,:],\
                                                        myres1_ref.bins[mychi1,:,:,:], myres2_ref.bins[mychi2,:,:,:], \
                                                        myres1_ref.chi_counts_sequential[:,mychi1,:],\
                                                        myres2_ref.chi_counts_sequential[:,mychi2,:],\
                                                        myres1.simbins[mychi1,:,:,:], myres2.simbins[mychi2,:,:,:],\
                                                        rp2.num_sims, rp2.nbins, numangles_bootstrap,\
                                                        myres1_ref.numangles, rp2.sigalpha, 0,\
                                                        bootstrap_choose, calc_variance=rp2.calc_variance,\
                                                        which_runs=rp2.which_runs,pair_runs=rp2.pair_runs,\
                                                        calc_mutinf_between_sims=rp2.calc_mutinf_between_sims,\
                                                        file_prefix=angle_str, plot_2d_histograms=False, adaptive_partitioning = adaptive_partitioning)
                                     

                         #now the next one is for the cross-term: pij ln (pij_star / (pi_star * pj_star)), 
                         #or we can think of this as the the pij-weighted ensemble average <ln (pij_star / (pi_star * pj_star))>_pij
                         #no permutations

                         #reset globals ... these are for efficiency of dihedral_mutent.py, I sacrificed transparency for efficiency
                         count_matrix = None
                         
                         numangles_bootstrap_matrix = None
                         numangles_bootstrap_vector = None

                         cross_mutinf_thisdof, cross_var_mi_thisdof, cross_mutinf_thisdof_different_sims, cross_var_thisdof_different_sims, \
                             mutinf_multinomial, mutinf_multinomial_sequential, pvalue_cross, cross_dKLtot_dKL1_dKL2  = \
                                     calc_mutinf_corrected(myres1_ref.chi_counts[:,mychi1,:],myres2_ref.chi_counts[:,mychi2,:],\
                                                        myres1_ref.bins[mychi1,:,:,:], myres2_ref.bins[mychi2,:,:,:], \
                                                        myres1_ref.chi_counts_sequential[:,mychi1,:],\
                                                        myres2_ref.chi_counts_sequential[:,mychi2,:],\
                                                        myres1_ref.simbins[mychi1,:,:,:], myres2_ref.simbins[mychi2,:,:,:],\
                                                        rp1.num_sims, rp1.nbins, numangles_bootstrap,\
                                                        myres1_ref.numangles, calc_variance=rp1.calc_variance,\
                                                        bootstrap_choose=bootstrap_choose, permutations=0, \
                                                        which_runs=rp2.which_runs,pair_runs=rp2.pair_runs,\
                                                        calc_mutinf_between_sims=rp2.calc_mutinf_between_sims,\
                                                        file_prefix=angle_str, plot_2d_histograms=False, \
                                                        adaptive_partitioning = adaptive_partitioning, \
                                                        chi_counts1_star = myres1.chi_counts[:,mychi1,:], chi_counts2_star = myres2.chi_counts[:,mychi2,:],\
                                                        bins1_star = myres1.bins[mychi1,:,:,:], bins2_star = myres2.bins[mychi2,:,:,:], \
                                                        chi_counts_sequential1_star = myres1.chi_counts_sequential[:,mychi1,:],\
                                                        chi_counts_sequential2_star = myres2.chi_counts_sequential[:,mychi2,:],\
                                                        bins1_sequential_star = myres1.simbins[mychi1,:,:,:], \
                                                        bins2_sequential_star = myres2.simbins[mychi2,:,:,:], \
                                                        numangles_star = myres1.numangles, \
                                                        numangles_bootstrap_star = numangles_bootstrap_star, bootstrap_choose_star = bootstrap_choose_star )
                                                        
                                                        # chi_counts1_star=None, chi_counts2_star=None, bins1_star=None, bins2_star=None,chi_counts_sequential1_star=None, chi_counts_sequential2_star=None, bins1_sequential_star=None, bins2_sequential_star=None, numangles_bootstrap_star=None ):
                         sigalpha = rp1.sigalpha
                         ########################
                         print "mutinf thisdof "+str(mutinf_thisdof)
                         print "crossinf thisdof"+str(cross_mutinf_thisdof[:,0])
                         mutdiv_thisdof = cross_mutinf_thisdof[:,0] - mutinf_thisdof # what the mutual divergence really is 
                         ########################
                         #pvalue_toolow1 = pvalue > sigalpha
                         #pvalue_toolow2 = pvalue_cross > sigalpha
                         #if(sum(pvalue_toolow2) > 0):
                         #    print "one or more values were not significant!"
                         #mutdiv_thisdof *= (1.0 - pvalue_toolow1 * pvalue_toolow2 * 1.0) #zeros elements with pvalues in either term that are below threshold
                         #mutdiv_thisdof *= (1.0 - pvalue_toolow2 * 1.0) #zeros elements with pvalues in either term that are below threshold
                         print "mutdiv thisdof: "+str(mutdiv_thisdof)
                     
                     if(res_ind1 == res_ind2 and mychi1 == mychi2):
                         #mut_div_res_matrix[:,res_ind1, res_ind2, mychi1, mychi2] = myres1.entropy[:,mychi1]
                         #mut_div_uncert_matrix[:,res_ind1, res_ind2, mychi1, mychi2] = myres1.var_ent[:,mychi1]
                         max_S = 0
                         ## still need to calc dKLdiv here
                         dKLtot_dresi_dresj_matrix[:,res_ind1, res_ind2] += 0 #myres1.dKLtot_dchis2[:,mychi1]

                     elif(res_ind1 == res_ind2 and mychi1 > mychi2):
                         mut_div_res_matrix[:,res_ind1, res_ind2, mychi1, mychi2] = mut_div_res_matrix[:,res_ind1, res_ind2, mychi2, mychi1]
                         mut_div_uncert_matrix[:,res_ind1, res_ind2, mychi1, mychi2] = mut_div_uncert_matrix[:,res_ind1, res_ind2, mychi2, mychi1]
                         max_S = 0
                         ## still need to calc dKLdiv here
                         dKLtot_dresi_dresj_matrix[:,res_ind1, res_ind2] += dKLtot_dKL1_dKL2
                     
                     else:
                         S = 0
                         mychi1 = int(mychi1)
                         mychi2 = int(mychi2)
                         blah = myres1.chi_pop_hist[:,mychi1,:]
                         blah = myres2.chi_pop_hist[:,mychi2,:]
                         blah = myres1.bins[mychi1,:,:,:]
                         blah =  myres2.bins[mychi2,:,:,:]
                         blah =  myres1.chi_pop_hist_sequential[:,mychi1,:]
                         dKLtot_dresi_dresj_matrix[:,res_ind1, res_ind2] += dKLtot_dKL1_dKL2
                         dKLtot_dresi_dresj_matrix[:,res_ind2, res_ind1] += dKLtot_dKL1_dKL2 #note res_ind1 neq res_ind2 here
                         mut_div_res_matrix[:,res_ind1 , res_ind2, mychi1, mychi2] = mutdiv_thisdof
                         mut_div_uncert_matrix[:,res_ind1, res_ind2, mychi1, mychi2] = var_mi_thisdof + cross_var_mi_thisdof # assume sum of variances
                         mut_div_res_matrix_different_sims[:,res_ind1, res_ind2, mychi1, mychi2] = mutinf_thisdof_different_sims #not changing this one yet
                         max_S = max([max_S,S])
                         mut_div_res_matrix[:,res_ind2, res_ind1, mychi2, mychi1] = mut_div_res_matrix[:,res_ind1, res_ind2, mychi1, mychi2] #symmetric matrix
                         #mut_div_uncert_matrix[res_ind1, res_ind2] = mut_div_uncert_matrix[res_ind1, res_ind2]
                         mut_div_uncert_matrix[:,res_ind2, res_ind1, mychi2, mychi1] = mut_div_uncert_matrix[:,res_ind1, res_ind2, mychi1, mychi2] #symmetric matrix
                         mut_div_res_matrix_different_sims[:,res_ind2, res_ind1, mychi2, mychi1] = mut_div_res_matrix_different_sims[:,res_ind1, res_ind2, mychi1, mychi2] #symmetric matrix
                         
            #print "mutinf=%.3f (uncert=%.3f; max(S)=%.3f" % (average((mut_div_res_matrix[:,res_ind1, res_ind2, : ,:]).flatten()), sum((mut_div_uncert_matrix[0,res_ind1, res_ind2, :, :]).flatten()), max_S),
            #if max_S > 0.26: print "#####",
            print
        
        return mut_div_res_matrix, mut_div_uncert_matrix, mut_div_res_matrix_different_sims, dKLtot_dresi_dresj_matrix


    #######################################################################################################################################
    ##### __init__: Set Up Kullback-Leibler Divergence Class for  Kullback-Leibler Divergence Expansion ################## 
    #########################################################################################################################################

    def __init__ (self,run_params1,run_params2,reslist1,reslist2,which_runs_ref = None, which_runs_ref_complement = None, do_bootstraps_with_complements = False):
        self.run_params1 = run_params1
        self.run_params2 = run_params2
        self.reslist1 = reslist1
        self.reslist2 = reslist2
        self.which_runs_ref = which_runs_ref
        self.which_runs_ref_complement = which_runs_ref_complement
        self.backbone_only = run_params1.backbone_only
        res1 = reslist1[0] # just to get some parameters
        res2 = reslist2[0] # just to get some parameters
        if(do_bootstraps_with_complements == False):
            self.bootstrap_sets = res2.chi_counts.shape[0] #since reference doesn't use bootstraps in this case
        else:
            self.subsample_choose_ref = int(self.run_params1.blocks_ref / 2.0) #subsample_choose_ref should be half of blocks_ref
            self.bootstrap_sets =  int(misc.comb(self.run_params1.blocks_ref,self.subsample_choose_ref)) #like bootstrap_sets, just for the reference, so it gets a 
        print "bootstrap_sets:"+str(self.bootstrap_sets)
        self.nbins = res1.chi_counts.shape[-1]
        self.numangles1 = sum(res1.numangles)
        self.numangles1_bootstrap = res1.numangles_bootstrap
        self.numangles2 = sum(res2.numangles)
        self.numangles2_bootstrap = res2.numangles_bootstrap
        self.min_numangles = min(self.numangles1, self.numangles2)
        self.kldiv = zeros((self.bootstrap_sets, 6), float64) # KLdiv for one residue: max torsions per res=6
        self.chisq = zeros((self.bootstrap_sets, 6), float64) # Chi-Squared for one residue: max torsions per res=6
        self.jsdiv = zeros((self.bootstrap_sets, 6), float64) # Jensen-Shannon Divergence (JSdiv) for one residue: max torsions per res=6
        self.kldiv_res = zeros((len(self.reslist1),self.bootstrap_sets, 6), float64) # KLdiv for all residues
        self.chisq_res = zeros((len(self.reslist1),self.bootstrap_sets, 6), float64) # Chi-Squared for all residues
        self.jsdiv_res = zeros((len(self.reslist1),self.bootstrap_sets, 6), float64) # JSdiv for all residues
        self.sigalpha = min(run_params1.sigalpha, run_params2.sigalpha)              # significance value for significance test
        counter = 0
        for res1, res2 in zip(self.reslist1, self.reslist2):
            #consider refactoring this into a more functional form
            #calculate KLdiv between these two residues... this routine may later include higher-order terms and coord. transforms.
            (self.kldiv_res[counter,:,:],  self.chisq_res[counter,:,:], self.jsdiv_res[counter,:,:]) = self.__KL_resi_resj__(res1,res2)  
            counter += 1
            
        self.kldiv_done = 1 # flag


#######################################################################################################################################
###
###  End of Class KLdiv
###
#######################################################################################################################################
#######################################################################################################################################






######################################################################################################################################
### run_kldiv:  initializes KLdiv class based on options passed in; runs the calculation and statistical filtering
#######################################################################################################################################

def run_kldiv(options, xvg_basedir1, xvg_basedir2, resfile_fn1, resfile_fn2):
    adaptive_partitioning = False #(options.adaptive == "yes")
    phipsi = 0
    backbone_only = 0
    if options.backbone == "phipsichi":
        phipsi = 2;
        backbone_only =0
    if options.backbone == "phipsi":
        phipsi = 2;
        backbone_only = 1
    if options.backbone == "coarse_phipsi":
        phipsi = -1
        backbone_only = 1
        options.binwidth = 90 #override
    bins=arange(-180,180,options.binwidth) #Compute bin edges
    nbins = len(bins)
    num_sims = options.num_sims
    assert(num_sims != None)
    xvgorpdb = "xvg"
    num_structs = None

    if options.bootstrap_set_size == None:
        options.bootstrap_set_size = num_sims
    print "number of simulations to take at a time as bootstrap samples:"+str(options.bootstrap_set_size)+"\n"
    which_runs1 = []
    pair_runs_list = []
    #perhaps should use xuniquePermutations instead 
    for myruns in xuniqueCombinations(range(num_sims), options.num_sims): #only one bootstrap for reference
        which_runs1.append(myruns)

    bootstrap_pair_runs_list = []
    for bootstrap in range((array(which_runs1)).shape[0]):
        pair_runs_list = []
        for myruns2 in xcombinations(which_runs1[bootstrap], 2):
            pair_runs_list.append(myruns2)
        bootstrap_pair_runs_list.append(pair_runs_list)
    pair_runs_array = array(bootstrap_pair_runs_list,int16)

    which_runs2 = []
    pair_runs_list2 = []
    for myruns in xuniqueCombinations(range(num_sims), options.bootstrap_set_size): 
        which_runs2.append(myruns)
            
    bootstrap_pair_runs_list2 = []
    for bootstrap in range((array(which_runs2)).shape[0]):
        pair_runs_list2 = []
        for myruns2 in xcombinations(which_runs2[bootstrap], 2):
            pair_runs_list2.append(myruns2)
        bootstrap_pair_runs_list2.append(pair_runs_list2)
    pair_runs_array2 = array(bootstrap_pair_runs_list2,int16)

    calc_stats  = False
    if(options.bootstrap_set_size <num_sims): 
        calc_stats = True
        
    which_runs_ref = []
    blocks_ref = num_sims  #hardcoded for now
    subsample_choose_ref = int(blocks_ref/2.0) #hardcoded for now
    for myruns in xuniqueCombinations(range(blocks_ref), subsample_choose_ref): 
                which_runs_ref.append(myruns)       #used as follows: which_runs[bootstrap,block_num]
    which_runs_ref = array(which_runs_ref)
    which_runs_ref_complement = []
    for mybootstrap in range(len(which_runs_ref)):
        runs_list = ones((blocks_ref),int16)
        for myblock in range(subsample_choose_ref):
            runs_list[which_runs_ref[mybootstrap,myblock]] = 0
        templist = []
        for myblock in range(blocks_ref):
            if (runs_list[myblock] > 0):
                templist.append(myblock )
        which_runs_ref_complement.append(templist)
    which_runs_ref_complement = array(which_runs_ref_complement)
    
    print "blocks of runs for reference:"
    print which_runs_ref
    print "complement of reference blocks of runs:"
    print which_runs_ref_complement
    
    ## FIRST ONE IS REFERENCE, DON'T USE BOOTSTRAPS HERE, WILL LOOK AT BOOTSTRAP SAMPLES OF REFERENCE LATER
    run_params1 = RunParameters(resfile_fn=resfile_fn1, phipsi=phipsi, backbone_only=backbone_only, nbins = nbins, permutations=0, adaptive_partitioning=adaptive_partitioning,
                                num_sims=num_sims, num_structs=num_structs, binwidth=options.binwidth, bins=bins, sigalpha=options.sigalpha, which_runs=which_runs1,
                                xvgorpdb=xvgorpdb, xvg_basedir=xvg_basedir1, calc_variance=False, xvg_chidir=options.xvg_chidir, pair_runs=pair_runs_array, skip=options.skip,
                                bootstrap_choose=options.num_sims, calc_mutinf_between_sims=False, load_matrices_numstructs=0, skip_over_steps=options.zoom_to_step, max_num_chis=options.max_num_chis, options=options, bootstrap_set_size=options.num_sims, pdbfile=options.pdbfile, xtcfile=options.xtcfile, blocks_ref = blocks_ref, mutual_divergence="no", output_timeseries=options.output_timeseries, backbone = options.backbone )
       
    ### SET REFERENCE = 1 
    options.reference = 1 #as this is the one with only one bootstrap, bootstrap_choose=num_sims
    
    #only calculate 2nd order terms for target wrt reference
    run_params2 = RunParameters(resfile_fn=resfile_fn2, phipsi=phipsi, backbone_only=backbone_only, nbins = nbins, permutations=0, adaptive_partitioning=adaptive_partitioning,
                                num_sims=num_sims, num_structs=num_structs, binwidth=options.binwidth, bins=bins, sigalpha=options.sigalpha, which_runs=which_runs2,
                                xvgorpdb=xvgorpdb, xvg_basedir=xvg_basedir2, calc_variance=False, xvg_chidir=options.xvg_chidir, pair_runs=pair_runs_array2, skip=options.skip,
                                bootstrap_choose=options.bootstrap_set_size, calc_mutinf_between_sims=False, load_matrices_numstructs=0, skip_over_steps=options.zoom_to_step, max_num_chis=options.max_num_chis, options=options,bootstrap_set_size=options.bootstrap_set_size, pdbfile=options.pdbfile, xtcfile=options.xtcfile,  blocks_ref = blocks_ref, mutual_divergence=options.mutual_divergence, output_timeseries = options.output_timeseries, backbone = options.backbone)  
    #but also use the subsets of the full reference ensemble to look at the variance of the local KL-divergence

    resfile_fn3 = resfile_fn1 #as reference is first one
    which_runs3 = which_runs2 #though with bootstraps just like the target
    resfile_fn4 = resfile_fn2 #need reference for target also for JSDiv
    which_runs4 = which_runs2 #again with bootstraps just like the target
    
    pair_runs_array3 = pair_runs_array2
    run_params3 =  RunParameters(resfile_fn=resfile_fn3, phipsi=phipsi, backbone_only=backbone_only, nbins = nbins, permutations=0, adaptive_partitioning=adaptive_partitioning,
                                num_sims=num_sims, num_structs=num_structs, binwidth=options.binwidth, bins=bins, sigalpha=options.sigalpha, which_runs=which_runs3,
                                xvgorpdb=xvgorpdb, xvg_basedir=xvg_basedir1, calc_variance=False, xvg_chidir=options.xvg_chidir, pair_runs=pair_runs_array3, skip=options.skip,
                                bootstrap_choose=options.bootstrap_set_size, calc_mutinf_between_sims=False, load_matrices_numstructs=0, skip_over_steps=options.zoom_to_step, max_num_chis=options.max_num_chis, options=options,bootstrap_set_size=options.bootstrap_set_size, pdbfile=options.pdbfile, xtcfile=options.xtcfile,  blocks_ref = blocks_ref, mutual_divergence="no", output_timeseries = options.output_timeseries, backbone = options.backbone )  


    pair_runs_array4 = pair_runs_array2
    run_params4 =  RunParameters(resfile_fn=resfile_fn4, phipsi=phipsi, backbone_only=backbone_only, nbins = nbins, permutations=0, adaptive_partitioning=adaptive_partitioning,
                                num_sims=num_sims, num_structs=num_structs, binwidth=options.binwidth, bins=bins, sigalpha=options.sigalpha, which_runs=which_runs4,
                                xvgorpdb=xvgorpdb, xvg_basedir=xvg_basedir1, calc_variance=False, xvg_chidir=options.xvg_chidir, pair_runs=pair_runs_array4, skip=options.skip,
                                bootstrap_choose=options.bootstrap_set_size, calc_mutinf_between_sims=False, load_matrices_numstructs=0, skip_over_steps=options.zoom_to_step, max_num_chis=options.max_num_chis, options=options,bootstrap_set_size=options.bootstrap_set_size, pdbfile=options.pdbfile, xtcfile=options.xtcfile,  blocks_ref = blocks_ref, mutual_divergence="no", output_timeseries="no")  
    
    
    ## LOAD DATA, GET A LIST OF CLASS ResidueChis ##

    if (options.minmax=="yes"):
	    reslist1, reslist2 = load_two_resfiles(run_params1, run_params2, load_angles=True)
    else:
	    reslist1 = load_resfile(run_params1, load_angles=True)
	    reslist2 = load_resfile(run_params2, load_angles=True)

    if(options.bootstrap_set_size < num_sims):  #if we're using statistical filtering
        print "\n loading bootstrap samples of reference for statistics\n"
        reslist3 = load_resfile(run_params3, load_angles=True)
        print "\n calculating KLdiv for reference ensemble\n"
        mykldiv_ref = KLdiv(run_params1, run_params3, reslist1, reslist3)
        mykldiv_ref2 = KLdiv(run_params2, run_params3, reslist2, reslist3)
        print "\n outputting reference kldiv's\n:"
        mykldiv_ref.output()
        ## CALC KLdiv ##
    else:
        if(options.num_sims > 1):
            print "\n using half-sized subsamples and complements of reference for Null Hypothesis distribution for statistics\n"
            mykldiv_ref = KLdiv(run_params1, run_params1, reslist1, reslist1, \
                                which_runs_ref=which_runs_ref, which_runs_ref_complement=which_runs_ref_complement, \
                                do_bootstraps_with_complements = True)
            mykldiv_ref2 = KLdiv(run_params2, run_params2, reslist2, reslist2, \
                                which_runs_ref=which_runs_ref, which_runs_ref_complement=which_runs_ref_complement, \
                                do_bootstraps_with_complements = True)
            #print "\n outputting reference kldiv's:\n"
            #mykldiv_ref.output()
            #mykldiv_ref2.output()
            
    
    mykldiv     = KLdiv(run_params1, run_params2, reslist1, reslist2)
    #myjsdiv2    = KLdiv(run_params2, run_params1, reslist2, reslist1)
    ## Statistical Filtering on mykldiv given reference ##
    #if(options.bootstrap_set_size < num_sims): #if we're using statistical filtering
    if(options.num_sims > 1):
        print "\n filtering KLdiv and JSdiv using Null distribution from reference\n"
        oldkldiv = copy(mykldiv)
        filter_div(mykldiv_ref, mykldiv_ref2, mykldiv)

    ## Output Results ##
    print "\n outputting target kldiv's\n"
    mykldiv.output()
    print "\n outputting pdbs using BioPython\n"
    ## Map onto B-factor of pdbfile
    mykldiv_ref.write_pdbs("_ref")
    mykldiv.write_pdbs()
    #oldkldiv.output_cov()
    
    ## Second Order stuff, if requested ##
    if options.mutual_divergence == "yes":
        if(options.bootstrap_set_size < num_sims):
            mykldiv.do_second_order()
        else: 
            mykldiv.do_second_order()

#######################################################################################################
### test_kldiv: a self-testing routine like a regression test. Basic functionality is tested.
### NEED TO TEST AGAINST ACTUALLY LOADING THE RESLIST THROUGH dihedral_mutent.py TO TEST AGAINST ITS REGRESSION
########################################################################################################

def test_kldiv(test_options, xvg_basedir1, xvg_basedir2, resfile_fn1, resfile_fn2):
    
    
    bins=arange(-180,180,test_options.binwidth) #Compute bin edges
    nbins = len(bins)

    chi_counts1 =   zeros((2, nbins), int64)
    chi_pop_hist1 = zeros((2, nbins),float64)
    chi_counts2 =   zeros((2, nbins), int64)
    chi_pop_hist2 = zeros((2, nbins),float64)
    
    angles = zeros((2,2,test_options.num_sims,12001),float64)
    numangles = 12001
    thisnumangles = [numangles, numangles]
    #read angles from xvg file, map to bins
    kdir = (xvg_basedir1, xvg_basedir2)
    for dataset in range(2): 
      for jchi in range(2):
          for n in range(test_options.num_sims):
                (data, title) = readxvg(str(kdir[dataset])+"run"+str(n+1)+"/chi"+str(jchi+1)+"PHE78.xvg",1,0)
                thisnumangles[dataset] = min(thisnumangles[dataset],len(data[:,1]))
                angles[dataset,jchi,n,:thisnumangles[dataset]] = data[:thisnumangles[dataset],1]
                anglestemp = angles[dataset,jchi,n,:thisnumangles[dataset]]
                anglestemp[anglestemp < 0] += 360 # wrap around
                angles[dataset,jchi,n,:thisnumangles[dataset]] = anglestemp
                if(jchi == 1): #correct PHE78 chi2 for symmetry
                    anglestemp = angles[dataset,jchi,n,:thisnumangles[dataset]]
                    anglestemp[anglestemp > 180] = anglestemp[anglestemp > 180] - 180
                    angles[dataset,jchi,n,:thisnumangles[dataset]] = anglestemp
                for anglenum in range(thisnumangles[dataset]): #map to bins
                    bin_num = binsingle(angles[dataset,jchi,n,anglenum], 1.0 / test_options.binwidth)
                    if(dataset == 0): 
                        chi_counts1[jchi,bin_num] += 1 
                    if(dataset == 1): 
                        chi_counts2[jchi,bin_num] += 1
                print "chi_counts1, chi: "+str(jchi)+ " :" + str(chi_counts1[jchi])
                print "chi_counts2, chi: "+str(jchi)+ " :" + str(chi_counts2[jchi])
         
    #now, calculate kldiv
    print "counts i1:"+str(chi_counts1[0,:])
    print "counts j1:"+str(chi_counts2[0,:])
    print "counts i2:"+str(chi_counts1[1,:])
    print "counts j2:"+str(chi_counts2[1,:])     
    

    for jchi in range(2):
        chi_pop_hist1[jchi,:] = chi_counts1[jchi,:] * 1.0 / (test_options.num_sims * thisnumangles[0])
        chi_pop_hist2[jchi,:] = chi_counts2[jchi,:] * 1.0 / (test_options.num_sims * thisnumangles[1])
    
    kldiv2 = zeros((2),float64)
    jsdiv2 = zeros((2),float64)
    pi1 = chi_pop_hist1[0,:]
    pi2 = chi_pop_hist1[1,:]
    pj1 = chi_pop_hist2[0,:]
    pj2 = chi_pop_hist2[1,:]

    assert(sum(pi1) > 0.99999 and sum(pi1) < 1.00001)
    assert(sum(pi2) > 0.99999 and sum(pi2) < 1.00001)
    assert(sum(pj1) > 0.99999 and sum(pj1) < 1.00001)
    assert(sum(pj2) > 0.99999 and sum(pj2) < 1.00001)
    
    
    print "pi1:"+str(chi_pop_hist1[0,:])
    print "pj1:"+str(chi_pop_hist2[0,:])
    print "pi2:"+str(chi_pop_hist1[1,:])
    print "pj2:"+str(chi_pop_hist2[1,:])     
    
    chi_pop_hist1[chi_pop_hist1==0] = SMALL * 0.5 # to avoid NaNs
    chi_pop_hist2[chi_pop_hist2==0] = SMALL * 0.5 # to avoid NaNs

    
    for mychi in range(2):            
        pi = chi_pop_hist1[mychi,:]
        pj = chi_pop_hist2[mychi,:]
        #Kullback-Leibler Divergence
        kldiv2[mychi] = sum(pj[pi > SMALL] * log (pj[pi > 0 + SMALL]/pi[pi > 0 + SMALL]), axis=-1)
        #Jensen-Shannon Divergence
        jsdiv2[mychi] = 0.5 * sum(pj * log (pj/(0.5*pi+0.5*pj+SMALL)), axis=-1) + 0.5 * \
            sum(pi * log (pi/(0.5*pi + 0.5*pj+SMALL)), axis=-1)
    
    #values to compare: 
    #kldiv chi1: 0.602409  chi2: 0.109821
    #jsdiv chi1: 0.188742  chi2: 0.030788
    kldiv_targ = [0.602408, 0.109821]
    jsdiv_targ = [0.188742, 0.030788]
    print "kldiv     : "+str(kldiv2)+"\n"
    print "kldiv targ: "+str(kldiv_targ)+"\n"
    print "jsdiv     : "+str(jsdiv2)+"\n"
    print "jsdiv targ: "+str(jsdiv_targ)+"\n"

    #assert(kldiv2[0] > 0.602408 and kldiv2[0] < 0.602410)
    #assert(kldiv2[1] > 0.109820 and kldiv2[1] < 0.109822)
    #assert(jsdiv2[0] > 0.188741 and jsdiv2[0] < 0.188743)
    #assert(jsdiv2[1] > 0.030788 and jsdiv2[1] < 0.030788)


###############################################################################
###############################################################################

###############################################################################
###  Main Program: sets options, runs calculation and self-testing routine
###############################################################################


if __name__ == "__main__":
    usage="%prog ref_xvg_basedir1 ref_xvg_basedir2 targ_resfile1 targ_resfile2  # where resfile is in the format <1-based-index> <aa type> <res num>"
    parser=OptionParser(usage)
    parser.add_option("-w", "--binwidth", default=15.0, type="float", help="width of the bins in degrees")
    parser.add_option("-n", "--num_sims", default=None, type="int", help="number of simulations")
    parser.add_option("-d", "--xvg_chidir", default = "/", type ="string", help="subdirectory under xvg_basedir/run# where chi angles are stored")
    parser.add_option("-a", "--minmax", default = "no", type ="string", help="adaptive min,max (yes|no)")
    parser.add_option("-b", "--backbone", default = "phipsichi", type = "string", help="chi: just sc  phipsi: just bb  phipsichi: bb + sc, coarse_phipsi: alpha, beta, turn, coil (4-bin discretization, use only with option -w 90)")
    parser.add_option("-s", "--sigalpha", default=0.1, type="float", help="p-value threshold for statistical filtering, lower is stricter")
    parser.add_option("-i", "--skip", default = 1, type = "int", help="interval between snapshots to consider, in whatever units of time snapshots were output in")
    parser.add_option("-y", "--abs", default = "no", type = "string", help="take absolute value of kldiv?")
    parser.add_option("-g", "--grassberger", default = "no", type = "string", help="use Grassberger estimate of pi and pj?")
    parser.add_option("-z", "--zoom_to_step", default = 0, type = "int", help="skips the first n snapshots in xvg files")
    parser.add_option("-m", "--max_num_chis", default = 99, type = "int", help="max number of sidechain chi angles per residue or ligand")
    parser.add_option("-o", "--bootstrap_set_size", default = None, type = "int", help="perform bootstrapping within this script; value is the size of the subsets to use")
    parser.add_option("-k", "--mutual_divergence", default = "no", type =  "string", help="calculate 2nd-order terms in KL-divergence expansion")
    parser.add_option("-f", "--pdbfile", default = None, type = "string", help="pdb structure file for additional 3-coord cartesian per residue")
    parser.add_option("-q", "--xtcfile", default = None, type = "string", help="gromacs xtc prefix in 'run' subdirectories for additional 3-coord cartesian per residue")
    parser.add_option("-e","--output_timeseries", default = "no", type = "string", help="output corrected dihedral timeseries (requires more memory) yes|no ")
    

    ## SETUP RUN PARAMETERS ##    
    run_params1 = 'None'
    run_params2 = 'None'
    run_params3 = 'None'
    (options,args)=parser.parse_args()
    options.adaptive = "no" # KLdiv first order cannot use adaptive partitioning unless dihedrals were ranked ordered for two systems (i.e. residue lists) together
    
    print "COMMANDS: ", " ".join(sys.argv)
    
    #adaptive_partitioning = (options.adaptive == "yes")  #want "yes" for second-order term
    
    if(options.backbone == "coarse_phipsi"):
        print "overridding binwidth, setting up four bins for backbone discretization"
        options.binwidth = 90.0

    xvg_basedir1, xvg_basedir2, resfile_fn1, resfile_fn2 = args
    run_kldiv(options, xvg_basedir1, xvg_basedir2, resfile_fn1, resfile_fn2)

    ## for self-testing procedure
    test_options=options
    test_options.num_sims = 5
    test_options.binwidth = 30
    test_options.bootstrap_set_size = 5
    test_options.grassberger = "no"
    test_options.xvg_chidir = "/"
    test_options.backbone = "chi"
    test_options.adaptive = "no"
    test_options.skip = 1
    test_options.zoom_to_step = 0
    test_options.mutual_divergence = "no"
    xvg_basedir1_test = str(sys.path[0]+"/test/PHE78_1M47_MD/")
    xvg_basedir2_test = str(sys.path[0]+"/test/PHE78_1M48_MD/")
    resfile_fn1_test = "PHE78.reslist"
    resfile_fn2_test = "PHE78.reslist"

    test_kldiv(test_options, xvg_basedir1_test, xvg_basedir2_test, resfile_fn1_test, resfile_fn2_test)
    

##############################################################################################
## END OF kl_diverge.py
##############################################################################################

