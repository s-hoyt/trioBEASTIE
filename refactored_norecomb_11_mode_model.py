#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE 
# MIT License
# Copyright (C)2022 William H. Majoros (bmajoros@alumni.duke.edu)
# and Copyright (C)2025 Stephanie H. Hoyt (stephanie.hoyt@duke.edu)
#=========================================================================
import sys
import ProgramName
import getopt
from Rex import Rex
rex=Rex()
from EssexParser import EssexParser
import numpy as np
import subprocess
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import ListVector
##these lines only need to be run once to set everything up in the python environment
# utils = importr('utils')
# utils.install_packages('rstan')
# utils.install_packages('codetools')
##
rstan = importr('rstan')
from math import exp
from scipy.special import logsumexp

MODES_nums = {1: [0, 0, 0, 0, 0, 0], #Null
              2: [0, 0, 0, 0, 1, 0], #Denovo
              3: [0, 0, 0, 0, 0, 1], #Denovo
              4: [0, 1, 0, 0, 0, 0], #Mother affected
              6: [0, 0, 0, 1, 0, 0], #Father affected
              8: [1, 0, 0, 0, 1, 0], #Mother affected, child inherits
              10: [0, 0, 1, 0, 0, 1]} #Father affected, child inherits

MODES = {1: "00 00 00 = all unaffected",
         2: "00 00 10 = child has a de novo in the causal variant",
         3: "00 00 01 = child has a de novo in the causal variant",
         4: "01 00 00 = mother affected, child doesn't inherit",
         6: "00 01 00 = father affected, child doesn't inherit",
         8: "10 00 10 = mother affected, child inherits",
         10: "00 10 01 = father affected, child inherits"}
#=========================================================================
class Site:
    def __init__(self,ID,phased):
        self.ID=ID
        self.phased=phased
        self.counts=np.zeros((3,2),int) # [indiv][haplotype]
        self.het=[0]*3 # [indiv]
#=========================================================================
class Gene:
    def __init__(self,ID):
        self.ID=ID
        self.sites=[]
    def addSite(self,site):
        self.sites.append(site)
#=========================================================================
def parseGene(sxGene):
    ID=sxGene[0];
    gene=Gene(ID);
    numSites=sxGene.numElements()-1
    for i in range(numSites):
        sxSite=sxGene[i+1];
        ID=sxSite[0]
        sxGenotypes=sxSite.findChild("genotypes")
        sxCounts=sxSite.findChild("counts")
        isPhased=sxSite.getAttribute("phased")
        site=Site(ID,isPhased)
        site.het[0]=isHet(sxGenotypes,"mother")
        site.het[1]=isHet(sxGenotypes,"father")
        site.het[2]=isHet(sxGenotypes,"child")
        site.counts[0]=getCounts(sxCounts,"mother")
        site.counts[1]=getCounts(sxCounts,"father")
        site.counts[2]=getCounts(sxCounts,"child")
        gene.addSite(site)
    return gene

def getCounts(sxCounts,label):
    child=sxCounts.findChild(label)
    return np.array([child[0],child[1]])

def isHet(sxGenotypes,label):
    child=sxGenotypes.findChild(label)
    return child[0]!=child[1]

def digit(c):
    return ord(c)-ord('0')

def initModes():
    array3D=[] # indexed as: [mode,indiv,haplotype]
    for i in MODES:
        modeString=MODES[i]
        # "00 00 10 = child has a de novo in the causal variant"
        rec=[] # indexed as: [indiv,haplotype]
        rec.append([digit(modeString[0]),digit(modeString[1])])
        rec.append([digit(modeString[3]),digit(modeString[4])])
        rec.append([digit(modeString[6]),digit(modeString[7])])
        array3D.append(rec)
    return array3D

def getHetsR(gene):
    hetMom = []
    hetDad = []
    hetChild = []
    for site in gene.sites:
        hetMom.append(site.het[0])
        hetDad.append(site.het[1])
        hetChild.append(site.het[2])
    hetsR = robjects.BoolVector(hetMom + hetDad + hetChild)
    hetsMatrixR = robjects.r.matrix(hetsR, ncol = 3, byrow = False) 
    return hetsMatrixR

def getCountsR(gene, numSites):
    counts = [[[], [], []], [[], [], []]]
    for site in gene.sites:
        counts[0][0].append(int(site.counts[0][0])) #mother allele 1
        counts[1][0].append(int(site.counts[0][1])) #mother allele 2
        counts[0][1].append(int(site.counts[1][0])) #father allele 1
        counts[1][1].append(int(site.counts[1][1])) #father allele 2
        counts[0][2].append(int(site.counts[2][0])) #child allele 1
        counts[1][2].append(int(site.counts[2][1])) #child allele 2
    countsR = robjects.IntVector(counts[0][0] + counts[0][1] + counts[0][2] + counts[1][0] + counts[1][1] + counts[1][2])
    countsMatrixR = robjects.r.array(countsR, dim = [numSites, 3, 2])
    return countsMatrixR

def getPhasingR(gene, numSites):
    phasing = []
    for site in gene.sites:
        phasing.append(int(site.phased))
    #casting it as a vector of dim = numSites matters when the bool vector is of length 1 (ie only 1 site), so it doesn't get interpreted as a scalar
    phasingR = robjects.r.array(robjects.BoolVector(phasing), dim = numSites) 
    return phasingR

def runNull(hets, counts, phasing, model, numSites, probAffected, probRecomb, probDenovo):
    nullMode = robjects.r.matrix(robjects.IntVector(MODES_nums[1]), ncol = 2)
    data = {"N_SITES": numSites, 
            "mode": nullMode, 
            "het": hets,
            "count": counts, 
            "isPhased": phasing,
            "probAffected": probAffected,
            "probRecomb": probRecomb,
            "probDenovo": probDenovo}
    named_list = ListVector(data)
    fitNull = rstan.sampling(object = model, 
                                data = named_list,
                                iter = 1,
                                chains = 1,
                                algorithm = "Fixed_param")
    null_posterior = rstan.get_posterior_mean(fitNull)[6] #accessed by index of numerator here instead of by name as in Rscript
    return null_posterior

def runAlt(hets, counts, phasing, model, numSites, probAffected, probRecomb, probDenovo, numSamples, null_posterior, outFile):
    theta_values = {1: 1.0} #initialize these w/ fixed point nulls and results from Null model
    theta_var_values = {1: 0.0}
    numerator_values = {1: null_posterior}
    rhat_values = {1: "NA"}
    for i in MODES_nums:
        if i == 1: continue #null mode done in another function
        altMode = robjects.r.matrix(robjects.IntVector(MODES_nums[i]), ncol = 2, byrow = True)
        data = {"N_SITES": numSites, 
                "mode": altMode, 
                "het": hets,
                "count": counts, 
                "isPhased": phasing,
                "probAffected": probAffected,
                "probRecomb": probRecomb,
                "probDenovo": probDenovo}
        named_list = ListVector(data)
        fit = rstan.sampling(object = model, 
                                data = named_list,
                                iter = numSamples,
                                init = 1)
        #check Rhat values
        num_params = len(robjects.r.names(fit))
        rhat_col_num = 9
        #number of params(rows - can change as we add/remove things from gen quant) * column number of rhat (9)
        rhat = robjects.r.summary(fit)[0][num_params*rhat_col_num] 
        rhat_values[i] = rhat
        if rhat > 1.05:
            with open(outFile, "a") as f:
                f.write("ERROR: DID NOT CONVERGE !!! MODE " + str(i+1) + "\n")
        #save results for this model
        #for theta, want to index into last column (avg'd over all chains)
        theta_col_num = 4
        theta_values[i] = rstan.get_posterior_mean(fit)[num_params * theta_col_num]
        theta_var_values[i] = robjects.r.var(rstan.extract(fit, "theta")[0])[0]
        # as long as num is the first thing in the gen quant block, and none of the other params change, 7 is correct index. otherwise...
        numerator = robjects.r.summary(fit)[0][7] #getting numerator aka posterior (likelihood * priors).
        numerator_values[i] = numerator
    return(theta_values, theta_var_values, numerator_values, rhat_values)
        
#=========================================================================
# main()
#=========================================================================

(options,args) = getopt.getopt(sys.argv[1:], "c:")
if(len(args)!=8):
    exit(ProgramName.get()+"[-c continue] <model> <input.essex> <#MCMC-samples> <firstGene-lastGene> <P(affected)> <P(recomb)> <P(denovo)> <outFile>\n  gene range is zero-based and inclusive\n")
(model,inputFile,numSamples,geneRange,probAffected,probRecomb,probDenovo,outFile)=args

numSamples = int(numSamples)
probAffected = float(probAffected)
probRecomb = float(probRecomb)
probDenovo = float(probDenovo)

if(not rex.find(r"(\d+)-(\d+)",geneRange)):
    exit(geneRange+": specify range of gene: first-last")
firstIndex=int(rex[1])
lastIndex=int(rex[2])

#If continuation = False, start at firstIndex and print header lines
#else, detect the last gene printed to the file and continue from there (no header lines)
continuation=False
for pair in options:
    (key,value)=pair
    if(key=="-c"): continuation=value

#read last gene ID and reset firstIndex
if continuation:
    with open(outFile, 'r') as f:
        last_line = f.readlines()[-12].strip() #changed from [-1] which reads actual last line to [-12] which reads last gene header if full gene output was printed
        print(last_line)
        if(not rex.find(r"GENE",last_line)):
            exit("Couldn't continue; invalid last line: " + last_line)
        firstIndex = int(last_line.split("GENE")[1]) + 1 #changed to be plus one because we already have output for this gene, so want to do the next one

modeArray=initModes()

#compile models once, instead of per gene or per alt model
null_m = rstan.stan_model(model + "_null.stan")
alt_m = rstan.stan_model(model + "_alt.stan")

# Process each gene
#geneIndex refers to position in essex file, not the labeled GENEID
geneIndex=0
parser=EssexParser(inputFile)

if not continuation:
    with open(outFile, "w") as f:
        f.write("Gene\tposteriorProb\tMode\ttheta\ttheta_var\tRhat\tModeDescrip\n")
while(True):
    elem=parser.nextElem()
    if(elem is None): break
    if(elem.getTag()!="gene"): raise Exception("Expecting 'gene' tag in essex")
    if(geneIndex<firstIndex):
        geneIndex+=1
        continue
    elif(geneIndex>lastIndex): break
    gene=parseGene(elem)
    if(gene is None): continue
    if not (continuation and geneIndex == firstIndex): #don't print the header for the first gene after continuing ; already printed
        with open(outFile, "a") as f:
            f.write(gene.ID + "\n")

    hetsR = getHetsR(gene)
    countsR = getCountsR(gene, len(gene.sites))
    phasingR = getPhasingR(gene, len(gene.sites))

    null_posterior = runNull(hetsR, countsR, phasingR, 
                                null_m,
                                len(gene.sites),
                                probAffected, 
                                probRecomb,
                                probDenovo)
        
    (theta_values, 
        theta_var_values,
        numerator_values, 
        rhat_values) = runAlt(hetsR, countsR, phasingR,
                                alt_m,
                                len(gene.sites),
                                probAffected, 
                                probRecomb,
                                probDenovo,
                                numSamples,
                                null_posterior,
                                outFile)
    
    
    #get posterior probs using bayes thm
    posterior_probs = {}
    denom = logsumexp(list(numerator_values.values()))
    for i in MODES_nums:
        posterior_probs[i] = exp(numerator_values[i] - denom)
    #sort based on which are most likely, relative to null
    sorted_posterior = dict(sorted(posterior_probs.items(), key=lambda item: item[1], reverse = True))
    #write to file
    for i in range(len(sorted_posterior)):
        m = list(sorted_posterior.keys())[i]
        this_rhat = rhat_values[m] if type(rhat_values[m]) is str else str(round(rhat_values[m], 3))
        with open(outFile, "a") as f:
            f.write("\t" + str(round(posterior_probs[m] * 100, 2)) + "%" + "\t" +\
                        "MODE " + str(m) + "\t" +\
                        str(round(theta_values[m], 2)) + "\t" +\
                        str(round(theta_var_values[m], 3)) + "\t" + \
                        this_rhat + "\t" +\
                        MODES[m] + "\n")
    geneIndex+=1