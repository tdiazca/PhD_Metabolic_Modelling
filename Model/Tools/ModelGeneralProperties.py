
######
########
##########   This module is used to define the general properties of the model
##########    
########        
######

import ScrumPy
import BuildLP
import SanityChecks
import BuildSepi
reload(BuildSepi)
from Util import Set
from Bioinf.PyoCyc import Compound

db=BuildSepi.UsrBuildOrg.orgdb

def GetTx(m): # move to a general module

    """pre: m = the model
       post: returns a list of transporters"""

    rv = filter(lambda s: "_tx" in s, m.smx.cnames)

    return rv

def BMtx(m): # move to a general module

    """pre: m = the model
       post: returns a list of biomass transporters"""

    rv = filter(lambda s: "_bm_tx" in s, m.smx.cnames)

    return rv

## Gene - reaction association

def RemoveSplit(reac):
    
    """pre: r = a reaction
       post: returns reaction name after removing -(NADP)/-(NAD)"""
    
    r = reac.rsplit('-(NADP)')[0].rsplit('-(NAD)')[0] #reverse split
    
    return r

def Reactions(m):

    """pre: True
       post: returns list of reactions excluding transporters"""
    
    rv = []
    
    tx = filter(lambda s: "_tx" in s, m.sm.cnames)
    
    for reac in Set.Complement(m.sm.cnames, tx):
            r = RemoveSplit(reac) # removes -NAD/NADP
            rv.append(r)

    return rv

def ReacAsso(m, reactions=[]):
    
    """pre: True
       post: returns dictionary with key as all reacs in model excluding transporters
              and values are lists of associated genes
              "No Gene" if reaction is not associated to genes
              "Reaction UID not found" if reac not present in database
               Note: reac names in keys would have -filter(lambda s: "_tx" in s, m.sm.cnames)NAD/NADP split from name"""
        
    rv = {}
        
    if reactions ==[]:
        reactions = Set.Complement(m.sm.cnames, GetTx(m)) # for all reacs except transporters

    for reac in reactions:
        r = RemoveSplit(reac)   # removes -NAD/NADP
        if db.has_key(r):
            genes=[]
            genes = db[r].GetGenes()
            if len(genes) != 0:
                rv[r] = genes
            else:
                rv[r] = "No genes"
        else:
            rv[r] = "Reaction UID not found"

    return rv

def Genes(m, reactions=[]):

    """pre: True
       post: returns list of genes associated with reactions in the model"""

    rv = []
    ra = ReacAsso(m, reactions)
    for r in ra.keys():
        if type(ra[r]) is list:
            for g in ra[r]:
                if g not in rv:
                    rv.append(g)
    return rv

def ReacToGene(m, reactions=[], reacasso={}):
    
    """pre: True
       post: returns dictionary with reacs with associated genes as keys
             and list of associated genes as values"""
        
    rv = {}
        
    if reacasso == {}:  # we could parse any other dict of reacs and genes if we wanted
        allreac = ReacAsso(m,reactions) # dict of reacs and their genes
        allreactions = m.sm.cnames
    else:
        #allreac = ReacAsso(m,reacasso.keys())
        allreac = ReacAsso(m,Set.Complement(reacasso.keys(), GetTx(m)))
        allreactions = reacasso.keys()

    NoGene = []
    NoUID = []
        
    for reac in Set.Complement(allreactions, GetTx(m)):
        r = RemoveSplit(reac)
        if r in allreac.keys(): # remember reacs names here have -NADP/NAD split!
            if allreac[r] == "No genes":
                NoGene.append(reac) # so add full reac name to list (before split occured)
            elif allreac[r] == "Reaction UID not found":
                NoUID.append(reac)
            else:
                rv[reac] = allreac[r]
                    
    print "Total no of reactions with gene association:",len(rv.keys())
    print "Total no of reactions with no gene association:",len(NoGene)
    print "Total no of reactions with no UID in db:",len(NoUID)

    return rv, NoGene, NoUID # make the dict and list available

## General model properties

def ModelProperty(m):

    """pre: m = model
       post: returns lists of total reactions, reactions excluding tx,
       transporters, metabolites, orphan metabolites, dead reactions,
       and reactions subsets"""

    TotalReacs = m.sm.cnames
    Transporters = filter(lambda s: "_tx" in s, m.smx.cnames)   # list
    Reactions = Set.Complement(m.sm.cnames, Transporters)       # list
    Metabolites = m.sm.rnames                                   # list
    OrphanMets = m.sm.OrphanMets()                              # list
    DeadReacs = Set.Complement(m.DeadReactions(), Transporters) # list
    ReacSubsets = reacsubs = m.EnzSubsets()                     # class

    print 'Total no of reactions including transporters:', len(TotalReacs)
    print 'Total no of transporters:', len(Transporters)
    print 'Total no of reactions excluding transporters:', len(Reactions)
    print 'Total no of internal metabolites:', len(Metabolites)
    print 'Total no of orphan metabolites:', len(OrphanMets)
    print 'Total no of dead reactions:', len(DeadReacs)
    print 'Total no of reaction subsets:', len(ReacSubsets) 

## ElMos of the ETC (analysis ETC stand alone module)

def ImportETCsubmodule(m):
    """pre: m = the model
       post: returns a submodel containing reacs involved in the ETC"""

    m1=ScrumPy.Model('ETC_TopLevel.spy') ; m1.Hide()
    return m1

def ElementaryModesETC(m):  #m = ScrumPy.Model('ETC_TopLevel.spy')

    """pre: ImportETCsubmodule(m)
            m1 = stand-alone module of the ETC ('ETC_TopLevel.spy')
       post: returns the net stoichiometries of the elementary modes of the ETC"""

    #k=m.sm.NullSpace()
    #k.IntegiseC()          # show entries in k matrix as integers
    #print k
    m1=ImportETCsubmodule(m); m1.Hide()
    elmo = m1.ElModes()
    elmo.Integise()         # show elmo coeficients as integers

    print 'These are the net stoichiometries of the elementary modes of the ETC stand-alone module:',"\n","\n",elmo.Stos()
    
    return elmo

    #print elmo.ReacsOf('ElMo_0')

## ElMos of a submodel

def ImportSubmodule(m):
    """pre: m = the model
       post: returns a submodel"""

    m1=ScrumPy.Model(m) ; m1.Hide()
    return m1

def ElementaryModesSubmodel(m):  #m = ScrumPy.Model('Name_submodelfile.spy')

    """pre: ImportSubmodule(m)
            m1 = submodel ('Name_submodelfile.spy')
       post: returns the net stoichiometries of the elementary modes of the submodel"""

    #k=m.sm.NullSpace()
    #k.IntegiseC()          # show entries in k matrix as integers
    #print k
    m1=ImportSubmodule(m); m1.Hide()
    elmo = m1.ElModes()
    elmo.Integise()         # show elmo coeficients as integers

    print 'These are the net stoichiometries of the elementary modes of the submodel:',"\n","\n",elmo.Stos()
    
    return elmo

    #print elmo.ReacsOf('ElMo_0')
