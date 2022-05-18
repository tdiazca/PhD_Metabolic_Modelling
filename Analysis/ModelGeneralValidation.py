
######
########
##########   This module is used to perform theoretical validation of the model
##########  
########        
######                

import BuildLP
import SanityChecks
import BuildSepi
reload(BuildSepi)
from Util import Set
from Bioinf.PyoCyc import Compound

db=BuildSepi.UsrBuildOrg.orgdb

## E conservation

def ATPCons(m): # lp object checks m for energy conservation

    """pre: m = the model
       post: if energy is not conseved, returns a solution for
       producing ATP out of nothing"""

    lp = BuildLP.BuildClosedLP(m)   # constraints tx flux to zero
    lp.SetFixedFlux({"ATPase":1}) # set flux through ATPase r to 1
    lp.Solve(False)
    
    if lp.IsStatusOptimal():
        sol=lp.GetPrimSol()
        return sol
    
    else:
        print "ATP is not produced when medium components are not available"

## RedOx conservation

def NADHCons(m): # lp object checks m for RedOx conservation

    """pre: NADHOxid reaction must be included in m (uncomment in top module)
            m = the model
       post: if RedOx balance is not conserved returns solution"""

    lp = BuildLP.BuildClosedLP(m)   # creates new LP, min flux as objective. Constraints tx flux to zero
    lp.SetFixedFlux({"NADHOxid":1}) # set flux through ATPase r to 1
    lp.Solve(False)
    
    if lp.IsStatusOptimal():
        sol=lp.GetPrimSol()
        return sol
    
    else:
        print  "NADH is not produced when medium components are not available"

def NADPHCons(m):

    """pre: reac NADPHOxid reaction must be included in m (uncomment in top module)
            m = the model
       post: if RedOx balance is not conserved returns solution"""

    lp = BuildLP.BuildClosedLP(m)   # constraints tx flux to zero
    lp.SetFixedFlux({"NADPHOxid":1}) # set flux through ATPase r to 1
    lp.Solve(False)
    
    if lp.IsStatusOptimal():
        sol=lp.GetPrimSol()
        return sol
    
    else:
        print  "NADPH is not produced when medium components are not available"

## Mass conservation

def BMtx(m): # move to a general module

    """pre: m = the model
       post: returns a list of biomass transporters"""
    
    rv = []

    for r in filter(lambda s: "_bm_tx" in s, m.smx.cnames):
        rv.append(r)

    return rv

def MassCons(m): # lp object checks m for mass conservation

    """pre: m = the modelfilter(lambda s: "_tx" in s, m.sm.cnames)
       post: returns a dictionary with keys as transporters for biomass components
       and values as solution for their production out of 'nothing' if mass is not
       conserved, or value {} if conserved"""

    MassDic={}

    lp = BuildLP.BuildClosedLP(m)   # constraints tx flux to zero

    for bmtx in BMtx(m):
        lp.ClearFluxConstraint(bmtx)
        lp.SetFixedFlux({bmtx:-1}) # set flux through BM tx to -1
        lp.Solve(False)
        
        sol = lp.GetPrimSol()
        MassDic[bmtx] = sol
        
        lp.SetFixedFlux({bmtx:0.0})

    UnconComp = []

    for comp in MassDic.keys():
        if MassDic[comp] != {}:
            UnconComp.append(comp)

    print "These compounds can be produced out of 'nothing': ", "\n", UnconComp

    return MassDic, UnconComp

## Atomic balance

def AtomicCheck(m):

    """pre: m is the model
       post: returns a dictionary with atomically unbalanced reactions as keys and
             their atomic imbalance as values. Returns a list of compounds with
             unknown empirical formulae"""

    # get PGDB used to build m (db=BuildSepi.orgdb) done out of fx (top of module)
    # get compound info from PGDB

    extras = Compound.DB("./","ExtraCompounds.dat") # add extra info on Unknown compounds (emp formulae)

    Unbals = SanityChecks.CheckImbals(m,db,extras)
    Unknowns = SanityChecks.FindUnknownCompounds(Unbals)

    print "This is the number of unbalanced reactions:", len(Unbals)
    print "This is the unbalance associated with the reactions:"
    for r in Unbals:
        print r,"\n", Unbals[r]
    print "This is the number of compounds with undefined empirical formulae:", len(Unknowns)
    print "These are the unknown compounds:", "\n", Unknowns
 
## Stoichiometry consistency

def StoiCons(m):

    """pre: m is the model
       post: returns a list of unconserved metabolites"""
    
    from Structural import StoCons
    
    Nt = m.smx.Copy(tx=1)                      # get transposed copy of external stoichiometry matrix (reacs in rows)
    UnconMets = StoCons.UnconservedMets(Nt)    # list of unconserved mets from conservation relationships

    print "This is the no of unconserved metabolites:", len(UnconMets)
    print "These metabolites are unconserved:",UnconMets
