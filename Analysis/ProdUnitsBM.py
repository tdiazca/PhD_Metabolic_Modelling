######
########
##########   This module is for LPA of a response for biomass production while
##########   meeting the cells energy demand (NGAM + GAM ATP demand)   
########    NOTE: remember, there is no NITRATE (NO3) or GLN in MMc media !!!
######              but there is GLN and NO3 in synovial fluid

import sys
from Util import Set
import BuildLP
from Bioinf import Biomass
import SepiBiomass
import DataForAnalysis # so it can recognise CatomsAll and NatomsAll

#DefaultBiomass=SepiBiomass.All

def Tx(m):
    """pre: m = the model
       post: returns a list of transporters"""
    rv=[]
    for r in filter(lambda s: "_tx" in s, m.smx.cnames):
        rv.append(r)
        
    return rv

def BMtx(m):
    """pre: m = the model
       post: returns a list of biomass transporters"""
    rv = []
    for r in filter(lambda s: "_bm_tx" in s, m.smx.cnames):
        rv.append(r)
        
    return rv

## Generate LP for production of whole biomass while meeting
    ## the ATP cell maintenace reqirements

def GetFluxDic(m, PropBF=0.0, AMaint=45.0): # set flux bm_tx for BMP and ATP demand
    """ lp object constrained to biomass fluxes defined by BM composition """
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
       post: generates a flux dictionary for biomass exporters according to description of
       biomass composition and sets flux through ATPase to the ATP maintenace cost"""
    ##bm = SepiBiomass.All.Copy() # Copy() doesnt work!
    bm =  Biomass.Composition(
        1.0,{
            "WholeCell": SepiBiomass.WholeCell,
            "Biofilm" : SepiBiomass.Biofilm
        }
    )
    
    bm.SetAmmount(PropBF, "Biofilm")
    bm.SetAmmount(1.0-PropBF, "WholeCell")
    
    fd = SepiBiomass.FluxDic(m, bm)
    fd["ATPase"] = AMaint

    return fd

def GetLP(m, PropBF=0.0, AMaint=45.0):
    """pre: GetFluxDic(m, PropBF=0.0, AMaint=45.0)     
       post: lp object constrained to biomass fluxes defined by BM composition
       and ATP maintenace demand"""

    lp = BuildLP.BuildLP(m) # lp = LP.lp(m) ; lp = m.GetLP() # lp.SetObjective(m.sm.cnames)
    fd = GetFluxDic(m, PropBF, AMaint)
    for r in fd.keys():
        if 'bm_tx' in r:
            lp.SetFixedFlux(fd)
    lp.SetFluxBounds({"ATPase": (AMaint, None)}) # If we dont want ATPase upper flux to be fixed

    return lp

def GetLPfreeBM(m, PropBF=0.0, AMaint=45.0):
    """pre: GetFluxDic(m, PropBF=0.0, AMaint=45.0)     
       post: lp object constrained to biomass fluxes defined by BM composition
       and ATP maintenace demand"""
    
    orig_rp = dict(m.sm.RevProps)
    for r in filter(lambda s: "bm_tx" in s, m.sm.cnames):
        m.sm.RevProps[r] = m.smx.RevProps[r] = "<>"

    lp = BuildLP.BuildLP(m) # lp = LP.lp(m) ; lp = m.GetLP() # lp.SetObjective(m.sm.cnames)

    fd = GetFluxDic(m, PropBF, AMaint)
        # if PropBF=1, still all other bm_tx are keys in dic, but flux in fd = 0
    
    for r in fd.keys():
        if 'bm_tx' in r:    # otherwise allows ATPase to go to zero flux! (since also in fd)
            lp.SetFluxBounds({r:(None,fd[r])}) # Bug = takes val of '0' if r not declared REV
#	lp.SetFluxBounds({r:(fd[r],None)})  # Bug = takes positive and negative values! if r not declared REV
    m.sm.RevProps = m.smx.RevProps = orig_rp
    
    lp.SetFluxBounds({"ATPase": (AMaint, None)}) # If we dont want ATPase upper flux to be fixed
                    
    return lp

## Analyse system for production of 1 unit of cell biomass (BMP)
    # while meeting ATP cell maintenance requirements

def ProdUnitBM(m, PropBF=0.0, AMaint=45.0, AA=None, o2=None, no3=None,
               targets=["GLN_AA_mm_tx"], FreeBM=None):
    # "GLN_AA_mm_tx", "NO3_mm_tx" are not present in MMc media
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            AA = flux constraint through AA tx
            o2 = flux constraint through O2 tx
            no3 = flux constraint through NO3 tx
            targets = targeted reactions to be blocked
            FreeBM = if != None, upper flux bound of BM_tx (biomass txs) is unbound
       post: returns a dictionary with solution for production of one unit of cell biomass"""

    sol = {}

    if FreeBM != None:
        lp = GetLPfreeBM(m, PropBF, AMaint)
    else:
        lp = GetLP(m, PropBF, AMaint)

    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})

    if o2 != None:
        lp.SetFixedFlux({"O2_mm_tx":o2})

    if no3 != None:
        lp.SetFixedFlux({"NO3_mm_tx":no3})

    for r in targets:
        lp.SetFixedFlux({r:0.0})

    lp.Solve(False)
    sol = lp.GetPrimSol()
    ObjVal = lp.GetObjVal()

    ByProds = {}
    for r in sol:
        if r in filter(lambda s: "_bp_tx" in s, m.smx.cnames):
            ByProds[r] = round(sol[r],3)

    CompUsed = {}
    for r in sol:
        if r in Set.Complement(Tx(m), BMtx(m)):
            CompUsed[r] = round(sol[r],3)

    TotalSolReac = []
    for r in sol:
        TotalSolReac.append(r)

    SolReac = []
    for r in sol:
        if r not in Tx(m):
            SolReac.append(r)

    if FreeBM != None:

        fd=GetFluxDic(m, PropBF, AMaint) # Fixed BM tx fluxes (1 gram BM)

        BMcompsExcreted = {}
        for r in filter(lambda s: "bm_tx" in s, sol.keys()):
            for r in Set.Intersect(sol.keys(),fd.keys()):
                # ATPase and (PIA1_bf_bm_tx, PIA2_bf_bm_tx with flux 0 if PropBF = 0.0)
                if sol[r] != fd[r]:
                    comp = r.split("_")[0]
                    NetExp = sol[r] - fd[r] 
                    BMcompsExcreted[comp+"_NetExport"] = NetExp

    print "Total no of reactions in solution: ", len(TotalSolReac)
    print "Total no of reactions in solution without transporters: ", len(SolReac)
    print "Objective value of the solution in mmol/gDW/h: ", round(ObjVal,3)
    if "GLC_mm_tx" in sol.keys():
        print "This is the ammount of GLC utilised (mmol): ", round(sol["GLC_mm_tx"],3)
    if "NH4_mm_tx" in sol.keys():
        print "This is the ammount of NH4 utilised/excreted (mmol): ", round(sol["NH4_mm_tx"],3)
    print "These are the ammounts of by-products exported (mmol): ", "\n",
    for bp in sorted(ByProds.keys()):
        print bp, ":  ", round(ByProds[bp],3)
    print "These are the amounts of compounds uptaken/excreted (mmol): ", "\n",
    for comp in sorted(CompUsed.keys()):
        print comp, ":  ", round(CompUsed[comp],3)

    if FreeBM != None:
        print "These are the net amounts of BM compounds excreted (mmol): ", "\n",
        for comp in sorted(BMcompsExcreted.keys()):
            print comp, ":  ", round(BMcompsExcreted[comp],6)

    if FreeBM != None:
        return sol, ObjVal, ByProds, CompUsed, TotalSolReac, SolReac, BMcompsExcreted

    else:
        return sol, ObjVal, ByProds, CompUsed, TotalSolReac, SolReac

def ProdBMcharact(m, PropBF=0.0, AMaint=45.0, AA=None, o2=None, no3=None,
                  targets=["GLN_AA_mm_tx"], FreeBM=None, 
                  CatomsAll=DataForAnalysis.Catoms_synovial, NatomsAll=DataForAnalysis.Natoms_synovial):

    """pre: m = the model
            Objective = reaction which flux must be minimised as LP objective
            AMaint = ATP maintenance demand (NGAM plus GAM),growth rate 1gDW/h = 45 mmol/gDW = 45.0
            AA = flux constraint through AA tx
            o2 = flux constraint through O2 tx
            no3 = flux constraint through NO3 tx
            targets = targeted reactions to be blocked
            CatomsAll = dic of n of C atoms in media components
            NatomsAll = dic of n of N atoms in media components
       post: returns a dictionary with solution for production of ATP"""

    res = ProdUnitBM(m, PropBF, AMaint, AA, o2, no3, targets, FreeBM)[0] 

    CompContribCuptake = ContribToTotalCuptake(res, CatomsAll)

    CompContribNuptake = ContribToTotalNuptake(res, NatomsAll)

    return res

def TotalCconsumed(res, CatomsAll=DataForAnalysis.Catoms_synovial):
    """pre: res = LP sol (e.g. ProdUnitBM(m))
            CatomsAll = dic of n of C atoms in media components
       post: returns a dic of total C atoms imported per medium components
            and total C atoms imported"""

    CperComp = {}     # C uptaken
    for comp in CatomsAll.keys():
        if comp in res.keys():
            if res[comp] > 0.0: 
                CperComp[comp] = CatomsAll[comp]*res[comp]
        else:
            CperComp[comp] = 0.0

    TotalC = sum(CperComp.values())

    print "This is the total C consumed (mmol):", round(TotalC,3), "\n",

    return CperComp, TotalC

def ContribToTotalCuptake(res, CatomsAll=DataForAnalysis.Catoms_synovial):
    """pre: res = LP sol (e.g. ProdUnitBM(m))
            CatomsAll = dic of n of C atoms in media components
       post: returns a dic of percentage of contribution of each medium
       component to the total C uptaken"""

    CperComp = TotalCconsumed(res, CatomsAll)[0]

    TotalC = TotalCconsumed(res, CatomsAll)[1]

    propTdicC = {}
    for comp in CatomsAll:
        if not res.has_key(comp): # if comp in medium is not uptaken   
            propTdicC[comp] = 0.0
        else:
            propTdicC[comp] = CperComp[comp]/TotalC
            
    TotalpercentC = round((sum(propTdicC.values())*100),3)

    print "\n", "These are the contributions of medium components to total C uptake (percentage): ", "\n"
    for comp in sorted(CatomsAll.keys()):
        print comp, "","","","", round(propTdicC[comp]*100,3)

    print "\n", "Sum: ", TotalpercentC, "\n"

    return propTdicC 

def TotalNconsumed(res, NatomsAll=DataForAnalysis.Natoms_synovial): # NatomsAll=NatomsAll
    """pre: res = LP sol (e.g. ProdUnitBM(m))
            NatomsAll = dic of n of N atoms in media components
       post: returns a dic of total N atoms imported per medium components
            and total N atoms imported"""

    NperComp = {}     # N uptaken
    for comp in NatomsAll.keys():
        if comp in res.keys():
            if res[comp] > 0.0: 
                NperComp[comp] = NatomsAll[comp]*res[comp]
        else:
            NperComp[comp] = 0.0

    TotalN = sum(NperComp.values())

    print "These is the total N consumed (mmol):", round(TotalN,3), "\n",

    return NperComp, TotalN

def ContribToTotalNuptake(res, NatomsAll=DataForAnalysis.Natoms_synovial): #NatomsAll=NatomsAll
    """pre: res = LP sol (e.g. ProdUnitBM(m))
            NatomsAll = dic of n of N atoms in media components
       post: returns a dic of percentage of contribution of each medium
       component to the total N uptaken"""

    NperComp = TotalNconsumed(res, NatomsAll)[0]

    TotalN = TotalNconsumed(res, NatomsAll)[1]

    propTdicN = {}

    for comp in NatomsAll:
        if not res.has_key(comp): # if comp in medium is not uptaken   
            propTdicN[comp] = 0.0
        else:
            if res[comp] > 0.0: 
                propTdicN[comp] = NperComp[comp]/TotalN
            
    TotalpercentN = round((sum(propTdicN.values())*100),3)

    print "\n", "These are the contributions of medium components to total N uptake (percentage): ", "\n"
    for comp in sorted(NatomsAll.keys()):
        if comp in propTdicN.keys():
            print comp, "","","","", round(propTdicN[comp]*100,3)
    print "\n", "Sum: ", TotalpercentN, "\n"

    return propTdicN

#################################################################
