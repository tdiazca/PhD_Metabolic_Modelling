######
########
##########   This module is for LPA of by-products utilisation for ATP production
##########    # For this, remember to block GLC import in function to study effect of GLC depletion
########    NOTE: remember, there is no NITRATE (NO3) in MMc media !!!
######

import sys
#sys.path.append('/home/stffmh2/eqd15ypu/Sepi_RP26A/Analysis/ByProd/')
from Util import Set
import BuildLP      # Need to load model object first (gives access to Tools dir)
reload(BuildLP)
##import ScrumPy    # Not needed if m is imported as a model object first.
##reload(ScrumPy)
#import LP          # Not needed if we import BuildLP (it imports it there)

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

def ATPprodMinFlux_ACET(m, AA=0.0, o2=None, no3=None, blocks=["GLN_AA_mm_tx","GLC_mm_tx"]):
    # Block "GLC_mm_tx" to study ACET util as C 
        # GLN is not in MMc media! ;  # No no3 in MMc media ;                                                                   
    """pre: m = the model
            AA = flux constraint through AA tx
            o2 = flux constraint through O2 tx
            no3 = flux constraint through NO3 tx
            blocks = reactions to be blocked
       post: returns a dictionary with solution for production of one mmol of ATP with
       minimisation of total flux as objective function and ACET available as C source"""

    sol = {}

    orig_rp = dict(m.sm.RevProps)

    #orig_rp_acet_tx = m.sm.RevProps["ACET_bp_tx"] # store original reversibility or "ACET_bp_tx"

    m.sm.RevProps["ACET_bp_tx"] = m.smx.RevProps["ACET_bp_tx"] = "<>" 

    lp = BuildLP.BuildLP(m)
    lp.SetFixedFlux({"ATPase":1})

    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})

    if o2 != None:
        lp.SetFixedFlux({"O2_mm_tx":o2})

    if no3 != None:
        lp.SetFixedFlux({"NO3_mm_tx":no3})

    for r in blocks:
            lp.SetFixedFlux({r:0.0})

    lp.Solve()

    m.sm.RevProps = m.smx.RevProps = orig_rp    # revert changes in directionallity of reacs

    #m.sm.RevProps = m.smx.RevProps = orig_rp_acet_tx # revert changes in directionallity of "ACET_bp_tx"

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

    if "ATPSynth" and "O2_mm_tx" in sol.keys():
        PO = sol["ATPSynth"] / (sol["O2_mm_tx"]*2) # Oxidative phosphorilation
    else:
        PO = 0.0

    if "ATPase" and "O2_mm_tx" in sol.keys():
        FullPO = sol["ATPase"] / (sol["O2_mm_tx"]*2)  # Oxidative plus substrate level phosphorilation
    else:
        FullPO = 0.0

    if "ACET_bp_tx" in sol.keys():
        ATPacet = sol["ATPase"] / sol["ACET_bp_tx"]  # ATP/ACET ratio
    else:
        ATPacet = 0.0
        
    if "GLC_mm_tx" in sol.keys():
        ATPglc = sol["ATPase"] / sol["GLC_mm_tx"]  # ATP/GLC ratio
    else:
        ATPglc = 0.0

    print "Total no of reactions in solution: ", len(TotalSolReac)
    print "Total no of reactions in solution without transporters: ", len(SolReac)
    if "ATPase" and "O2_mm_tx" in sol.keys():
        print "This is the P/O ratio (total): ", round(FullPO,3) 
    if "ATPSynth" and "O2_mm_tx" in sol.keys():
        print "This is the P/O ratio for oxidative phosphorilation: ", round(PO,3)
    if "ACET_bp_tx" in sol.keys():
        print "This is the ATP/ACET ratio: ", round(ATPacet,3)
    if "GLC_mm_tx" in sol.keys():
        print "This is the ATP/GLC ratio: ", round(ATPglc,3)
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

    return sol, ObjVal, ByProds, CompUsed, TotalSolReac, SolReac, PO, FullPO, ATPacet, ATPglc

def ATPprodMinFlux_ByProd(m, AA=0.0, o2=None, no3=None, blocks=["GLN_AA_mm_tx"]):
    # Block "GLC_mm_tx" to study GLC depletion conditions;
        # GLN is not in MMc media! ;  # No no3 in MMc media ;                                                                   
    """pre: m = the model
            AA = flux constraint through AA tx
            o2 = flux constraint through O2 tx
            no3 = flux constraint through NO3 tx
            blocks = reactions to be blocked
       post: returns a dictionary with solution for production of one mmol of ATP with
       minimisation of total flux as objective function and by-prods available as C sources"""

    sol = {}

    orig_rp = dict(m.sm.RevProps)

    #orig_rp_acet_tx = m.sm.RevProps["ACET_bp_tx"] # store original reversibility or "ACET_bp_tx"

    for r in filter(lambda s: "_bp_tx" in s, m.sm.cnames):
        m.sm.RevProps[r] = m.smx.RevProps[r] = "<>" 

    lp = BuildLP.BuildLP(m)
    lp.SetFixedFlux({"ATPase":1})

    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})

    if o2 != None:
        lp.SetFixedFlux({"O2_mm_tx":o2})

    if no3 != None:
        lp.SetFixedFlux({"NO3_mm_tx":no3})

    for r in blocks:
            lp.SetFixedFlux({r:0.0})

    lp.Solve()

    m.sm.RevProps = m.smx.RevProps = orig_rp    # revert changes in directionallity of reacs

    #m.sm.RevProps = m.smx.RevProps = orig_rp_acet_tx # revert changes in directionallity of "ACET_bp_tx"

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

    if "ATPSynth" and "O2_mm_tx" in sol.keys():
        PO = sol["ATPSynth"] / (sol["O2_mm_tx"]*2) # Oxidative phosphorilation
    else:
        PO = 0.0

    if "ATPase" and "O2_mm_tx" in sol.keys():
        FullPO = sol["ATPase"] / (sol["O2_mm_tx"]*2)  # Oxidative plus substrate level phosphorilation
    else:
        FullPO = 0.0

    if "ACET_bp_tx" in sol.keys():
        ATPacet = sol["ATPase"] / sol["ACET_bp_tx"]  # ATP/ACET ratio
    else:
        ATPacet = 0.0
        
    if "GLC_mm_tx" in sol.keys():
        ATPglc = sol["ATPase"] / sol["GLC_mm_tx"]  # ATP/GLC ratio
    else:
        ATPglc = 0.0

    print "Total no of reactions in solution: ", len(TotalSolReac)
    print "Total no of reactions in solution without transporters: ", len(SolReac)
    if "ATPase" and "O2_mm_tx" in sol.keys():
        print "This is the P/O ratio (total): ", round(FullPO,3) 
    if "ATPSynth" and "O2_mm_tx" in sol.keys():
        print "This is the P/O ratio for oxidative phosphorilation: ", round(PO,3)
    if "ACET_bp_tx" in sol.keys():
        print "This is the ATP/ACET ratio: ", round(ATPacet,3)
    if "GLC_mm_tx" in sol.keys():
        print "This is the ATP/GLC ratio: ", round(ATPglc,3)
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

    return sol, ObjVal, ByProds, CompUsed, TotalSolReac, SolReac, PO, FullPO, ATPacet, ATPglc
