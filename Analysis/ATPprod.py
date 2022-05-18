######
########
##########   This module is for LPA of ATP production under different conditions
##########   without considering BMP (theoretical maximum ATP yields)   
########    NOTE: remember, there is no NITRATE (NO3) or GLN (glutamine) in MMc medium.
######              BUT both are present in synovial fluid. Block their tx accordingly.

from Util import Set
import BuildLP
import DataForAnalysis # so it can recognise CatomsAll and NatomsAll

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

## Analyse system for production of ATP maintenace cost (45.0 mmol) or 1 mol ATP
    ## in presence/ absence O2 and NO3

def ProdATP(m, Objective=None, ATPdemand=1.0, AA=None, o2=None, no3=None,
            blocks=[]):
    # "GLN_AA_mm_tx", "NO3_mm_tx" are not present in MMc media but are present in synovial fluid

    """pre: m = the model
            Objective = reaction which flux must be minimised as LP objective
            ATPdemand = if ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            AA = flux constraint through AA tx
            o2 = flux constraint through O2 tx
            no3 = flux constraint through NO3 tx
            blocks = reactions to be blocked
       post: returns a dictionary with solution for production of ATP"""

    sol = {}

    if Objective!= None:
        lp = BuildLP.BuildMinObjective(m, Objective) # lp = m.GetLP(); # lp.SetObjective(Objective)
    else:
        lp = BuildLP.BuildLP(m) # lp = m.GetLP(); # lp.SetObjective(m.sm.cnames)

    lp.SetFixedFlux({"ATPase": ATPdemand})

    if AA != None:
        AAtx=filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames)
        if Objective!= None:
            for aa in Set.Complement(AAtx, Objective):
                lp.SetFixedFlux({aa:AA})
        else:
            for aa in AAtx:
                lp.SetFixedFlux({aa:AA})

    if o2 != None:
        lp.SetFixedFlux({"O2_mm_tx":o2})

    if no3 != None:
        lp.SetFixedFlux({"NO3_mm_tx":no3})

    for r in blocks:
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

    return sol, ObjVal, ByProds, CompUsed, TotalSolReac, SolReac

def PrintValsProdATP(m, Objective=None, ATPdemand=1.0, AA=None, o2=None, no3=None, blocks=[]):
    # "GLN_AA_mm_tx", "NO3_mm_tx" are not present in MMc media but are present in synovial fluid

    """pre: m = the model
            Objective = reaction which flux must be minimised as LP objective
            ATPdemand = if ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            AA = flux constraint through AA tx
            o2 = flux constraint through O2 tx
            no3 = flux constraint through NO3 tx
            blocks = reactions to be blocked
       post: returns a dictionary with solution for production of ATP"""

    sol, ObjVal, ByProds, CompUsed, TotalSolReac, SolReac = ProdATP(m, Objective, ATPdemand, AA, o2, no3, blocks)

    print "Total no of reactions in solution: ", len(TotalSolReac)
    print "Total no of reactions in solution without transporters: ", len(SolReac)
    print "Objective value of the solution in mmol/gDW/h: ", round(ObjVal,3)
    if "GLC_mm_tx" in sol.keys():
        print "This is the ammount of GLC utilised (mmol): ", round(sol["GLC_mm_tx"],3), "\n",
    else:
        print "No GLC is consumed", "\n",
    if "NH4_mm_tx" in sol.keys():
        print "This is the ammount of NH4 utilised/excreted (mmol): ", round(sol["NH4_mm_tx"],3), "\n",
    else:
        print "No NH4 is uptaken/excreted", "\n",
    print "These are the ammounts of by-products exported (mmol): ", "\n",
    for bp in sorted(ByProds.keys()):
        print bp, ":  ", round(ByProds[bp],3)
    print "These are the amounts of compounds uptaken/excreted (mmol): ", "\n",
    for comp in sorted(CompUsed.keys()):
        print comp, ":  ", round(CompUsed[comp],3)

def TotalCconsumed(res, CatomsAll=DataForAnalysis.Catoms_synovial):
    """pre: res = LP sol (e.g. ProdATP(m))
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
    """pre: res = LP sol (e.g. ProdATP(m))
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
    """pre: res = LP sol (e.g. ProdATP(m))
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
    """pre: res = LP sol (e.g. ProdATP(m))
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

def ATPglucose(res):
    """pre: res = LP sol (e.g. ProdATP(m))
       post: returns the ATP/GLC ratio for res"""

    if "GLC_mm_tx" in res.keys():
        ATPglc = res["ATPase"] / res["GLC_mm_tx"]  # ATP/GLC ratio
    else:
        ATPglc = 0.0

    print "This is the ATP/GLC ratio: ", round(ATPglc,3), "\n"

    return ATPglc

def ATPcarbon(res, CatomsAll=DataForAnalysis.Catoms_synovial):
    """pre: res = LP sol (e.g. ProdATP(m))
            CatomsAll = dic of n of C atoms in media components
       post: returns the ATP/C ratio for res"""

    TotalC = TotalCconsumed(res, CatomsAll)[1]
    #if "GLC_mm_tx" in res.keys():
    ATPperC = res["ATPase"] / TotalC  # ATP/C ratio
        
    print "This is the ATP/C ratio: ", round(ATPperC,3), "\n"

    return ATPperC

def OxidPOratio(res):     # Oxidative phosphorilation
    """pre: res = LP sol (e.g. ProdATP(m))
       post: returns the PO ratio for res"""

    if "O2_mm_tx" in res.keys() or "NO3_mm_tx" in res.keys():
        if "ATPSynth" in res.keys():
            if "O2_mm_tx" in res.keys() and "NO3_mm_tx" in res.keys(): 
                PO = res["ATPSynth"] / (res["O2_mm_tx"]*2 + res["NO3_mm_tx"]*3) 
            elif "O2_mm_tx" in res.keys() and "NO3_mm_tx" not in res.keys():
                PO = res["ATPSynth"] / (res["O2_mm_tx"]*2) 
            elif "NO3_mm_tx" in res.keys() and "O2_mm_tx" not in res.keys():
                PO = res["ATPSynth"] / (res["NO3_mm_tx"]*3)
            else:
                PO = 0.0
    else:
        PO = 0.0

    print "This is the P/O ratio for oxidative phosphorilation: ", round(PO,3), "\n"

    return PO

def OxidPlusSubsPOratio(res): # Oxidative plus substrate level phosphorilation
    """pre: res = LP sol (e.g. ProdATPm))
       post: returns the FullPO ratio for res"""

    if "O2_mm_tx" in res.keys() and "NO3_mm_tx" in res.keys():
        FullPO = res["ATPase"] / (res["O2_mm_tx"]*2 + res["NO3_mm_tx"]*3)
    elif "O2_mm_tx" in res.keys() and "NO3_mm_tx" not in res.keys():
        FullPO = res["ATPase"] / (res["O2_mm_tx"]*2) 
    elif "NO3_mm_tx" in res.keys() and "O2_mm_tx" not in res.keys():
        FullPO = res["ATPase"] / (res["NO3_mm_tx"]*3) 
    else:
        FullPO = 0.0

    print "This is the P/O ratio (total): ", round(FullPO,3),

    return FullPO

## Metabolic responses for ATP production under a range of conditions
    ## and minimising total flux or import of C source (Objective)
        ## Thesis chapter 3

def ATPprodResponse(m, Objective=None, ATPdemand=1.0, AA=None, o2=None, no3=None,
    blocks=[], CatomsAll=DataForAnalysis.Catoms_synovial):
    """pre: m = the model
            Objective = reaction which flux must be minimised as LP objective
            ATPdemand = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW = 45.0
            AA = flux constraint through AA tx
            o2 = flux constraint through O2 tx
            no3 = flux constraint through NO3 tx
            blocks = reactions to be blocked
            CatomsAll = dic of n of C atoms in media components
            NatomsAll = dic of n of N atoms in media components
       post: returns a dictionary with solution for production of ATP"""

    res = ProdATP(m, Objective, ATPdemand, AA, o2, no3, blocks)[0]

    PrintVals = PrintValsProdATP(m, Objective, ATPdemand, AA, o2, no3, blocks)

    ATPperGLCratio = ATPglucose(res)

    ATPperCratio = ATPcarbon(res, CatomsAll)

    POratioOxidPhosp = OxidPOratio(res)

    POratioOxidPlusSubstPhosp = OxidPlusSubsPOratio(res)

    return res
    

def ProdATPcharact(m, Objective=None, ATPdemand=45.0, AA=None, o2=None, no3=None, blocks=[],
    CatomsAll=DataForAnalysis.Catoms_synovial, NatomsAll=DataForAnalysis.Natoms_synovial):

    """pre: m = the model
            Objective = reaction which flux must be minimised as LP objective
            ATPdemand = if ATP maintenance demand (NGAM plus GAM),growth rate 1gDW/h = 45 mmol/gDW = 45.0
            AA = flux constraint through AA tx
            o2 = flux constraint through O2 tx
            no3 = flux constraint through NO3 tx
            blocks = reactions to be blocked
            CatomsAll = dic of n of C atoms in media components
            NatomsAll = dic of n of N atoms in media components
       post: returns a dictionary with solution for production of ATP"""

    res = ProdATP(m, Objective, ATPdemand, AA, o2, no3, blocks)[0] 

    CompContribCuptake = ContribToTotalCuptake(res, CatomsAll)

    CompContribNuptake = ContribToTotalNuptake(res, NatomsAll)

    ATPperGLCratio = ATPglucose(res)

    ATPperCratio = ATPcarbon(res, CatomsAll)

    POratioOxidPhosp = OxidPOratio(res)

    POratioOxidPlusSubstPhosp = OxidPlusSubsPOratio(res)

    return res
