######
########
##########   This module uses LPA for studying AA utilisation as C sources for BMP
##########
########    NOTE: remember, "GLN_AA_mm_tx" is blocked cause is not present in MMc media.
######


import sys
# sys.path.append('../Analysis/AAutil/')
from Util import Set
import BuildLP
from Bioinf import Biomass
import SepiBiomass
import math
import AAessential
import DataForAnalysis

######################### VERY IMPORTANT: check that we have the same CTxs list in module AAessential.py !!!!!
##
### C transporters potentially available to our model
####CTxs = [
####    'THR_AA_mm_tx',
####    'LEU_AA_mm_tx',
####    'PHE_AA_mm_tx',
####    'TYR_AA_mm_tx',
####    'HIS_AA_mm_tx',
####    'GLN_AA_mm_tx', ## reintated for completeness - turn it off if not in media
####    'GLT_AA_mm_tx',
####    'SER_AA_mm_tx',
####    'MET_AA_mm_tx',
####    'ASN_AA_mm_tx',
####    'PRO_AA_mm_tx',
####    'ALA_AA_mm_tx',
####    'ILE_AA_mm_tx',
####    'ASP_AA_mm_tx',
####    'GLY_AA_mm_tx',
####    'CYS_AA_mm_tx',
####    'LYS_AA_mm_tx',
####    'TRP_AA_mm_tx',
####    'VAL_AA_mm_tx',
####    'ARG_AA_mm_tx',
####    'GLC_mm_tx'
####]
##
#### If medium is synovial fluid use this list =
##
##CTxs = [
##    'THR_AA_mm_tx',
##    'LEU_AA_mm_tx',
##    'PHE_AA_mm_tx',
##    'TYR_AA_mm_tx',
##    'HIS_AA_mm_tx',
##    'GLN_AA_mm_tx', ## reintated for completeness - turn it off if not in media
##    'GLT_AA_mm_tx',
##    'SER_AA_mm_tx',
##    'MET_AA_mm_tx',
##    'ASN_AA_mm_tx',
##    'PRO_AA_mm_tx',
##    'ALA_AA_mm_tx',
##    'ILE_AA_mm_tx',
##    'ASP_AA_mm_tx',
##    'GLY_AA_mm_tx',
##    'CYS_AA_mm_tx',
##    'LYS_AA_mm_tx',
##    'TRP_AA_mm_tx',
##    'VAL_AA_mm_tx',
##    'ARG_AA_mm_tx',
##    'GLC_mm_tx',
##    'CITRULLINE_mm_tx', # "L-CITRULLINE" # C 6
##    'ORNIT_mm_tx', # "L-ORNITHINE" # C 5
##    '4-AMINO-BUTYRATE_mm_tx', # "4-AMINO-BUTYRATE" # C 4
##    'UREA_mm_tx' # "UREA" # C 1
##]

## Extra n sources in synovial fluid vs MM medium =
## ExtraNTxs = [
##    'NO3_mm_tx', # "NITRATE" # C 0
##    'CITRULLINE_mm_tx', # "L-CITRULLINE" # C 6
##    'ORNIT_mm_tx', # "L-ORNITHINE" # C 5
##    '4-AMINO-BUTYRATE_mm_tx', # "4-AMINO-BUTYRATE" # C 4
##    'UREA_mm_tx' # "UREA" # C 1
##]

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
            m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
       post: lp object constrained to biomass fluxes defined by BM composition
       and ATP maintenace demand"""

    lp = BuildLP.BuildLP(m) # lp = LP.lp(m) ; lp = m.GetLP() # lp.SetObjective(m.sm.cnames)
    fd = GetFluxDic(m, PropBF, AMaint)
    lp.SetFixedFlux(fd)

    return lp


def UtilSingleComps_C(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"],
    CTxs=DataForAnalysis.CTxs_MM, rm=[]):
## This fx checks for potential of each individual AA as sole C sources
    # for BMP upon standar conditions in MMc medium(19 AA and NH4 available)
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions which flux we want to block
            CTxs = list of tx of C sources available in media
            rm = tx reacs for C sources not available to remove from CTxs
        post: returns a dictionary with AAs as keys and values are tuples for solutions
        for BMP on the AA as single N source with [0] = LPA sol, [1] = the objective value;
        Also a list AAs as valid single C sources"""
    
    SingleAA = {}
    ValidSources = []
    NonValidSources = []

    lp = GetLP(m, PropBF, AMaint)
    
    for r in blocks:
        lp.SetFixedFlux({r:0.0})

    for comp in rm:
        CTxs.remove(comp)
    
    block = CTxs[:]
    BlockDic = dict(zip(block, [0.0]*len(block)))
    lp.SetFixedFlux(BlockDic)

    for c in Set.Complement(block, blocks):
        lp.ClearFluxConstraint(c)            
        lp.Solve(False)
        if lp.IsStatusOptimal():
            sol = lp.GetPrimSol()
            ObjVal = lp.GetObjVal()
            SingleAA[c] = sol, ObjVal           # return lp sol, ObjVal (tuple)
            lp.SetFixedFlux({c:0.0})
        else:
            SingleAA[c] = "Nosol", 0.0, 0.0      
            lp.SetFixedFlux({c:0.0})        

    for c in SingleAA.keys():
        if SingleAA[c][0] != "Nosol":
            ValidSources.append(c)

    for c in SingleAA.keys():
        if SingleAA[c][0] == "Nosol":
            NonValidSources.append(c)

    return SingleAA, ValidSources, NonValidSources

def PrintUtilSingleAAs_C(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"],
    CTxs=DataForAnalysis.CTxs_MM, rm=[]):
## This fx checks for potential of each individual AA as sole C sources for BMP
    ## and prints lists of valid AAs and non-valid AAs
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium

    """pre: UtilSingleAA(m...)
            m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions which flux we want to block
            CTxs = list of tx of C sources available in media
            rm = tx reacs for C sources not available to remove from CTxs
        post: prints values for sols for BMP in single AAs as C sources:
        objective value, glucose consumption, NH4 excreted"""

    ValsUtilSingleAAs_C = UtilSingleComps_C(m, PropBF, AMaint, blocks, CTxs, rm)

    print "These are the single medium compounds that can be utilised as sole C sources: ", sorted(ValsUtilSingleAAs_C[1]), "\n"
    print "These are the single medium compounds that cannot be utilised as sole C sources: ", sorted(ValsUtilSingleAAs_C[2]), "\n"

def PrintValsUtilSingleAAs_C(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"],
    CTxs=DataForAnalysis.CTxs_MM, rm=[]):
## This fx checks for potential of each individual AA as sole C sources for BMP
    ## and print relevant values for the solutions found
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium
    """pre: UtilSingleAA(m...)
            m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions which flux we want to block
            CTxs = list of tx of C sources available in media
            rm = tx reacs for C sources not available to remove from CTxs
        post: prints values for sols for BMP in single AAs as C sources:
        amounts consumed and objective values"""

    ValsUtilSingleAAs_C = UtilSingleComps_C(m, PropBF, AMaint, blocks, CTxs, rm)

    print "\n", "These are the amounts of medium compounds consumed (mmol) for BMP when utilising a single medium compounds as C sources: ", "\n"
    for n in sorted(ValsUtilSingleAAs_C[1]):
        print n, "","","","", round(ValsUtilSingleAAs_C[0][n][0][n],3)
    print "\n", "These are the objective values (mmol/gDW/h) for BMP utilising single medium compounds as C sources: ", "\n"
    for n in sorted(ValsUtilSingleAAs_C[1]):
        print n, "","","","", round(ValsUtilSingleAAs_C[0][n][1],3)


def ContribCuptakeAAs(m, PropBF=0.0, AMaint=45.0, AA=None,
    CTxs=DataForAnalysis.CTxs_MM, rm=[], blocks=["GLN_AA_mm_tx"],
    Catoms=DataForAnalysis.Catoms_MM_AA):
## Ratio uptake vs demand and contribution to total C uptake
    # for 1 unit BMP upon standar conditions in MMc medium(19 AA and NH4 available)
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium

    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
            CTxs = list of tx of C sources available in media
            rm = tx reacs for C sources not available to remove from CTxs
            Catoms = n C atoms in AAs in medium
        post: returns dict mapping AA txs to their ratio of uptake/demand
        and proportional contribution to total C uptake"""

    lp = GetLP(m, PropBF, AMaint)
    
    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})

    for r in blocks:
        lp.SetFixedFlux({r:0.0})

    for comp in rm: # rm = list of N sources not available
        CTxs.remove(comp)

    lp.Solve(False)
    sol = lp.GetPrimSol()

    AAReqs = AARequirement(m, PropBF, AMaint, CTxs) # dic of fluxes required in the AA biomass output

    sol, CuptakeGLC, CuptakeAAs, TotalCuptaken = TotalCuptakeAAs(m, PropBF, AMaint, AA, blocks, Catoms)
    
    TotalC = TotalCuptaken


    CperAA = {}     # C uptaken as AAs
    for aa in Catoms.keys():
        if aa in sol.keys():
            CperAA[aa] = Catoms[aa]*sol[aa]
        else:
            CperAA[aa] = 0.0

    rv = {}
    for aa in AAReqs.keys():           # for those aa required in the AA biomass output  
        if aa not in sol.keys(): # if AA not uptaken
            propUD = propT = 0.0
        else:
            if AAReqs[aa] == 0.0:
                propUD = 0.0  # ratio uptake/demand = flux in sol (uptaken) / flux for bm output (demand for BMP)
                propT = CperAA[aa]/TotalC   # prop contrib to total C uptake = flux in sol*number C atoms / total C demand
                                            # and / total C demand
            else:
                propUD = sol[aa]/AAReqs[aa]  # ratio uptake/demand = flux in sol (uptaken) / flux for bm output (demand for BMP)
                propT = CperAA[aa]/TotalC   # prop contrib to total C uptake = flux in sol*number C atoms / total C demand
                                            # and / total C demand 
        rv[aa] = propUD, propT

    if "GLC_mm_tx" in sol.keys():
        GLCcontrib = round(sol["GLC_mm_tx"]*6/TotalC *100,3)
    else: GLCcontrib = 0.0

    propTdic = {}
    for aa in AAReqs:
        if aa not in sol.keys(): # if AA not uptaken
            propTdic[aa] = 0.0
        else:
            propTdic[aa] = CperAA[aa]/TotalC

    Totalpercent = round(((sum(propTdic.values())*100) + GLCcontrib),3)

    print "\n", "These are the AA uptake/demand ratios for BMP (1 gram): ", "\n"
    for aa in sorted(AAReqs.keys()):
        print aa, "","","","", round(rv[aa][0],3)
    print "\n", "This is the total C demand under these conditions (mmol): ", "\n", round(TotalC,3)    
    print "\n", "These are the AA contributions to total C uptaken from AAs (percentages): ", "\n"
    for aa in sorted(AAReqs.keys()):
        print aa, "","","","", round(rv[aa][1]*100,3)
    print "\n", "This is the GLC contribution to total C uptake (percentage): ",GLCcontrib, "\n", 
    print "\n", "Sum: ", Totalpercent
    
    return rv

def ContribCuptakeALL(m, PropBF=0.0, AMaint=45.0, AA=None,
    CTxs=DataForAnalysis.CTxs_MM, rm=[], blocks=["GLN_AA_mm_tx"],
    Catoms=DataForAnalysis.Catoms_MM, FreeBM=None):
## Ratio uptake vs demand and contribution to total C uptake
    # for 1 unit BMP upon standar conditions in MMc medium(19 AA and NH4 available)
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium

    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
            CTxs = list of tx of C sources available in media
            rm = tx reacs for C sources not available to remove from CTxs
            Catoms = n C atoms in AAs in medium
            FreeBM = if != None, upper flux bound of BM_tx (biomass txs) is unbound
        post: returns dict mapping AA txs to their ratio of uptake/demand
        and proportional contribution to total C uptake"""

    if FreeBM != None:
        lp = GetLPfreeBM(m, PropBF, AMaint)
    else:
        lp = GetLP(m, PropBF, AMaint)
    
    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})

    for r in blocks:
        lp.SetFixedFlux({r:0.0})

    for comp in rm: # rm = list of N sources not available
        CTxs.remove(comp)

    lp.Solve(False)
    sol = lp.GetPrimSol()

    AAReqs = AARequirement(m, PropBF, AMaint, CTxs) # dic of fluxes required in the AA biomass output

    sol, CuptakeGLC, TotalCuptaken = TotalCuptakeALL(m, PropBF, AMaint, AA, blocks, Catoms)

    TotalC = TotalCuptaken

    CperComp = {}     # C uptaken
    for comp in Catoms.keys():
        if comp in sol.keys():
            CperComp[comp] = Catoms[comp]*sol[comp]
        else:
            CperComp[comp] = 0.0

    rv = {}
    for comp in Catoms.keys():   # for those comp in medium containing C atoms 
        if comp not in sol.keys():  # if AA not uptaken
            propUD = propT = 0.0
        else:
            if comp not in AAReqs.keys():
                propUD = 0.0                    # ratio uptake/demand = flux in sol (uptaken) / flux for bm output (demand for BMP)
                propT = CperComp[comp]/TotalC   # prop contrib to total C uptake = flux in sol*number C atoms / total C uptaken
            else:
                if AAReqs[comp] == 0.0:         # if comp is not part of the BM (zero demand)
                    propUD = 0.0              
                    propT = CperComp[comp]/TotalC   
                else:
                    propUD = sol[comp]/AAReqs[comp]     # ratio uptake/demand = flux in sol (uptaken) / flux for bm output (demand for BMP)
                    propT = CperComp[comp]/TotalC    
        rv[comp] = propUD, propT

    if "GLC_mm_tx" in sol.keys():
        GLCcontrib = round(sol["GLC_mm_tx"]*6/TotalC *100,3)
    else: GLCcontrib = 0.0

    propTdic = {}
    for comp in Catoms:
        if comp not in sol.keys():  # if comp in medium is not uptaken 
            propTdic[comp] = 0.0
        else:
            propTdic[comp] = CperComp[comp]/TotalC

    Totalpercent = round((sum(propTdic.values())*100),3)

    print "\n", "These are the AA uptake/demand ratios for BMP (1 gram): ", "\n"
    for aa in sorted(AAReqs.keys()):
        print aa, "","","","", round(rv[aa][0],3)
    print "\n", "This is the total C demand under these conditions (mmol): ", "\n", round(TotalC,3)  
    print "\n", "These are the contributions of medium components to total C uptake (percentage): ", "\n"
    for comp in sorted(Catoms.keys()):
        print comp, "","","","", round(rv[comp][1]*100,3)
    print "\n", "Sum: ", Totalpercent
    
    return rv

####################################################################  FREE BM (no upper flux bounds)
## Analyse system for production of 1 unit of cell biomass (BMP)
    # while meeting ATP cell maintenance requirements
        # and allowing BM prods to be exported at a higher proportion
            # than needed to produce BM!

def GetLPfreeBM(m, PropBF=0.0, AMaint=45.0):

    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
       post: lp object constrained to biomass fluxes defined by BM composition
       and ATP maintenace demand"""

    orig_rp = dict(m.sm.RevProps)
    for r in filter(lambda s: "bm_tx" in s, m.sm.cnames):
        m.sm.RevProps[r] = m.smx.RevProps[r] = "<>"

    lp = BuildLP.BuildLP(m) # lp = LP.lp(m) ; lp = m.GetLP() # lp.SetObjective(m.sm.cnames)

    fd = GetFluxDic(m, PropBF, AMaint)
    
    for r in fd.keys():
        if 'bm_tx' in r:    # otherwise allows ATPase to go to zero flux! (since also in fd)
            lp.SetFluxBounds({r:(None,fd[r])}) # Bug = takes val of '0' if r not declared REV
#	lp.SetFluxBounds({r:(fd[r],None)})  # Bug = takes positive and negative values! if r not declared REV
    m.sm.RevProps = m.smx.RevProps = orig_rp
    
    lp.SetFluxBounds({"ATPase": (AMaint, None)}) # If we dont want ATPase upper flux to be fixed
                    
    return lp

####################################################################### ENDS FREE BM 

def AARequirement(m, PropBF=0.0, AMaint=45.0, CTxs=DataForAnalysis.CTxs_MM):
## This function matches the fluxes through AA importers available to the flux   
    ## requiered through the AA exporters for BMP of 1 gram of planktonic biomass
    """pre: mm = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            CTxs = list of tx of C sources available in media
            rm = tx reacs for C sources not available to remove from CTxs
        post: returns dict mapping AA importers in media to the flux
        through the AA exporters for production of BM (1 gram)"""
    rv = {}
    fd = GetFluxDic(m, PropBF, AMaint)

    CInputs = CTxs[:]   # list of all possible C sources (media txs)

    COutputs = map(lambda a: a.replace("_AA_mm_","_AA_bm_"), CTxs) # replace in the list the _AA_mm_ suffix by _bm_
    Out2In = dict(zip(COutputs, CInputs))                          # dic mapping txs:  AA_bm_tx to _AA_mm_tx
    UsedOutputs = Set.Intersect(COutputs, fd)                   # list with those AA_bm_tx generated from CTxs
                                                                # that are also in fd for BMP
    
    for out in UsedOutputs:             # for AA_bm_tx in fd for which AA_mm_tx in media (comps in medium with bm_tx in model)
        rv[Out2In[out]] = -fd[out]      # rv is a dic with their AA_mm_tx as keys
                                        # and flux is same flux as for AA_bm_tx for BMP (but positive)
    return rv

def TotalCuptakeAAs(m, PropBF=0.0, AMaint=45.0, AA=None, blocks=["GLN_AA_mm_tx"],
        Catoms=DataForAnalysis.Catoms_MM_AA):
## This function calculates total C uptaken for BMP (1 gram) from AAs in medium
   ## and ATP maintenance upon standar conditions, MMc medium
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            Catoms = n C atoms C in AAs in medium 
        post: returns solution for BMP in standar MM (GLC and AAs present);
        Plus values for C uptaken"""

    lp = GetLP(m, PropBF, AMaint)
    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})
    for r in blocks:
        lp.SetFixedFlux({r:0.0})
        
    lp.Solve(False)
    sol = lp.GetPrimSol()

    CperAA = {}     # C uptaken as AAs
    for aa in Catoms.keys():
        if aa in sol.keys():
            CperAA[aa] = Catoms[aa]*sol[aa]
        else:
            CperAA[aa] = 0.0

    print "\n", "These are the mmol of C uptaken as AAs for BMP: ", "\n"
    for aa in sorted(CperAA.keys()):
        print aa, "","","","", round(CperAA[aa],3)
        
    if "GLC_mm_tx" in sol.keys():
        CuptakeGLC = round(((sol["GLC_mm_tx"])*6),3)   # C in equals C out (as BM and by-products)
    else: CuptakeGLC = 0.0
    CuptakeAAs = round(sum(CperAA.values()),3)
    TotalCuptaken = round(CuptakeGLC + CuptakeAAs,3)
    print "\n", "This is the total C uptaken as GLC (mmol): ", CuptakeGLC
    print "\n", "This is the total C uptaken as AAs (mmol): ", CuptakeAAs
    print "\n", "This is the total C uptaken (mmol): ", TotalCuptaken

    return sol, CuptakeGLC, CuptakeAAs, TotalCuptaken

def TotalCuptakeALL(m, PropBF=0.0, AMaint=45.0, AA=None, blocks=["GLN_AA_mm_tx"],
        Catoms=DataForAnalysis.Catoms_MM, FreeBM=None):
## This function calculates total C uptaken for BMP (1 gram) from all C copntaining comps in medium
   ## and ATP maintenance upon standar conditions, MMc medium
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium

    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            Catoms = n C atoms in ncompounds in medium
            FreeBM = if != None, upper flux bound of BM_tx (biomass txs) is unboun
        post: returns solution for BMP in standar MM (GLC and AAs present);
        Plus values for C uptaken"""

    if FreeBM != None:
        lp = GetLPfreeBM(m, PropBF, AMaint)
    else:
        lp = GetLP(m, PropBF, AMaint)

    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})
    for r in blocks:
        lp.SetFixedFlux({r:0.0})
        
    lp.Solve(False)
    sol = lp.GetPrimSol()
    
    CperComp = {}     # C uptaken as comps from medium
    for comp in Catoms.keys():
        if comp in sol.keys():
            CperComp[comp] = Catoms[comp]*sol[comp]
        else:
            CperComp[comp] = 0.0

    print "\n", "These are the mmol of C uptaken as medium compounds for BMP: ", "\n"
    for comp in sorted(CperComp.keys()):
        print comp, "","","","", round(CperComp[comp],3)
        
    if "GLC_mm_tx" in sol.keys():
        CuptakeGLC = round(((sol["GLC_mm_tx"])*6),3)   # C in equals C out (as BM and by-products)
    else: CuptakeGLC = 0.0

    TotalCuptaken = round(sum(CperComp.values()),3)
    print "\n", "This is the total C uptaken as GLC (mmol): ", CuptakeGLC
    print "\n", "This is the total C uptaken (mmol): ", TotalCuptaken

    return sol, CuptakeGLC, TotalCuptaken

def TotalCexcret(m, PropBF=0.0, AMaint=45.0, AA=None, blocks=["GLN_AA_mm_tx"],
        CatomsBM=DataForAnalysis.Catoms_BM_AA, CatomsBPs=DataForAnalysis.CatomByProds):
## This function calculates total C excreted as AAs in BM and by-products
    ## during BMP (1 gram) in MM media (GLC plus 19 AAs)
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions which flux we want to block
            CatomsBM = n C atoms in amino acids in the biomass
            CatomsBPs = n C atoms in compounds potentially excreted as BPs
        post: returns solution for BMP in standar MM (GLC and AAs present);
        Plus values for C excreted as AAs and by-products"""

    lp = GetLP(m, PropBF, AMaint)
    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})
    for r in blocks:
        lp.SetFixedFlux({r:0.0})
    lp.Solve(False)
    sol = lp.GetPrimSol()
    
    CperAA = {}     # C excreted as AAs in BM
    for aa in CatomsBM.keys():
        if aa in sol.keys():
            CperAA[aa] = CatomsBM[aa]*sol[aa]
        else :
            CperAA[aa] = 0.0

    print "\n", "These are the mmol of C excreted as AAs in 1 gram of BM: ", "\n"
    for aa in sorted(CperAA.keys()):
        print aa, "","","","", round(CperAA[aa],3)
    
    CperByprod = {}     # C excreted as by-products when producing 1 gram of BM
    for bp in CatomsBPs.keys():
        if bp in sol.keys():
            CperByprod[bp] = CatomsBPs[bp]*sol[bp]

    print "\n", "These are the mmol of C excreted as by-products during BMP: ", "\n"
    for bp in sorted(CperByprod.keys()):
        print bp, "","","","", round(CperByprod[bp],3)

    CexcretedAAs = round(sum(CperAA.values()),3)
    CexcretedBPs = round(sum(CperByprod.values()),3)
    TotalCexcreted = round(CexcretedAAs + CexcretedBPs,3)
    print "\n", "This is the total C excreted as AAs (mmol): ", CexcretedAAs
    print "\n", "This is the total C exreted as by-products (mmol): ", CexcretedBPs
    print "\n", "This is the total C excreted as AAs and by-products (mmol): ", TotalCexcreted
    return sol, CexcretedAAs, CexcretedBPs, TotalCexcreted
