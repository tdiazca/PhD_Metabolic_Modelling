######
########
##########   This module uses LPA for studying AA utilisation as N sources for BMP
##########
########    NOTE: remember, "GLN_AA_mm_tx" is blocked because are not present in MMc media
######          
####                available in synovial fluid so do not block in this case (also NO3)
##          NOTE: select adequate NTxs list (uncomment list of N containing compounds depending on medium)


import sys
# sys.path.append('../Analysis/AAutil/')
from Util import Set
import BuildLP
from Bioinf import Biomass
import SepiBiomass
import math
import AAessential
import DataForAnalysis

###### N and C transporters potentially available to our model

## select right list from DataForAnalysis file depending on media:

    # NTxs_MinMed ; NTxs_MM ; NTxs_synovial

    # CTxs_MinMed ; CTxs_MM ; CTxs_synovial

## # AA defined with AAutilN.UtilSingleAA(m) as non-viable single N sources.

    # NonViabAA

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

def UtilAA(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"],
           NTxs=DataForAnalysis.NTxs_MM, rm=[]):
## This fx calculates the overall uptake or excretion of AA
    # for 1 unit BMP upon standar conditions in MMc medium(19 AA and NH4 available)
            # "GLN_AA_mm_tx" is blocked cause is not present in MMc medium
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
        post: returns a dictionary with AAs in media as keys and values are
        total amount of AA uptaken/excreted during BMP (mmol/gram cDW/hour)"""

    StdBMP = AAessential.StandarBMPinMM(m, PropBF, AMaint, blocks)
    fd = GetFluxDic(m, PropBF, AMaint)  # fd for BMP (1 unit) and AMaint of m

    AAavail = NTxs[:]

    for comp in rm:
        AAavail.remove(comp)

    AAuptake = {}
    for aa in AAavail:
        if aa not in StdBMP[0].keys():
            AAuptake[aa] = 0.0
        else:
            AAuptake[aa] = StdBMP[0][aa]

    AAuptake2 = {}
    for aa in AAuptake.keys():
        if '_AA_mm_tx' in aa:
            AAuptake2[aa.replace('_AA_mm_tx','')] = AAuptake[aa]
        else:
            if '_mm_tx' in aa:
                AAuptake2[aa.replace('_mm_tx','')] = AAuptake[aa]
                
    AAexcret = {}
    for aa in filter(lambda s: "_AA_bm_tx" in s, m.smx.cnames):
        if aa not in StdBMP[0].keys():
            AAexcret[aa] = 0.0
        else:
            AAexcret[aa] = StdBMP[0][aa]

    AAexcret2 = {}
    for aa in AAexcret.keys():
        AAexcret2[aa.replace('_AA_bm_tx','')] = AAexcret[aa]

    AAtotalFlux = {}
    
    NewKeys = []
    for aa in AAavail:
        if '_AA_mm_tx' in aa:
            NewKeys.append(aa.replace('_AA_mm_tx',''))
        else:
            if '_mm_tx' in aa:
                NewKeys.append(aa.replace('_mm_tx',''))

    for aa in NewKeys:
        if aa in AAexcret2.keys():
            AAtotalFlux[aa] = AAuptake2[aa] + AAexcret2[aa]
        else:
            AAtotalFlux[aa] = AAuptake2[aa]

    print "Total (overall) amounts of AA uptaken/excreted during BMP (mmol/gram cDW/hour): ", "\n"
    for aa in sorted(AAtotalFlux):
        print aa, ":  ", round(AAtotalFlux[aa],3) 

    return AAtotalFlux, StdBMP

##########################################################  FREE BM (no upper flux bounds)
##############
####

def UtilAAfreeBM(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"],
                 NTxs=DataForAnalysis.NTxs_MM, rm=[]):
## This fx calculates the overall uptake or excretion of AA
    # for 1 unit BMP upon standar conditions in MMc medium(19 AA and NH4 available)
            # "GLN_AA_mm_tx" is blocked cause is not present in MMc medium

    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
        post: returns a dictionary with AAs in media as keys and values are
        total amount of AA uptaken/excreted during BMP (mmol/gram cDW/hour)"""

    StdBMP = AAessential.StandarBMPinMMfreeBM(m, PropBF, AMaint, blocks)
    fd = GetFluxDic(m, PropBF, AMaint)  # fd for BMP (1 unit) and AMaint of m

    AAavail = NTxs[:]
##    rm = ["NH4_mm_tx","NO3_mm_tx","CITRULLINE_mm_tx","ORNIT_mm_tx",
##          "4-AMINO-BUTYRATE_mm_tx","UREA_mm_tx"]
    for comp in rm:
        AAavail.remove(comp)

    AAuptake = {}
    for aa in AAavail:
        if aa not in StdBMP[0].keys():
            AAuptake[aa] = 0.0
        else:
            AAuptake[aa] = StdBMP[0][aa]

    AAuptake2 = {}
    for aa in AAuptake.keys():
        if '_AA_mm_tx' in aa:
            AAuptake2[aa.replace('_AA_mm_tx','')] = AAuptake[aa]
        else:
            if '_mm_tx' in aa:
                AAuptake2[aa.replace('_mm_tx','')] = AAuptake[aa]

    AAexcret = {}
    for aa in filter(lambda s: "_AA_bm_tx" in s, m.smx.cnames):
        if aa not in StdBMP[0].keys():
            AAexcret[aa] = 0.0
        else:
            AAexcret[aa] = StdBMP[0][aa]

    AAexcret2 = {}
    for aa in AAexcret.keys():
        AAexcret2[aa.replace('_AA_bm_tx','')] = AAexcret[aa]

    AAtotalFlux = {}
    
    NewKeys = []
    for aa in AAavail:
        if '_AA_mm_tx' in aa:
            NewKeys.append(aa.replace('_AA_mm_tx',''))
        else:
            if '_mm_tx' in aa:
                NewKeys.append(aa.replace('_mm_tx',''))

    for aa in NewKeys:
        if aa in AAexcret2.keys():
            AAtotalFlux[aa] = AAuptake2[aa] + AAexcret2[aa]
        else:
            AAtotalFlux[aa] = AAuptake2[aa]

    print "Total (overall) amounts of AA uptaken/excreted during BMP (mmol/gram cDW/hour): ", "\n"
    for aa in sorted(AAtotalFlux):
        print aa, ":  ", round(AAtotalFlux[aa],3) 

    return AAtotalFlux, StdBMP

####################################################################### ENDS FREE BM

def UtilSingleAA(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"],
                 NTxs=DataForAnalysis.NTxs_MM, rm=[]):
## This fx checks for potential of each individual AA as sole N sources for BMP
    # for BMP upon standar conditions in MMc medium(19 AA and NH4 available)
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium

    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions which flux we want to block
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
        post: returns a dictionary with AAs as keys and values are tuples for solutions
        for BMP on the AA as single N source with [0] = LPA sol, [1] = the objective value,
        [2] = the glucose consumption;[3] = the NH4 excreted;
        also a list AAs as valid single N sources also a list AAs as valid single N sources"""

    SingleAA = {}
    ValidSources = []
    NonValidSources = []

    lp = GetLP(m, PropBF, AMaint)
    #lp = GetLPfreeBM(m, PropBF, AMaint) # needed for Met util as sole N source 
    for r in blocks:
        lp.SetFixedFlux({r:0.0})

    for comp in rm:
        NTxs.remove(comp)

    block = NTxs[:]
    block.remove("NH4_mm_tx")                   # block is NTxs list with NH4 removed
    BlockDic = dict(zip(block, [0.0]*len(block))) # make dir for N tx with 0 as values
    lp.SetFixedFlux(BlockDic)
    lp.SetFluxBounds({"NH4_mm_tx":(None,0.0)})

    for n in Set.Complement(block, blocks):
        lp.ClearFluxConstraint(n)               # remove constrain (1 at a time)
        lp.Solve(False)
        if lp.IsStatusOptimal():
            sol = lp.GetPrimSol()
            ObjVal = lp.GetObjVal()
            GLCcon = sol["GLC_mm_tx"]
            SingleAA[n] = sol, ObjVal, GLCcon # return lp sol, ObjVal, sol["GLC_mm_tx"] (tuple)
            lp.SetFixedFlux({n:0.0})
        else:
            SingleAA[n] = "Nosol", 0.0, 0.0      # ObjVal = 0.0, and no solution possible
            lp.SetFixedFlux({n:0.0})            # reinstate flux constraint for that specific NTx

    for n in SingleAA.keys():
        if SingleAA[n][0] != "Nosol":
            ValidSources.append(n)

    for n in SingleAA.keys():
        if SingleAA[n][0] == "Nosol":
            NonValidSources.append(n)

    return SingleAA, ValidSources, NonValidSources

def PrintUtilSingleAAs(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"],
    NTxs=DataForAnalysis.NTxs_MM, rm=[]):
## This fx checks for potential of each individual AA as sole N sources for BMP
    ## and prints lists of valid AAs and non-valid AAs
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium

    """pre: UtilSingleAA(m...)
            m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions which flux we want to block
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
        post: prints values for sols for BMP in single AAs as N sources:
        objective value, glucose consumption, NH4 excreted"""

    ValsUtilSingleAAs = UtilSingleAA(m, PropBF, AMaint, blocks, NTxs, rm)

    print "These are the single AAs that can be utilised as sole N sources: ", sorted(ValsUtilSingleAAs[1]), "\n"
    print "These are the single AAs that cannot be utilised as sole N sources: ", sorted(ValsUtilSingleAAs[2]), "\n"

def PrintValsUtilSingleAAs(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"],
    NTxs=DataForAnalysis.NTxs_MM, rm=[]):
## This fx checks for potential of each individual AA as sole N sources for BMP
    ## and print relevant values for the solutions found
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium

    """pre: UtilSingleAA(m...)
            m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions which flux we want to block
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
        post: prints values for sols for BMP in single AAs as N sources:
        objective value, glucose consumption, NH4 excreted"""

    ValsUtilSingleAAs = UtilSingleAA(m, PropBF, AMaint, blocks, NTxs, rm)

    print "\n", "These are the amounts of AA consumed (mmol) for BMP when utilising a single AA as N source: ", "\n"
    for n in sorted(ValsUtilSingleAAs[1]):
        print n, "","","","", round(ValsUtilSingleAAs[0][n][0][n],3)
    print "\n", "These are the objective values (mmol/gDW/h) for BMP utilising a single AA as N source: ", "\n"
    for n in sorted(ValsUtilSingleAAs[1]):
        print n, "","","","", round(ValsUtilSingleAAs[0][n][1],3)
    print "\n", "These are the mmol of glucose consumed for BMP utilising a single AA as N source: ", "\n"
    for n in sorted(ValsUtilSingleAAs[1]):
        print n, "","","","", round(ValsUtilSingleAAs[0][n][2],3)
    print "\n", "These are the mmol of NH4 exported: ", "\n"
    for n in sorted(ValsUtilSingleAAs[1]):
        if "NH4_mm_tx" in ValsUtilSingleAAs[0][n][0].keys():
            print n, "","","","", round(ValsUtilSingleAAs[0][n][0]["NH4_mm_tx"],3)
    print "\n", "These are the number of reactions in the solutions: ", "\n"
    for n in sorted(ValsUtilSingleAAs[1]):
        print n, "","","","", len(ValsUtilSingleAAs[0][n][0])

def AACheckProds (m, blocks=["GLN_AA_mm_tx"], NTxs=DataForAnalysis.NTxs_MM, rm=[]):
## This function checks which AAs can be produced when utilising
        ## one single AA at a time as sole N source.
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium

    """pre: m = the model
            blocks = reactions which flux we want to block
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
        post: returns a dic with AA avail as keys and solutions for
        individual AA production as values"""

    AAsfromAA = {}

    lp = BuildLP.BuildLP(m) # lp = m.GetLP() # lp.SetObjective(m.sm.cnames)

    for r in blocks:
        lp.SetFixedFlux({r:0.0})

    for comp in rm:
        NTxs.remove(comp)

    block = NTxs[:]                                 # NTxs = all N sources in MM media
    block.remove("NH4_mm_tx")                       # block is NTxs list with NH4 removed
    BlockDic = dict(zip(block, [0.0]*len(block)))   # make dir for N tx with 0 as values
    lp.SetFixedFlux(BlockDic)
    lp.SetFluxBounds({"NH4_mm_tx":(None,0.0)})      # so NH4 can be exported but not imported as N source

    #lp.SetFixedFlux({"GLC_mm_tx":0.0}) # If we wanted to assess AA as sole C source

    for AA in Set.Complement(block, blocks):
        res = {}
        lp.ClearFluxConstraint(AA)
        for prod in filter(lambda s: "AA_bm_tx" in s, m.smx.cnames): #for AA in bm
            lp.SetFixedFlux({prod:-1})  # export one unit of prod
            lp.Solve(False)
            sol=lp.GetPrimSol()
            if lp.IsStatusOptimal():  # stores solution if it is feasible (as values)
                res[prod] = sol
            else:
                res[prod] = "Not produced"
            lp.ClearFluxConstraint(prod) # Clear constraint so m back to normal
                
        AAsfromAA[AA] = res
        lp.SetFixedFlux({AA:0.0}) # reinstate flux constraint for that specific AAtx

    for AA in AAsfromAA.keys():
        print "With ",AA, "as sole N source available",
        prodsNotProduced = []
        for prod in AAsfromAA[AA].keys():
            if AAsfromAA[AA][prod] == "Not produced":
                prodsNotProduced.append(prod)
        print "the following AAs cannot be synthesised: ", "\n", prodsNotProduced, "\n"

    return AAsfromAA

def EffectBMPinSingleAA(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"],
        NTxs=DataForAnalysis.NTxs_MM, rm=[]):
## This fx calculates the effect of BMP and ATP maintenace
    ##  with single individual AA/NH4 as sole N sources in ObjVal and GLC consumption
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium

    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            NTxs = list of tx of N sources available in media
            blocks = reactions to be blocked
            rm = tx reacs for N sources not available to remove from NTxs
        post: returns dicts with changes in objective value and GLC consumption
        expressed as percentages and has percentage of variation"""

    SingleAA, ValidSources, NonValidSources = UtilSingleAA(m, PropBF, AMaint, blocks, NTxs, rm)
    StdBMP = AAessential.StandarBMPinMM(m, PropBF, AMaint, blocks)
    BMPnh4 = AAessential.BMPonlyNH4(m, PropBF, AMaint, blocks, NTxs, rm)

    effObjValperc = {}

    effGLCperc = {}

    effObjValvar = {}

    effGLCvar = {}

    effObjValpercNH4 = {}

    effGLCpercNH4 = {}

    effObjValvarNH4 = {}

    effGLCvarNH4 = {}

    for n in ValidSources:
        effObjValperc[n] = (SingleAA[n][1] *100) / StdBMP[1]

    for n in ValidSources:
        effObjValvar[n] =  effObjValperc[n] - 100

    for n in ValidSources:
        effGLCperc[n] = (SingleAA[n][2] *100) / StdBMP[2]

    for n in ValidSources:
        effGLCvar[n] =  effGLCperc[n] - 100

    effObjValpercNH4["AllAArm"] = (BMPnh4[1] *100) / StdBMP[1]
    effGLCpercNH4["AllAArm"] = (BMPnh4[2] *100) / StdBMP[2]
    effObjValvarNH4["AllAArm"] = effObjValpercNH4["AllAArm"] - 100
    effGLCvarNH4["AllAArm"] = effGLCpercNH4["AllAArm"] - 100

    print  "This is the effect of BMP on single AAs on the objective value: ", "\n"
    for n in sorted(ValidSources):
        print n, "", "", "%", round(effObjValperc[n],3), "", "", "variation", "","", round(effObjValvar[n],3)
    print "\n", "This is the effect of BMP on single AAs on the GLC consumption: ", "\n"
    for n in sorted(ValidSources):
        print n, "", "","","", "%", round(effGLCperc[n],3), "", "","","", "variation", "","", round(effGLCvar[n],3)
    print "\n", "This is the effect of BMP on just NH4 (no AAs) on the objective value: "
    print "\n", "AllAArm", "", "", "%", round(effObjValpercNH4["AllAArm"],3), "", "", "variation", "","", round(effObjValvarNH4["AllAArm"],3)

    print "\n", "This is the effect of BMP on just NH4 (no AAs) upon the GLC consumption: "
    print "\n", "AllAArm", "", "", "%", round(effGLCpercNH4["AllAArm"],3), "", "", "variation", "","", round(effGLCvarNH4["AllAArm"],3)

    return effObjValperc, effObjValvar, effGLCperc, effGLCvar, effObjValpercNH4, effGLCpercNH4, effObjValvarNH4, effGLCvarNH4

###########

def ForceAAuptake(m, PropBF=0.0, AMaint=45.0, AA=["GLY_AA_mm_tx"],
        blocks=["GLN_AA_mm_tx"]):
## This fx investigates LP sol obtained upon forcing uptake
    ## of certain amino acids (AA); Here GLY since it is not uptaken in standar sol
        ## for BMP upon standar conditions in MMc medium (19 AA and NH4 available)
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            AA = amino acid in MM media forced to be uptaken
            blocks = reactions which flux we want to block
        post: returns a tuple with [0] = solution to FBA (dic), [1] = its objective value,
        [2] = its glucose consumption"""

    ForcedUptake = ()  # tuple
    fd = GetFluxDic(m, PropBF, AMaint)  # fd for BMP (1 unit) and AMaint of m
    lp = GetLP(m, PropBF, AMaint)

    for r in blocks:
        lp.SetFixedFlux({r:0.0})

    aaBM = map(lambda aa: aa.replace("_AA_mm_","_AA_bm_"), AA)

    BMflux = {}
    for aa in aaBM:
        BMflux[aa] = -fd[aa]

    MapTxs = {}
    for aa in AA:
        MapTxs[aa] = aa.replace("_AA_mm_","_AA_bm_")

    SetFlux = {}
    for aa in AA:
        SetFlux[aa] = BMflux[MapTxs[aa]]
    
    for aa in SetFlux.keys():
        lp.SetFixedFlux({aa:SetFlux[aa]})
        
    lp.Solve(False)
    sol = lp.GetPrimSol()
    ObjVal = lp.GetObjVal()
    GLCcomp = sol["GLC_mm_tx"]
    
    ForcedUptake = sol, ObjVal, GLCcomp # return ObjVal and lp sol

    print "This is the objective value (mmol/gDW/h) for BMP in MM forcing uptake of", aa, ": ", round(ForcedUptake[1],3), "\n"
    print "This is the number of reactions in the solution for BMP in MM forcing uptake of", aa, ": ", len(ForcedUptake[0]), "\n"
    print "This is the GLC consumed (mmol) for BMP in MM forcing uptake of", aa, ": ", round(ForcedUptake[2],3),"\n"
    print "This is the NH4 imported/excreted (mmol): ", round(ForcedUptake[0]["NH4_mm_tx"],3), "\n"
    print "These are by-products excreted (mmol): " ,"\n"
    for r in sorted(ForcedUptake[0]):
        if r in filter(lambda s: "_bp_tx" in s, m.smx.cnames):
            print r, round(ForcedUptake[0][r],3), "\n"

    return ForcedUptake

def ForceAAuptakeUtilAA(m, PropBF=0.0, AMaint=45.0, AA=["GLY_AA_mm_tx"],
        blocks=["GLN_AA_mm_tx"], NTxs=DataForAnalysis.NTxs_MM, rm=[]):
## This fx calculates the overall uptake or excretion of AAs
    ## upon forcing uptake of a certain amino acid/s (AA)
        ## for BMP upon standar conditions in MMc medium (19 AA, NH4 avail)
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium

    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
        post: returns a dictionary with AAs in media as keys and values are
        total amount of AA uptaken/excreted during BMP (mmol/gram cDW/hour)"""

    ForceAA = ForceAAuptake(m, PropBF, AMaint, AA, blocks)
    fd = GetFluxDic(m, PropBF, AMaint)  # fd for BMP (1 unit) and AMaint of m

    AAavail = NTxs[:]
    
    for comp in rm:
        AAavail.remove(comp)

    AAuptake = {}
    for aa in AAavail:
        if aa not in ForceAA[0].keys():
            AAuptake[aa] = 0.0
        else:
            AAuptake[aa] = ForceAA[0][aa]   # this fx forces import of GLY

    AAuptake2 = {}
    for aa in AAuptake.keys():
        if '_AA_mm_tx' in aa:
            AAuptake2[aa.replace('_AA_mm_tx','')] = AAuptake[aa]
        else:
            if '_mm_tx' in aa:  # e.g. for comps in NTxs in synovial fluid (no AA_mm_tx)
                AAuptake2[aa.replace('_mm_tx','')] = AAuptake[aa]

    AAexcret = {}
    for aa in filter(lambda s: "_AA_bm_tx" in s, m.smx.cnames):
        if aa not in ForceAA[0].keys():
            AAexcret[aa] = 0.0
        else:
            AAexcret[aa] = ForceAA[0][aa]

    AAexcret2 = {}
    for aa in AAexcret.keys():
        AAexcret2[aa.replace('_AA_bm_tx','')] = AAexcret[aa]

    AAtotalFlux = {}
    
    NewKeys = []
    for aa in AAavail:
        if '_AA_mm_tx' in aa:
            NewKeys.append(aa.replace('_AA_mm_tx',''))
        else:
            if '_mm_tx' in aa:
                NewKeys.append(aa.replace('_mm_tx',''))

    for aa in NewKeys:
        if aa in AAexcret2.keys():
            AAtotalFlux[aa] = AAuptake2[aa] + AAexcret2[aa]
        else:
            AAtotalFlux[aa] = AAuptake2[aa]

    print "Total amounts of AA uptaken/excreted during BMP (mmol/gram cDW/hour): ", "\n"
    for aa in sorted(AAtotalFlux):
        print aa, ":  ", round(AAtotalFlux[aa],3) 

    return AAtotalFlux, ForceAA

def Nassimilation(m, prods =["GLT_AA_bm_tx"], AA=0.0, blocks=[]):
## This fx investigates N assimilation from NH4 and solves a LP for GLT production
    ## with NH4 as sol N source and GLC as sole C source
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium
    """pre: m = the model
            prods = compounds to be produced and excreted (1mmol)
            AA = flux through AA importers
            blocks = reactions to be blocked
        post: lp solution for production and excretion of 1mmol of AAs (GLT)
        with NH4 as sole N source and GLC as sole C source"""

    lp = BuildLP.BuildLP(m)

    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})

    for prod in prods:
        lp.SetFixedFlux({prod:-1})

    for r in blocks:
        lp.SetFixedFlux({r:0.0})

    lp.Solve(False)
    sol = lp.GetPrimSol()
    
    print "These are the compounds imported/exported in the solution: ", "\n"
    for r in sorted(sol):
        if "_tx" in r:
            print round(sol[r]), r
    print "\n", "These are the reactions present in the solution for N assimilation from NH4 and production of 1 mmol of GLT: ", "\n"
    print len(sol), "reactions", "\n"
    for r in sorted(sol):
        print round(sol[r],2), "\n", m.smx.ReacToStr(r)

    return sol

##############
    
def ContribNuptakeAAs(m, PropBF=0.0, AMaint=45.0, AA=None,
    NTxs=DataForAnalysis.NTxs_MM, rm=[], blocks=["GLN_AA_mm_tx"],
    Natoms=DataForAnalysis.Natoms_MM_AA):
## ratio uptake vs demand and contribution to total N uptake
    # for 1 unit BMP upon standar conditions in MMc medium(19 AA and NH4 available)
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium

    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
            Natoms = n atoms in AAs in medium
        post: returns dict mapping AA txs to their ratio of uptake/demand
        and proportional contribution to total N uptake"""

    lp = GetLP(m, PropBF, AMaint)
    
    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})

    for r in blocks:
        lp.SetFixedFlux({r:0.0})

    for comp in rm: # rm = list of N sources not available
        NTxs.remove(comp)

    lp.Solve(False)
    sol = lp.GetPrimSol()

    AAReqs = AARequirement(m, PropBF, AMaint, NTxs, rm)         # dic of fluxes required in the AA biomass output
    
    if "NO2_mm_tx" in sol.keys() and "NH4_mm_tx" in sol.keys(): # only comps containing N that can be excreted appart from AAs in BM are NO2 and NH4 as by-products.
            TotalN = TotalNDemmand(m, PropBF, AMaint, NTxs, rm, blocks) - sol["NH4_mm_tx"] - sol["NO2_mm_tx"]    # Now, total N consumed is:
                                                        # Demand calculated from NH4 when sole N source (= N excreted = BMP)
                                                        # plus the N exctreted as NH4 if AAs avail! -- = +
                                                        # so N in BMP plus N excreted as any other byprod (NH4, NO2 etc.).
    if "NO2_mm_tx" in sol.keys() and "NH4_mm_tx" not in sol.keys():
            TotalN = TotalNDemmand(m, PropBF, AMaint, NTxs, rm, blocks) - sol["NO2_mm_tx"]
            
    if "NH4_mm_tx" in sol.keys() and "NO2_mm_tx" not in sol.keys():
            TotalN = TotalNDemmand(m, PropBF, AMaint, NTxs, rm, blocks) - sol["NH4_mm_tx"]
    else:
        TotalN = TotalNDemmand(m, PropBF, AMaint, NTxs, rm, blocks)
                                                            

    NperAA = {}     # N uptaken as AAs
    for aa in Natoms.keys():
        if aa in sol.keys():
            NperAA[aa] = Natoms[aa]*sol[aa]
    
    rv = {}
    for aa in AAReqs.keys():           # for those aa required in the AA biomass output  
        if aa not in sol.keys(): # if AA not uptaken
            propUD = propT = 0.0
        else:
            if AAReqs[aa] == 0.0:
                propUD = 0.0  # ratio uptake/demand = flux in sol (uptaken) / flux for bm output (demand for BMP)
                propT = NperAA[aa]/TotalN   # prop contrib to total N uptake = flux in sol*number N atoms / total N demand
                                            # and / total N demand (from NH4)
            else:
                propUD = sol[aa]/AAReqs[aa]  # ratio uptake/demand = flux in sol (uptaken) / flux for bm output (demand for BMP)
                propT = NperAA[aa]/TotalN   # prop contrib to total N uptake = flux in sol*number N atoms / total N demand
                                            # and / total N demand (from NH4) 
        rv[aa] = propUD, propT

    propTdic = {}
    for aa in AAReqs:
        if aa not in sol.keys(): # if AA not uptaken
            propTdic[aa] = 0.0
        else:
            propTdic[aa] = NperAA[aa]/TotalN
    Totalpercent = round((sum(propTdic.values())*100),3)

    print "\n", "These are the AA uptake/demand ratios for BMP (1 gram): ", "\n"
    for aa in sorted(AAReqs.keys()):
        print aa, "","","","", round(rv[aa][0],3)
    print "\n", "This is the total N demand under these conditions (mmol): ", "\n", round(TotalN,3)    
    print "\n", "These are the AA contributions to the total N uptaken (percentage): ", "\n"
    for aa in sorted(AAReqs.keys()):
        print aa, "","","","", round(rv[aa][1]*100,3)

    print "\n", "Sum: ", Totalpercent
    
    return rv

def ContribNuptakeALL(m, PropBF=0.0, AMaint=45.0, AA=None,
    NTxs=DataForAnalysis.NTxs_synovial, rm=[], blocks=["GLN_AA_mm_tx"],
    Natoms=DataForAnalysis.Natoms_synovial_noNH4, FreeBM=None):
## ratio uptake vs demand and contribution to total N uptake
    # for 1 unit BMP upon standar conditions in MMc medium(19 AA and NH4 available)
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium

    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
            blocks = reactions to be blocked
            Natoms = n atoms in N containing comps in medium excluding NH4 (since it is excreted)
            FreeBM = if != None, upper flux bound of BM_tx (biomass txs) is unbound
        post: returns dict mapping AA txs to their ratio of uptake/demand
        and proportional contribution to total N uptake"""

    if FreeBM != None:
        lp = GetLPfreeBM(m, PropBF, AMaint)
    else:
        lp = GetLP(m, PropBF, AMaint)

    if AA != None:
        for aa in filter(lambda s: "_AA_mm_tx" in s, m.sm.cnames):
            lp.SetFixedFlux({aa:AA})
    for r in blocks:
        lp.SetFixedFlux({r:0.0})
##    Block = Set.Complement(NTxs, Avail)             # For AA in NTxs not in Avail (so we can block their uptake)
##    BlockDic = dict(zip(Block, [0] * len(Block)))   # flux dictionary to block uptake of AA if not available
##    lp.SetFixedFlux(BlockDic)                       # blocks flux
    lp.Solve(False)
    sol = lp.GetPrimSol()

    for comp in rm: # rm = list of N sources not available
        NTxs.remove(comp)

    AAReqs = AARequirement(m, PropBF, AMaint, NTxs, rm)      # dic of fluxes required in the AA biomass output

    if "NO2_mm_tx" in sol.keys() and "NH4_mm_tx" in sol.keys(): # only comps containing N that can be excreted appart from AAs in BM are NO2 and NH4 as by-products.
            TotalN = TotalNDemmand(m, PropBF, AMaint, NTxs, rm, blocks) - sol["NH4_mm_tx"] - sol["NO2_mm_tx"]    # Now, total N consumed is:
                                                        # Demand calculated from NH4 when sole N source (= N excreted = BMP)
                                                        # plus the N exctreted as NH4 if AAs avail! -- = +
                                                        # so N in BMP plus N excreted as any other byprod (NH4, NO2 etc.).
    if "NO2_mm_tx" in sol.keys() and "NH4_mm_tx" not in sol.keys():
            TotalN = TotalNDemmand(m, PropBF, AMaint, NTxs, rm, blocks) - sol["NO2_mm_tx"]
            
    if "NH4_mm_tx" in sol.keys() and "NO2_mm_tx" not in sol.keys():
            TotalN = TotalNDemmand(m, PropBF, AMaint, NTxs, rm, blocks) - sol["NH4_mm_tx"]
    else:
        TotalN = TotalNDemmand(m, PropBF, AMaint, NTxs, rm, blocks)

    NperComp = {}     # N uptaken as AAs
    for comp in Natoms.keys():
        if comp in sol.keys():
            NperComp[comp] = Natoms[comp]*sol[comp]
    
    rv = {}
    for comp in Natoms:           # for those comps with N available
        if comp not in sol.keys():  # if AA not uptaken
            propUD = propT =0.0
        else:
            if comp not in AAReqs.keys():
                propUD = 0.0                    # ratio uptake/demand = flux in sol (uptaken) / flux for bm output (demand for BMP)
                propT = NperComp[comp]/TotalN   # prop contrib to total N uptake = flux in sol*number N atoms / total N uptaken
            else:
                if AAReqs[comp] == 0.0:         # if comp is not part of the BM (zero demand)
                    propUD = 0.0              
                    propT = NperComp[comp]/TotalN  
                else:
                    propUD = sol[comp]/AAReqs[comp]     # ratio uptake/demand = flux in sol (uptaken) / flux for bm output (demand for BMP)
                    propT = NperComp[comp]/TotalN    
        rv[comp] = propUD, propT

    propTdic = {}
    for comp in Natoms:
        if comp not in sol.keys():  # if comp in medium is not uptaken
            propTdic[comp] = 0.0
        else:
            propTdic[comp] = NperComp[comp]/TotalN

    Totalpercent = round((sum(propTdic.values())*100),3)

    print "\n", "These are the AA uptake/demand ratios for BMP (1 gram): ", "\n"
    for aa in sorted(AAReqs.keys()):
        print aa, "","","","", round(rv[aa][0],3)
    print "\n", "This is the total N demand under these conditions (mmol): ", "\n", round(TotalN,3)    
    print "\n", "These are the contributions of medium components to total N uptake (percentage): ", "\n"
    for comp in sorted(Natoms.keys()):
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

def TotalNuptakeALL(m, PropBF=0.0, AMaint=45.0, AA=None,
        blocks=["GLN_AA_mm_tx"], Natoms=DataForAnalysis.Natoms_synovial_noNH4,
        FreeBM=None):
## This function calculates total N uptaken for BMP (1 gram)
   ## and ATP maintenance upon standar conditions, MMc medium
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            Natoms = n atoms in N containing comps in medium excluding NH4 (since it is excreted)
            FreeBM = if != None, upper flux bound of BM_tx (biomass txs) is unbound
        post: returns solution for BMP in standar MM (GLC and AAs present);
        Plus values for N uptaken"""
    
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

    NperComp = {}     # N uptaken as comps from medium
    for comp in Natoms.keys():
        if comp in sol.keys():
            NperComp[comp] = Natoms[comp]*sol[comp]
        else:
            NperComp[comp] = 0.0

    print "\n", "These are the mmol of N uptaken as medium compounds for BMP: ", "\n"
    for comp in sorted(NperComp.keys()):
        print comp, "","","","", round(NperComp[comp],3)

    TotalNuptaken = round(sum(NperComp.values()),3)
    print "\n", "This is the total N uptaken (mmol): ", TotalNuptaken
    return TotalNuptaken

####################################################################### ENDS FREE BM 

def AARequirement(m, PropBF=0.0, AMaint=45.0, NTxs=DataForAnalysis.NTxs_MM, rm=[]):
## This function matches the fluxes through AA importers available to the flux   
    ## requiered through the AA exporters for BMP of 1 gram of planktonic biomass
    """pre: mm = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
        post: returns dict mapping AA importers in media to the flux
        through the AA exporters for production of BM (1 gram)"""
    rv = {}
    fd = GetFluxDic(m, PropBF, AMaint)

    NInputs = NTxs[:]   # list of all possible N sources (media txs)

    NOutputs = map(lambda a: a.replace("_AA_mm_","_AA_bm_"), NTxs) # replace in the list the _AA_mm_ suffix by _bm_
    Out2In = dict(zip(NOutputs, NInputs))                          # dic mapping txs:  AA_bm_tx to _AA_mm_tx
    UsedOutputs = Set.Intersect(NOutputs, fd)                   # list with those AA_bm_tx generated from NTxs
                                                                # that are also in fd for BMP
    
    for out in UsedOutputs:             # for AA_bm_tx in fd for which AA_mm_tx in media (comps in medium with bm_tx in model)
        rv[Out2In[out]] = -fd[out]      # rv is a dic with their AA_mm_tx as keys
                                        # and flux is same flux as for AA_bm_tx for BMP (but positive)
    return rv

def NexportBMasAA(sol, NatomsBM=DataForAnalysis.Natoms_BM_AA):
## This function calculates total amount of N excreted as amino acids
    ## for BMP of 1 gram of cDW (per hour)
    """pre: sol = LP solution
            Natoms = n atoms of N in AAs in BM!!!
        post: returns total amount of N excreted as amino acids (mmol)
        during BMP per gram cDW/h"""


    NperAAbm = {}     # N excreted as AAs
    for aa in NatomsBM.keys():
        if aa in sol.keys():
            NperAAbm[aa] = NatomsBM[aa]*sol[aa]

    TotalNescretedAsAA = sum(NperAAbm.values())
            
    print "\n", "This is the total N excreted as amino acids during BMP (per 1 gram cDW): ", round(TotalNescretedAsAA,3) 

    return TotalNescretedAsAA

def TotalNDemmand(m, PropBF=0.0, AMaint=45.0, NTxs=DataForAnalysis.NTxs_MM, rm=[],
                  blocks=["GLN_AA_mm_tx"]): 
## This function calculates total N demand as NH4 imported for BMP (1 gram) and ATP maintenance
    ## ONLY true if no other N sources are available
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
            blocks = reactions to be blocked
        post: returns flux value for NH4 tx in solution for BMP were
        NH4 is the only available N source"""
    
    lp = GetLP(m, PropBF, AMaint)

    for r in blocks:
        lp.SetFixedFlux({r:0.0})

    block = NTxs[:]
    block.remove("NH4_mm_tx")
    blockdic = dict(zip(block, [0]*len(block))) 
    lp.SetFixedFlux(blockdic)
    lp.Solve(False)
    res = lp.GetPrimSol()
    NdemandNH4 = res["NH4_mm_tx"]
    print "\n", "This is the total N demand when NH4 is the only N source available (mmol): ", round(NdemandNH4,3)
    # Only comps excreted containing N are AAs so N uptake as NH4 = total N amount in AAs  in the biomass (steady state)
    return NdemandNH4

def TotalGLCDemmand(m, PropBF=0.0, AMaint=45.0, NTxs=DataForAnalysis.NTxs_MM, rm=[],
                    blocks=["GLN_AA_mm_tx"]):
## This function calculates total GLC demand for BMP and ATP maintenace
    ## ONLY true if no other C sources are available
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
            blocks = reactions to be blocked
        post: returns flux value for GLC importer in solution for BMP were
        GLC is the only available C source"""
    
    lp = GetLP(m, PropBF, AMaint)
    
    for r in blocks:
        lp.SetFixedFlux({r:0.0})

    block = NTxs[:]

    block.remove("NH4_mm_tx")   # Needed as N source

    blockdic = dict(zip(block, [0]*len(block))) 
    lp.SetFixedFlux(blockdic)
    lp.Solve(False)
    sol = lp.GetPrimSol()
    GLCdemand = sol["GLC_mm_tx"]

    print "\n", "This is the total GLC demand when GLC is the only C source available (mmol): ", round(GLCdemand,3) 
    
    return GLCdemand

def TotalCDemmand(m, PropBF=0.0, AMaint=45.0, NTxs=DataForAnalysis.NTxs_MM, rm=[],
                    blocks=["GLN_AA_mm_tx"]):
## This function calculates total C demand for BMP (1 gram) from GLC
    ## ONLY true if no other C sources are available
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
            blocks = reactions to be blocked
        post: returns flux value for GLC importer in solution for BMP were
        GLC is the only available C source multiplied by 6 = C demand"""

    GLCdemand=TotalGLCDemmand(m, PropBF, AMaint, NTxs, rm, blocks)

    Cdemand = (GLCdemand*6)
    
    print "\n", "This is the total C demand when GLC is the only C source available (mmol): ", round(Cdemand,3)

    return Cdemand

#################################

def CompareSols(m, PropBF=0.0, AMaint=45.0, AA=["GLY_AA_mm_tx"],
                blocks=["GLN_AA_mm_tx"],NTxs=DataForAnalysis.NTxs_MM, rm=[]):
## This fx compares solutions when flux through AA importer/s is unbound or fixed
        ## for 1 unit BMP upon standar conditions in MMc medium (19 AA, NH4 avail)
            # "GLN_AA_mm_tx" is blocked because is not present in MMc medium

    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
        post: print reactions of interest with to compare their fluxes between two solutions"""

    sol1 = ff = ForceAAuptakeUtilAA(m, PropBF, AMaint, AA, blocks, NTxs, rm)
    sol2 = st = UtilAA(m, PropBF, AMaint, blocks, NTxs, rm)[1]    # = AAessential.StandarBMPinMM(m) 

    print "\n", "Reactions in standar solution not present when flux through AA importer/s is forced: ", "\n"
    for r in st[0]:
        if r not in ff[1][0]:
	    print round(st[0][r],3), "\n", m.smx.ReacToStr(r)

    print "\n", "Reactions in new solution (AA importer/s flux forced) not present in standar solution: ", "\n"
    for r in ff[1][0]:
        if r not in st[0]:
	    print round(ff[1][0][r],3), "\n", m.smx.ReacToStr(r)

    print "\n", "Reactions involved with ATP present in both solutions which flux differs in > 1e-3 mmol: ", "\n"
    for r in m.sm.InvolvedWith("ATP"):
        if r in ff[1][0]:
            if r in st[0]:
                if (abs(ff[1][0][r]) - abs(st[0][r]) ) > 1e-3:
                    print round(ff[1][0][r],3),"","", "Reaction when flux is forced", "\n", m.smx.ReacToStr(r)
                    print round(st[0][r],3),"","", "Reaction in standar solution", "\n", m.smx.ReacToStr(r)
            else:
                pass
        else:
            pass

    print "\n", "Reactions involved with NADH present in both solutions which flux differs in > 1e-3 mmol: ", "\n"
    for r in m.sm.InvolvedWith("NADH"):
        if r in ff[1][0]:
            if r in st[0]:
                if (abs(ff[1][0][r]) - abs(st[0][r]) ) > 1e-3:
                    print round(ff[1][0][r],3),"","", "Reaction when flux is forced", "\n", m.smx.ReacToStr(r)
                    print round(st[0][r],3),"","", "Reaction in standar solution", "\n", m.smx.ReacToStr(r)
            else:
                pass
        else:
            pass

    print "\n", "Reactions involved with NADPH present in both solutions which flux differs in > 1e-3 mmol: ", "\n"
    for r in m.sm.InvolvedWith("NADPH"):
        if r in ff[1][0]:
            if r in st[0]:
                if (abs(ff[1][0][r]) - abs(st[0][r]) ) > 1e-3:
                    print round(ff[1][0][r],3),"","", "Reaction when flux is forced", "\n", m.smx.ReacToStr(r)
                    print round(st[0][r],3),"","", "Reaction in standar solution", "\n", m.smx.ReacToStr(r)
            else:
                pass
        else:
            pass

####################################
