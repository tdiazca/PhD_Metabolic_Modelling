######
########
##########   This module uses LPA for identification of essential AA
##########      and study the effect of amino acid deprivation upon ObjVal and GLC consumption
########
#######    NOTE: remember, "GLN_AA_mm_tx", "NO3_mm_tx" are blocked 
######              cause are not present in MMc media
####                    they are available in synovial fluid so do not block in this case
##          NOTE: select adequate NTxs list (uncomment list of N containing compounds depending on medium)

import sys
# sys.path.append('../Analysis/AAutil/')
from Util import Set
import BuildLP
import DataForAnalysis
from Bioinf import Biomass
import SepiBiomass
import math

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

def StandarBMPinMM(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"]):
## This fx calculates objective function and GLC consumption values
    # for BMP and ATP maintenace upon standar conditions in MMc medium(19 AA and NH4 available)
            # "GLN_AA_mm_tx", "NO3_mm_tx" are blocked cause are not present in MMc medium

    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
        post: returns a tuple with [0] = solution to FBA for BMP in standar media (dic),
        [1] = its objective value, [2] = its glucose consumption"""

    StdBMP = ()  # tuple
    lp = GetLP(m, PropBF, AMaint)
    for r in blocks:
        lp.SetFixedFlux({r:0.0})
    lp.Solve(False)
    sol = lp.GetPrimSol()
    ObjVal = lp.GetObjVal()
    GLCcomp = sol["GLC_mm_tx"]
    StdBMP= sol, ObjVal, GLCcomp # return ObjVal and lp sol

    return StdBMP

def PrintValsStdBMP(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"]):
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
        post: prints main parameters for sol for BMP in standar medium"""

    StdBMP = StandarBMPinMM(m, PropBF, AMaint, blocks)

    print "This is the objective value (mmol/gDW/h) for BMP in MM (19 AA plus NH4 available): ", round(StdBMP[1],3), "\n"
    print "This is the number of reactions in the solution for BMP in MM (19 AA plus NH4 available): ", len(StdBMP[0]), "\n"
    print "This is the GLC consumed (mmol) for BMP in MM (19 AA plus NH4 available): ", round(StdBMP[2],3),"\n"
    print "This is the NH4 imported/excreted (mmol): ", round(StdBMP[0]["NH4_mm_tx"],3), "\n"
    print "These are by-products excreted (mmol): " ,"\n"
    for r in sorted(StdBMP[0]):
        if r in filter(lambda s: "_bp_tx" in s, m.smx.cnames):
            print r, round(StdBMP[0][r],3), "\n"

#############################################################  FREE BM (no upper flux bounds)
## Analyse system for production of 1 unit of cell biomass (BMP)
    # while meeting ATP cell maintenance requirements
        # and allowing BM prods to be exported at a higher proportion
            # than needed to produce BM!

def GetLPfreeBM(m, PropBF=0.0, AMaint=45.0):
    """pre: GetFluxDic(m, PropBF=0.0, AMaint=45.0)
            m = the model
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

def StandarBMPinMMfreeBM(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"]):
## This fx calculates objective function and GLC consumption values
    # for BMP and ATP maintenace upon standar conditions in MMc medium(19 AA and NH4 available)
        # relaxed upper flux bounds for BM tx reactions
            # "GLN_AA_mm_tx", "NO3_mm_tx" are blocked cause are not present in MMc medium
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
        post: returns a tuple with [0] = solution to FBA for BMP in standar media (dic),
        [1] = its objective value, [2] = its glucose consumption ; when BM tx have relaxed
        upper flux bounds"""

    StdBMPfreeBM = ()  # tuple
    lp = GetLPfreeBM(m, PropBF, AMaint)
    #lp = GetLP(m, PropBF, AMaint)
    for r in blocks:
        lp.SetFixedFlux({r:0.0})
    lp.Solve(False)
    sol = lp.GetPrimSol()
    ObjVal = lp.GetObjVal()
    GLCcomp = sol["GLC_mm_tx"]
    StdBMPfreeBM= sol, ObjVal, GLCcomp # return ObjVal and lp sol

    return StdBMPfreeBM

def PrintValsStdBMPfreeBM(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"]):
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
        post: pprints main parameters for sol for BMP in standar medium when BM tx have relaxed
        upper flux bounds"""

    StdBMPfreeBM = StandarBMPinMMfreeBM(m, PropBF, AMaint, blocks)

    print "This is the objective value (mmol/gDW/h) for BMP in MM (19 AA plus NH4 available): ", round(StdBMPfreeBM[1],3), "\n"
    print "This is the number of reactions in the solution for BMP in MM (19 AA plus NH4 available): ", len(StdBMPfreeBM[0]), "\n"
    print "This is the GLC consumed (mmol) for BMP in MM (19 AA plus NH4 available): ", round(StdBMPfreeBM[2],3),"\n"
    print "This is the NH4 imported/excreted (mmol): ", round(StdBMPfreeBM[0]["NH4_mm_tx"],3), "\n"
    print "These are by-products excreted (mmol): " ,"\n"
    for r in sorted(StdBMPfreeBM[0]):
        if r in filter(lambda s: "_bp_tx" in s, m.smx.cnames):
            print r, round(StdBMPfreeBM[0][r],3), "\n"

##################################### ENDS FREE BM

def EssentialAA(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"],
                NTxs=DataForAnalysis.NTxs_MM, rm=[]):
## Essentiality of single amino acids/NH4 for BMP
    # while meeting the E demand for cell maintenance (NGAM + GAM ATP demand) 
        # "GLN_AA_mm_tx" is blocked cause is not present in MMc media

    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
        post: returns a dic with single AA transporter names as keys and either essentiality
        or feasible sol as values; returns a list of essential and a list of non essential AA for BMP"""

    EssAA = {}

    lp = GetLP(m, PropBF, AMaint)

    for r in blocks:
        lp.SetFixedFlux({r:0.0})

##    targets = NTxs
##    targets.remove('NH4_mm_tx')
##    lp.SetFluxBounds({'NH4_mm_tx':(None,0)}) # NH4 export is allowed

    for comp in rm:
        NTxs.remove(comp)

    for n in NTxs: # we are just blocking one at a time so NH4 export is allowed
        if n in m.smx.cnames:
            lp.SetFixedFlux({n:0.0})                              
            lp.Solve(False)     
            if lp.IsStatusOptimal():
                EssAA[n] = lp.GetPrimSol() # return ObjVal and lp sol
            else:
                EssAA[n] = "Essential"
            if n in Set.Complement(NTxs, blocks):   # so for blocks, constraint is not removed!
                lp.ClearFluxConstraint(n)

    EssentialAAs = []
    for n in EssAA.keys():
        if EssAA[n] == "Essential":
            EssentialAAs.append(n)

    NonEssentialAAs = []
    for n in EssAA.keys():
        if EssAA[n] != "Essential":
            NonEssentialAAs.append(n)
    
    return EssAA, EssentialAAs, NonEssentialAAs

def PrintListEssentialAAs(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"],
                          NTxs=DataForAnalysis.NTxs_MM, rm=[]):
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
        post: prints a list of essential and a list of non essential AA for BMP"""

    DefineEssentialAAs = EssentialAA(m, PropBF, AMaint, blocks, NTxs, rm)

    EssentialAAs = DefineEssentialAAs[1]

    NonEssentialAAs =  DefineEssentialAAs[2]
    
    print "These AA/N sources are essential for biomass production:", "\n", sorted(EssentialAAs), "\n"
    print "These AA/N sources are not essential for biomass production:", "\n", sorted(NonEssentialAAs), "\n"

def BMPonlyNH4(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"],
               NTxs=DataForAnalysis.NTxs_MM, rm=[]):
## This fx calculates objective function and GLC consumption values
    # for BMP in MM medium with NH4 as sole N source
            # "GLN_AA_mm_tx", "NO3_mm_tx" are blocked cause are not present in MMc medium

    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTx
        post: returns a tuple with [0] = solution to FBA (dic), [1] = its objective value,
        [2] = its glucose consumption"""

    NonEssentialAAs = EssentialAA(m, PropBF, AMaint, blocks, NTxs, rm)[2]

    BMPnh4 = ()  # tuple
    
    lp = GetLP(m, PropBF, AMaint)
    
    for r in blocks:
        lp.SetFixedFlux({r:0.0})

    targets = NonEssentialAAs # because all AAs are non-essential
    targets.remove('NH4_mm_tx')
    for n in targets:
        lp.SetFixedFlux({n:0.0})
    lp.Solve(False)
    sol = lp.GetPrimSol()
    ObjVal = lp.GetObjVal()
    GLCcomp = sol["GLC_mm_tx"]
    BMPnh4 = sol, ObjVal, GLCcomp # return ObjVal and lp sol

    return BMPnh4

def PrintValsBMPinNH4(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"],
    NTxs=DataForAnalysis.NTxs_MM, rm=[]):
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
        post: prints main parameters for sol for BMP in NH4 as sole N
        source"""

    BMPnh4 = BMPonlyNH4(m, PropBF, AMaint, blocks, NTxs, rm)

    print "This is the objective value (mmol/gDW/h) for BMP in MM with NH4 as sole N source: ", round(BMPnh4[1],3), "\n"
    print "This is the number of reactions in the solution for BMP in MM with NH4 as sole N source: ", len(BMPnh4[0]), "\n"
    print "This is the GLC consumed (mmol) for BMP in MM with NH4 as sole N source: ", round(BMPnh4[2],3),"\n"
    print "This is the NH4 imported/excreted (mmol): ", round(BMPnh4[0]["NH4_mm_tx"],3), "\n"
    print "These are by-products excreted (mmol): " ,"\n"
    for r in sorted(BMPnh4[0]):
        if r in filter(lambda s: "_bp_tx" in s, m.smx.cnames):
            print r, round(BMPnh4[0][r],3), "\n"

def BMPwithAAselection(m, PropBF=0.0, AMaint=45.0,
    SetAAavail=["GLT_AA_mm_tx","ALA_AA_mm_tx","ARG_AA_mm_tx"], blocks=["GLN_AA_mm_tx"],
    NTxs=DataForAnalysis.NTxs_MM, rm=[]):
## This fx checks for feasibility of BMP with a selection of AAs as N sources
        # "GLN_AA_mm_tx", "NO3_mm_tx" are blocked by default cause are not present in MMc medium

    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            SetAAavail = list of AA available in media
            blocks = reactions which flux we want to block
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
        post: if feasible, returns solution for BMP under these conditions"""

    lp = GetLP(m, PropBF, AMaint)
    for r in blocks:
        lp.SetFixedFlux({r:0.0})

    for comp in rm:
        NTxs.remove(comp)

    block = NTxs[:]                                 # NTxs = all N sources in MM media
    block.remove("NH4_mm_tx")                       # block is NTxs list with NH4 removed
    BlockDic = dict(zip(block, [0.0]*len(block)))   # make dir for N tx with 0 as values
    lp.SetFixedFlux(BlockDic)
    lp.SetFluxBounds({"NH4_mm_tx":(None,0.0)})

    Unlock = Set.Intersect(NTxs, SetAAavail)        # For AA in NTxs and SetAAavail
    for aa in Unlock:
        lp.ClearFluxConstraint(aa)
    lp.Solve(False)

    if lp.IsStatusOptimal():
        sol = lp.GetPrimSol()
        LenSol = len(sol)
        ObjVal = lp.GetObjVal()
        GLCcon = sol["GLC_mm_tx"]
        if "NH4_mm_tx" in sol.keys():
            NH4exp = sol["NH4_mm_tx"]
        else:
            NH4exp = 0.0
        return sol, LenSol, ObjVal, GLCcon, NH4exp 

    else:
        print "BMP not possible"

def PrintValsBMPinAAselect(m, PropBF=0.0, AMaint=45.0,
    SetAAavail=["GLT_AA_mm_tx","ALA_AA_mm_tx","ARG_AA_mm_tx"], blocks=["GLN_AA_mm_tx"],
    NTxs=DataForAnalysis.NTxs_MM, rm=[]):
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            SetAAavail = list of AA available in media
            blocks = reactions which flux we want to block
        post: if feasible, print values for BMP under these conditions"""

    sol = BMPwithAAselection(m, PropBF, AMaint, SetAAavail, blocks, NTxs, rm)[0]
    LenSol = len(sol)
    ObjVal = BMPwithAAselection(m, PropBF, AMaint, SetAAavail, blocks, NTxs, rm)[2]
    GLCcon = sol["GLC_mm_tx"]

    if "NH4_mm_tx" in sol.keys():
        NH4exp = sol["NH4_mm_tx"]
    else:
        NH4exp = 0.0

    print "\n", "This is the lenth of the solution: ", LenSol,
    print "\n", "This is the objective value of the solution (mmol/gDW/h): ", round(ObjVal,3),
    print "\n", "This is the amount of glucose consumed (mmol): ",round(GLCcon,3),
    print "\n", "This is the amount of NH4 exported (mmol): ",round(NH4exp,3),
    print "\n", "These are the amounts of AA consumed (mmol): ", "\n"
    for comp in sorted(BuildLP.AAmediatx(m)):
        if comp in sol.keys():
            print round(sol[comp],3), comp
        else:
            pass
    print "\n", "These are the amounts of by-products produced (mmol): ", "\n"
    for comp in sorted(BuildLP.BPtx(m)):
        if comp in sol.keys():
            print round(sol[comp],3), comp
        else:
            pass

def AArmBMP(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"],
            NTxs=DataForAnalysis.NTxs_MM, rm=[]):
## This fx calculates solutions for BMP upon removal of single AA/NH4
    ## prints objective values and flux through GLC importer
            ## "GLN_AA_mm_tx", "NO3_mm_tx" are blocked by default
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
        post: returns a dictionary with AA/NH4 as keys and values are tuples for solutions
        upon rm of each individual AA/NH4 with [0] = LPA sol, [1] = the objective value,
        [2] = the glucose consumption"""

    NonEssentialAAs = EssentialAA(m, PropBF, AMaint, blocks, NTxs, rm)[2]

    AArm = {}

    lp = GetLP(m, PropBF, AMaint)

    for r in blocks:
        lp.SetFixedFlux({r:0.0})

    for n in NonEssentialAAs:   # we are just blocking one at a time so NH4 export is allowed
        if n == 'NH4_mm_tx':
            lp.SetFluxBounds({"NH4_mm_tx":(None,0.0)}) # NH4 not available for import but can still be exported
            lp.Solve(False)
            if lp.IsStatusOptimal():
                sol = lp.GetPrimSol()
                ObjVal = lp.GetObjVal()
                GLCcomp = sol["GLC_mm_tx"]
                AArm[n] = sol, ObjVal, GLCcomp # return ObjVal and lp sol
            if n in Set.Complement(NonEssentialAAs, blocks):
                lp.ClearFluxConstraint(n)
        else:
            lp.SetFixedFlux({n:0.0})                              
            lp.Solve(False)
            if lp.IsStatusOptimal():
                sol = lp.GetPrimSol()
                ObjVal = lp.GetObjVal()
                GLCcomp = sol["GLC_mm_tx"]
                AArm[n] = sol, ObjVal, GLCcomp # return ObjVal and lp sol
            if n in Set.Complement(NonEssentialAAs, blocks):
                lp.ClearFluxConstraint(n)

    return AArm

def PrintValsUponAArm(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"]):
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
        post: prints the objective values and glucose consumption for BMP upon removal
        of single AA/NH4"""

    AArm = AArmBMP(m, PropBF, AMaint, blocks)

    print "These are the objective values (mmol/gDW/h) upon removal of single AA/NH4 for BMP:", "\n"
    for n in sorted(AArm.keys()):
        print n, "","","","", round(AArm[n][1],3)
    print "\n", "These are the mmol of glucose consumed upon removal of single AA/NH4 for BMP:", "\n"
    for n in sorted(AArm.keys()):
        print n, "","","","", round(AArm[n][2],3)
    
def EffectAArmBMP(m, PropBF=0.0, AMaint=45.0, blocks=["GLN_AA_mm_tx"],
    NTxs=DataForAnalysis.NTxs_MM, rm=[]):
## This fx calculates the effect of single AA/NH4 rm upon BMP 
    ## in objective value and GLC consumption
            # "GLN_AA_mm_tx", "NO3_mm_tx" are blocked cause are not present in MMc media
    """pre: m = the model
            PropBF = biofilm proportion in biomass composition
            AMaint = ATP maintenance demand (NGAM plus GAM), growth rate 1gDW/h = 45 mmol/gDW
            blocks = reactions to be blocked
            NTxs = list of tx of N sources available in media
            rm = tx reacs for N sources not available to remove from NTxs
        post: returns dicts with changes in objective value and GLC consumption
        expressed as percentages and has percentage of variation"""

    NonEssentialAAs = EssentialAA(m, PropBF, AMaint, blocks, NTxs, rm)[2]
    AArm = AArmBMP(m, PropBF, AMaint, blocks)
    ObjValStdBMP = StandarBMPinMM(m, PropBF, AMaint, blocks)[1]
    glcStdBMP = StandarBMPinMM(m, PropBF, AMaint, blocks)[2]
    ObjValBMPnh4 = BMPonlyNH4(m, PropBF, AMaint, blocks, NTxs, rm)[1]
    glcBMPnh4 = BMPonlyNH4(m, PropBF, AMaint, blocks, NTxs, rm)[2] 

    effObjValperc = {}

    effGLCperc = {}

    effObjValvar = {}

    effGLCvar = {}

    effObjValpercNH4 = {}

    effGLCpercNH4 = {}

    effObjValvarNH4 = {}

    effGLCvarNH4 = {}

    for n in NonEssentialAAs:
        effObjValperc[n] = (AArm[n][1] *100) / ObjValStdBMP

    for n in NonEssentialAAs:
        effObjValvar[n] =  effObjValperc[n] - 100

    for n in NonEssentialAAs:
        effGLCperc[n] = (AArm[n][2] *100) / glcStdBMP

    for n in NonEssentialAAs:
        effGLCvar[n] =  effGLCperc[n] - 100

    effObjValpercNH4["AllAArm"] = (ObjValBMPnh4 *100) / ObjValStdBMP
    effGLCpercNH4["AllAArm"] = (glcBMPnh4 *100) / glcStdBMP
    effObjValvarNH4["AllAArm"] = effObjValpercNH4["AllAArm"] - 100
    effGLCvarNH4["AllAArm"] = effGLCpercNH4["AllAArm"] - 100

    print  "This is the effect of removal of single AA/NH4 upon the objective value for BMP: ", "\n"
    for n in sorted(NonEssentialAAs):
        print n, "", "","","", "%", round(effObjValperc[n],3), "", "","","", "variation", round(effObjValvar[n],3)

    print  "This is the effect of removal of single AA/NH4 upon the GLC consumption for BMP: ", "\n"
    for n in sorted(NonEssentialAAs):
        print n, "", "","","", "%", round(effGLCperc[n],3), "", "","","", "variation", round(effGLCvar[n],3)

    print  "This is the effect of removal of all AA upon the objective value for BMP: ", "\n"
    print "AllAArm", "", "","","", "%", round(effObjValpercNH4["AllAArm"],3), "", "","","", "variation", round(effObjValvarNH4["AllAArm"],3), "\n"

    print  "This is the effect of removal of all AA upon the GLC consumption for BMP: ", "\n"
    print "AllAArm", "", "","","", "%", round(effGLCpercNH4["AllAArm"],3), "", "","","", "variation", round(effGLCvarNH4["AllAArm"],3)

    return effObjValperc, effObjValvar, effGLCperc, effGLCvar, effObjValpercNH4, effGLCpercNH4, effObjValvarNH4, effGLCvarNH4

