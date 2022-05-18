import LP  # The ScrumPy linear programming module

import SepiBiomass
reload(SepiBiomass) # description of biomass (BM) composition

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

def BPtx(m):

    """pre: m = the model
       post: returns a list of by-product transporters"""
    
    rv = []

    for r in filter(lambda s: "_bp_tx" in s, m.smx.cnames):
        rv.append(r)

    return rv

def AAmediatx(m):

    """pre: m = the model
       post: returns a list of media amino acid transporters"""
    
    rv = []

    for r in filter(lambda s: "_AA_mm_tx" in s, m.smx.cnames):
        rv.append(r)

    return rv

def AAbiomasstx(m):

    """pre: m = the model
       post: returns a list of biomass amino acid transporters"""
    
    rv = []

    for r in filter(lambda s: "_AA_bm_tx" in s, m.smx.cnames):
        rv.append(r)

    return rv

def BuildLP(m):
    
    """ build lp from model m minimising total flux """
    
    lp = LP.lp(m) # lp = m.GetLP()  
    lp.SetObjective(m.sm.cnames)   # minimise total flux

    return lp

##def BuildMinGLCLP(m):
##    
##    """ lp object that minimises GLC import """
##    
##    lp = LP.lp(m) # lp = m.GetLP()  
##    lp.SetObjective("GLC_mm_tx")   # minimise GLC import = optimise usage
##
##    return lp

def BuildMinObjective(m, Objective='GLC_mm_tx'):
    
    """ lp object that minimises flux through reaction defined as Objective """
    
    lp = LP.lp(m) # lp = m.GetLP()  
    lp.SetObjective(Objective)   # minimise GLC import = optimise usage

    return lp

def BuildBiomassLP(m, fd=None):
    
    """ lp object constrained to biomass fluxes defined by BM composition """

    if fd==None:
        fd = SepiBiomass.FluxDic(m)
        # generates a dic of fluxes according to BM description
            # for producing 1 arbitrary unit of biomass (default constraint val of 1)

    lp = BuildLP(m)
    lp.SetFixedFlux(fd) # sets output flux to biomass comp parameters
                        # deffault value is 1! (flux constrain)
    return lp

def BuildBiomassATPmaintLP(m, AMaint=45.0, fd=None):
    # AMaint = GA+NGA maintenance cost = 40 + 5
        #(assuming growth rate 1gh-1)
    
    """ lp object constrained to biomass fluxes defined by BM composition
        and ATP maintance cost """

    if fd==None:
        fd = SepiBiomass.FluxDic(m)
        #fd["ATPase"] = AMaint
        # generates a dic of fluxes according to BM description
            # for 1 arbitrary unit of BM (default constraint val is 1)

    lp = BuildLP(m)
    lp.SetFixedFlux(fd) # sets output flux to biomass comp parameters ; deffault value is 1!
    lp.SetFluxBounds({'ATPase':(AMaint, None)})
    return lp

def BuildPropBiomassATPmaintLP(m, PropBF=0.0, AMaint=45.0, fd=None): ## just add prop of BF
    
    """ lp object constrained to prop of BF in BM and biomass fluxes
        defined by BM composition and ATP maintance cost """

    if fd==None:
        bm = SepiBiomass.All.Copy()
        bm.SetAmmount(PropBF, "Biofilm")
        bm.SetAmmount(1.0-PropBF, "Cells")
        fd = SepiBiomass.FluxDic(m, bm)
        #fd["ATPase"] = AMaint
        
    lp = BuildLP(m)
    lp.SetFixedFlux(fd)
    lp.SetFluxBounds({'ATPase':(AMaint, None)})
    return lp

##def BuildBiomassLP(m, PropBiofilm = 0.0): ## Doesnt account for ATP demand!!
##    """ pre: 0.0 <=PropBiofilm <= 1.0  """
##
##    lp = BuildLP(m)
##    fd = SepiBiomass.GetFD(m, PropBiofilm)
##    lp.SetFixedFlux(fd)
##    
##    return lp

def BuildClosedLP(m):

    """ lp object contraining flux through tx reacs to zero """

    lp = BuildLP(m)
    
    for tx in filter(lambda s: "_tx" in s, m.sm.cnames): # all tx in model
        lp.SetFixedFlux({tx:0.0}) # set flux to 0.0

    return lp
