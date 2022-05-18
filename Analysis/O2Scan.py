######
########
##########   This module uses LPA for studying reactions that responde to variation
###########     on O2 concentration for production of biomass or ATP
##########
########
######    
####

##    NOTE: remember, "GLN_AA_mm_tx", "NO3_mm_tx" to be blocked 
##              if medium is MMc media but available in synovial fluid

import ScrumPy
from Util import Set
from Data import DataSets
import BuildLP
#from Util import Sci
#from Util import Set
from Bioinf import Biomass
import SepiBiomass
    # need to specify proportion of biofilm comps later

### Reacs that create multiple optima =
    
Ignore = ['GMP-SYN-NH3-RXN','NAD-SYNTH-NH3-RXN','CYTIDINEKIN-RXN',
          'L-GLN-FRUCT-6-P-AMINOTRANS-RXN']

    # Reacs causing multiple optima =

    # 'GMP-SYN-NH3-RXN','GMP-SYN-GLUT-RXN', rm = 'GMP-SYN-NH3-RXN'
    # 'NAD-SYNTH-NH3-RXN','NAD-SYNTH-GLN-RXN', rm = 'NAD-SYNTH-NH3-RXN'
    # 'CYTIDINEKIN-RXN','CYTIKIN-RXN' rm = 'CYTIDINEKIN-RXN'
    # 'GLUCOSAMINE-6-P-DEAMIN-RXN','L-GLN-FRUCT-6-P-AMINOTRANS-RXN',
            # rm = 'L-GLN-FRUCT-6-P-AMINOTRANS-RXN'

    # 'D-ALANINE-AMINOTRANSFERASE-RXN','ALANINE-AMINOTRANSFERASE-RXN',
            # not always, leave in for now
    # 'PYRUVFORMLY-RXN','PYRUVDEH-RXN', # OK, no mul optima!


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
    bm = Biomass.Composition(
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
    #lp.SetFixedFlux(fd)
    for r in fd.keys():
        if 'bm_tx' in r:
            lp.SetFixedFlux(fd)
    lp.SetFluxBounds({"ATPase": (AMaint, None)}) # If we dont want ATPase upper flux to be fixed

    return lp

    ##lp.SetFluxBounds({"ATPase": (AMaint, None)})
    ### If we dont want ATPase upper fluxb bound to be fixed

def GetLPfreeBM(m, PropBF=0.0, AMaint=45.0):
    """pre: GetFluxDic(m, PropBF=0.0, AMaint=45.0)     
       post: lp object constrained to biomass fluxes defined by BM composition
       and ATP maintenace demand, with upper flux bounds for BM transporters relaxed"""
    orig_rp = dict(m.sm.RevProps)
    for r in filter(lambda s: "bm_tx" in s, m.sm.cnames):
        m.sm.RevProps[r] = m.smx.RevProps[r] = "<>" # Otherise, tx could take random values (bug in ScrumPy!)

    lp = BuildLP.BuildLP(m) # lp = LP.lp(m) ; lp = m.GetLP() # lp.SetObjective(m.sm.cnames)
    
    fd = GetFluxDic(m, PropBF, AMaint)
    
    for r in fd.keys():
        if 'bm_tx' in r:    # otherwise allows ATPase to go to zero flux! (since also in fd)
            lp.SetFluxBounds({r:(None,fd[r])}) # so bm_tx products can be exported at higher prop as in BM
    m.sm.RevProps = m.smx.RevProps = orig_rp
    
    lp.SetFluxBounds({"ATPase": (AMaint, None)}) # If we dont want ATPase upper flux to be fixed
                    
    return lp

def FRange(Start, Stop, nPoints):
    """pre: Start > Stop
       post: returns a list of n points increasing value at a fixed rate
               between Start and Stop"""
    Start = float(Start)
    incr  = (Stop-Start)/(nPoints -1)
    return map(lambda x: Start+incr*x, range(nPoints))
        # for each point in the range nPoints, calculate the corresponding increment
            # lambda repeats a fx over an over for each item on a list

class ScanResult(DataSets.DataSet):

    def Changers(self, thresh=1e-3):
        """pre: self= O2Scan(...)
           post: returns list of reactions in self which flux changes by more than thresh
                   caution, includes O2Lim and ObjVal"""
        changers=[]
        changers=filter(lambda c: self.Range(c) > thresh, self.cnames)
        return changers

    def ChangingReacs(self,m,thresh=1e-3):
        """pre: self= O2Scan(...)
           post: returns list of reactions in self which flux changes by more than thresh"""
        changing=[]
        changers=ScanResult.Changers(self, thresh)
        changing=Set.Intersect(changers, m.smx.cnames)
        print "These are the reactions which flux responds to the variation in the O2 concentration: ",len(changing), "\n", changing, "\n"
        return changing

    def ChangingTx(self,m,thresh=1e-3):
        """pre: self= O2Scan(...)
           post: returns list of transporting reacs
           in self which flux changes by more than thresh"""
        changingTx=[]
        changing=ScanResult.ChangingReacs(self,m,thresh)
        changingTx=Set.Intersect(filter(lambda s: '_tx' in s, m.smx.cnames), changing)
        print "These are the transporters which flux responds to the variation in the O2 concentration: ",len(changingTx), "\n", changingTx, "\n"
        return changingTx

    def Switchers(self,thresh=1e-3):
        """pre: self= O2Scan(...)
           post: returns list of reactions that have at least one value < than thresh"""
        switchers = []
        for c in self.cnames:
            if min(map(abs, self.GetCol(c))) < thresh: # if at any point val in column is < than
                switchers.append(c)
        return switchers

    def ReacsInDataset(self,m):
        """pre: self= O2Scan(...)
           post: returns list of reactions active at any point of the scan"""
        AllActiveReacs=[]
        AllActiveReacs=Set.Intersect(self.cnames,m.sm.cnames)
        print "These are the reactions active at any point of the scan: ",len(AllActiveReacs), "\n", AllActiveReacs, "\n"
        return AllActiveReacs

    def ActiveInRow(self,r):
        """pre: self= O2Scan(...)
           post: returns a dic with reacs active in sol in dataset number 'r' as keys
                   and their fluxes as values"""
        d = self.RowDict(r) # r = index for rows in dataset
        return dict(filter(lambda x: d.items(), x[1]!=0)) 

    def CalcAABal(self): # net flux in media tx for AAs
        """pre: m = self = O2Scan(m),
           post: add new col to self (dataset) with net flux for AA txs"""
        for mtx in filter(lambda s: "AA_mm" in s, self.cnames): # AA media transp
            btx = mtx.replace("_mm_", "_bm_")
            bal = map(lambda m,b: m-b, self.GetCol(mtx), self.GetCol(btx))
            aa = mtx.split("_")[0]
            self.NewCol(bal, aa+"_BAL")

    def PlotChangingReacs(self,m,thresh=1e-3):
        """pre: self= O2Scan(...)
           post: returns a plot with O2Lim values in the X exe and flux through
           reactions changing flux in the Y exe"""
        changing=ScanResult.ChangingReacs(self,m,thresh)
        return self.AddToPlot(changing)

    def PlotChangingTxs(self,m,thresh=1e-3):
        """pre: self= O2Scan(...)
           post: returns a plot with O2Lim values in the X exe and flux through
           transporting reactions changing flux in the Y exe"""
        changingTx=ScanResult.ChangingTx(self,m,thresh)
        return self.AddToPlot(changingTx)

    def MakeSubModel(self,m,thresh=1e-3,File='SubModel.spy'):
        """pre: res = O2Scan(m), thresh>0;
           post: sub-model containing reacs extracted from Changers in dataset"""
        reacs = Set.Intersect(self.Changers(thresh),m.sm.cnames) # removes 'ObjVal' and 'O2Lim'
        sm = m.smx.SubSys(reacs)
        sm.ToScrumPy(m.md.QuoteMap,File,Externs=m.Externals())
                # there would probably be unconserved mets (to be considered as externals here)

    def MakeSubModelWholeDataSet(self,m,File='SubModel.spy'):
        """pre: self = O2Scan(m),
           post: sub-model containing all reacs extracted from dataset"""
        reacs = Set.Intersect(self.cnames,m.sm.cnames) # removes 'ObjVal' and 'O2Lim'
        sm = m.smx.SubSys(reacs)
        sm.ToScrumPy(m.md.QuoteMap,File,Externs=m.Externals())
        # writes a File describing a model containing reacs

    def ElementaryModesSubmodel(self,m,File='SubModel.spy'):
        """pre: self= O2Scan(...)
                MakeSubModelWholeDataSet(...) = generates is a file containing all reactions active at any point of the scan
           post: returns the net stoichiometries of the elementary modes of the submodel"""
        ScanResult.MakeSubModelWholeDataSet(self,m,File='SubModel.spy')
        mm=ScrumPy.Model(File)
        #k=m.sm.NullSpace()
        #k.IntegiseC()          # show entries in k matrix as integers
        #print k

        elmo = mm.ElModes()
        elmo.Integise()         # show elmo coeficients as integers

        print 'There are n elementary modes in the subnetwork: ',len(elmo),"\n"
        print 'These are the net stoichiometries of the elementary modes of the submodel: ',"\n","\n",elmo.Stos()   
        return elmo


def O2Scan(m, Lo=0.0, Hi=12, nPoints = 100, PropBF=0.0,
        AMaint=45.0, BlockAA=False, BlockOther=[], FreeBM=None):
    """
	pre: m = the model
	     Lo = lower limit for upper flux bound of the O2 importer (O2_mm_tx); Lo <=0.0
	     Hi = upper limit for upper flux bound of the O2 importer; Hi>Lo
	     nPoints = number of points to sample ; >1
	     PropBF = proportion of biofilm in biomass; 0.0 <= PropBF <= 1.0
             AMaint = ATP demand (NGAM,GAM) growth rate 1gDW/h = 45 mmol/gDW;
                      0.0 <= AMaint <= 45.0
	     BlockAA = blocks flux through AA importers (AA_mm_tx)
	     BlockOther = blocks flux through other reacs
             FreeBM = if != None, upper flux bound of BM_tx (biomass txs) is unbound
	post: returns a dataset (ds) containing the scan results. ds.lp will be the lp object
    """

    rv = ScanResult()

    if FreeBM != None:
        lp = GetLPfreeBM(m, PropBF, AMaint)
    else:
        lp = GetLP(m, PropBF, AMaint)

    AAtx = []
    if BlockAA:
        AAtx = filter(lambda s: "AA_mm_tx" in s, m.sm.cnames)

    BlockReacs = Ignore+AAtx+BlockOther
    lp.SetFixedFlux(dict(zip(BlockReacs, [0.0]*len(BlockReacs))))

    O2_tx = 'O2_mm_tx'
    SetO2 = lambda O2Lim: lp.SetFluxBounds({O2_tx:(0.0,O2Lim)})   
    
    for O2Lim in FRange(Lo, Hi, nPoints):
        SetO2(O2Lim)
        lp.Solve(False)
        print ".", 
        if lp.IsStatusOptimal():
            sol = lp.GetPrimSol()
            sol['ObjVal'] = lp.GetObjVal()
            sol['O2Lim'] = O2Lim
            
            if FreeBM != None:
                fd=GetFluxDic(m, PropBF, AMaint) 
                for r in filter(lambda s: "bm_tx" in s, sol.keys()):
                    for r in Set.Intersect(sol.keys(),fd.keys()):
                        # ATPase and (PIA1_bf_bm_tx, PIA2_bf_bm_tx with flux 0 if PropBF = 0.0)
                        if sol[r] != fd[r]:
                            comp = r.split("_")[0]
                            NetExp = sol[r] - fd[r] 
                            sol[comp+"_NetExport"] = NetExp
                rv.UpdateFromDic(sol)
            else:
                for r in filter(lambda s: "AA_bm_tx" in s, sol.keys()):
                    aa = r.split("_")[0]
                    NetExp = sol[r]
                    sol[aa+"_Export"] = NetExp
                rv.UpdateFromDic(sol)

    rv.SetPlotX('O2Lim')

    if not (PropBF == 1.0): #  if not (PropBF == 1.0 or ATP only):
        rv.CalcAABal()

    rv.lp = lp
    return rv
