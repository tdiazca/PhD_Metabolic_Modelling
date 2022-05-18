######
########
##########   This module uses LPA for studying reactions that responde to variation
###########     on O2 concentration for production of ATP
##########      (either GAM and NGAM costs = 45 mmol/gDW/h  or 1 mol ATP)
########
######    
####

##    NOTE: remember, "GLN_AA_mm_tx", "NO3_mm_tx" to be blocked 
##              if medium is MMc media but available in synovial fluid
                  
##    NOTE: select adequate NTxs list
##              (uncomment list of N containing compounds depending on medium)

import ScrumPy
from Util import Set
from Data import DataSets
import BuildLP
#from Util import Sci

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

def FRange(Start, Stop, nPoints):
    """pre: Start > Stop
       post: returns a list of n points increasing value at a fixed rate
               between Start and Stop"""
    Start = float(Start)
    incr  = (Stop-Start)/(nPoints -1)
    return map(lambda x: Start+incr*x, range(nPoints))
        # for each point in the range nPoints, calculate the corresponding increment
            # lambda repeats a fx over an over for each item on a list

class ScanResult(DataSets.DataSet): # creates a class which has a dataset stucture

    def Changers(self, thresh=1e-3):
        """pre: self= O2ScanATPonly(...)
           post: returns list of reactions in self which flux changes by more than thresh
                   caution, includes O2Lim and ObjVal"""
        changers=[]
        changers=filter(lambda c: self.Range(c) > thresh, self.cnames)
        return changers

    def ChangingReacs(self,m,thresh=1e-3):
        """pre: self= O2ScanATPonly(...)
           post: returns list of reactions in self which flux changes by more than thresh"""
        changing=[]
        changers=ScanResult.Changers(self, thresh)
        changing=Set.Intersect(changers, m.smx.cnames)
        print "These are the reactions which flux responds to the variation in the O2 concentration: ",len(changing), "\n", changing, "\n"
        return changing

    def ChangingTx(self,m,thresh=1e-3):
        """pre: self= O2ScanATPonly(...)
           post: returns list of transporting reacs
           in self which flux changes by more than thresh"""
        changingTx=[]
        changing=ScanResult.ChangingReacs(self,m,thresh)
        changingTx=Set.Intersect(filter(lambda s: '_tx' in s, m.smx.cnames), changing)
        print "These are the transporters which flux responds to the variation in the O2 concentration: ",len(changingTx), "\n", changingTx, "\n"
        return changingTx

    def ReacsInDataset(self,m):
        """pre: self= O2ScanATPonly(...)
           post: returns list of reactions active at any point of the scan"""
        AllActiveReacs=[]
        AllActiveReacs=Set.Intersect(self.cnames,m.sm.cnames)
        print "These are the reactions active at any point of the scan: ",len(AllActiveReacs), "\n", AllActiveReacs, "\n"
        return AllActiveReacs 

    def PlotChangingReacs(self,m,thresh=1e-3):
        """pre: self= O2ScanATPonly(...)
           post: returns a plot with O2Lim values in the X exe and flux through
           reactions changing flux in the Y exe"""
        changing=ScanResult.ChangingReacs(self,m,thresh)
        return self.AddToPlot(changing)

    def PlotChangingTxs(self,m,thresh=1e-3):
        """pre: self= O2ScanATPonly(...)
           post: returns a plot with O2Lim values in the X exe and flux through
           transporting reactions changing flux in the Y exe"""
        changingTx=ScanResult.ChangingTx(self,m,thresh)
        return self.AddToPlot(changingTx)

    def MakeSubModelWholeDataSet(self,m,File='SubModel.spy'):
        """pre: self = O2ScanATPonly(m),
           post: sub-model containing all reacs extracted from dataset"""
        reacs = Set.Intersect(self.cnames,m.sm.cnames) # removes 'ObjVal' and 'O2Lim'
        sm = m.smx.SubSys(reacs)
        sm.ToScrumPy(m.md.QuoteMap,File,Externs=m.Externals())
        # writes a File describing a model containing reacs

    def ElementaryModesSubmodel(self,m,File='SubModel.spy'):
        """pre: self= O2ScanATPonly(...)
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

        #print elmo.ReacsOf('ElMo_0')   

def O2ScanATPonly(m, Lo=0.0, Hi=12, nPoints = 100, AMaint=45.0, BlockAA=False,
                  BlockOther=[]):
    """
	pre: m = the model
	     Lo = lower limit for upper flux bound of the O2 importer (O2_mm_tx); Lo <=0.0
	     Hi = upper limit for upper flux bound of the O2 importer; Hi>Lo
	     nPoints = number of points to sample ; >1
             AMaint = ATP demand (NGAM,GAM) growth rate 1gDW/h = 45 mmol/gDW;
                      0.0 <= AMaint <= 45.0
	     BlockAA = blocks flux through AA importers (AA_mm_tx)
	     BlockOther = blocks flux through other reacs
	post: returns a dataset (ds) containing the scan results. ds.lp will be the lp object
    """

    rv = ScanResult() # creates object of the class

    lp = BuildLP.BuildLP(m) # lp = LP.lp(m) ; lp = m.GetLP() # lp.SetObjective(m.sm.cnames)

    lp.SetFixedFlux({"ATPase": AMaint})

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
            for r in filter(lambda s: "AA_bm_tx" in s, sol.keys()):
                aa = r.split("_")[0]
                NetExp = sol[r]
                sol[aa+"_Export"] = NetExp
            rv.UpdateFromDic(sol)

    rv.SetPlotX('O2Lim')

##    if not (PropBF == 1.0 or ATPOnly):
#    if not (PropBF == 1.0):
#    rv.CalcAABal()

    rv.lp = lp # - make it available to the user
    return rv
