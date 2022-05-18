
from Bioinf import PyoCyc
import BuildOrg, CompartmentDic
# CompartmentDic contains the class that will be populated by CompDic

    
orgdb = None

class NulComp: # default CompDic if CompartModName is None
        CompDic  = {}
CompDic = NulComp() # overwritten by Init() if CompartModName is not None

def ShowCorrec():
    BuildOrg.CorrecMod.Show()

def HideCorrec():
    BuildOrg.CorrecMod.Hide()

    
def Init(OrgDB, ExtraMetsDB=None,CompartModName=None ):
        """
        OrgDB : name of the organsim db
        ExtraMetsDB : name of a database with additional compound data
        CompartModName : Name of a module assigning db IDs to compartments
        """
       
        global orgdb
        global CompDic
        
        if orgdb == None:
                orgdb = PyoCyc.Organism(data=OrgDB)

        if ExtraMetsDB !=None:              
                updatedb = PyoCyc.Compound.DB(".", ExtraMetsDB)
                for k in updatedb.keys():
                        orgdb.Compound[k] = updatedb[k]
                       
  
        if CompartModName != None:
                exec "import %s as cd" % CompartModName
                CompDic = cd

def ReloadCD():
        try:
                reload(CompDic)
        except:
                pass
    
    
def BuildModel():
    ReloadCD()
    reacs = orgdb.Reaction.keys()
    cd = CompartmentDic.CompartmentDic(CompDic.CompDic)
    BuildOrg.BuildModel(reacs, cd, orgdb)

def RebuildModel(m, Hide=True):
        ReloadCD()
        reacs = orgdb.Reaction.keys()

        cd = CompartmentDic.CompartmentDic(CompDic.CompDic)
        m = BuildOrg.RebuildModel(m,reacs,cd, orgdb)
        if Hide:
            m.Hide()
        return m

