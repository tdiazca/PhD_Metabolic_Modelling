import os, sys
import ScrumPy

from ScrumPy.Util import Set
from ScrumPy.Structural import StoMat

from Bioinf import PyoCyc

import Unwanted, Substitutes # These must be in the cwd of the importing module


Corrections = "Corrections.spy"
CorrecMod = ScrumPy.Model(Corrections)
CorrecMod.Hide()

DirectionMap= {
    "IRREVERSIBLE-LEFT-TO-RIGHT" : StoMat.t_Irrev,
    "IRREVERSIBLE-RIGHT-TO-LEFT" : StoMat.t_BackIrrev,
    "PHYSIOL-LEFT-TO-RIGHT"      : StoMat.t_Irrev,
    "PHYSIOL-RIGHT-TO-LEFT"      : StoMat.t_BackIrrev,
    "LEFT-TO-RIGHT"              : StoMat.t_Rever,
    "RIGHT-TO-LEFT"              : StoMat.t_BackIrrev,
    "REVERSIBLE"                 : StoMat.t_Rever
}
DefaultDirec = "IRREVERSIBLE-LEFT-TO-RIGHT"
NoCompID = 'NoComp'

ReacDirKey = PyoCyc.Tags.ReacDir
       

def RemovePipes(sm):
    sm.rnames  = map(lambda x: x.replace("|",""),sm.rnames)
    sm.Externs = map(lambda x: x.replace("|",""),sm.Externs)

def RemoveBadMets(sm):

    BadMatches = filter(lambda met: Unwanted.BadMetMatches(met), sm.rnames)

    for met in Set.Intersect(Unwanted.Mets, sm.rnames)+BadMatches:
                
        for reac in sm.InvolvedWith(met):
            sm.DelReac(reac)
        try:
            sm.DelRow(met)
        except:
            print met, " not found"
            # this is a bug, Set.Intersect should ensure this doesn't happen



def RemoveBadReacs(sm):

    BadMatches = filter(lambda met: Unwanted.BadReacMatches(met), sm.cnames)
    for r in Set.Intersect(Unwanted.Reacs,sm.cnames)+BadMatches:
        try:        sm.DelReac(r)
        except: pass # lazy 


def AddCorrections(sm):

    CorrecMod.Reload()
    # TODO: check for syntax errors in CorrecMod
    csm = CorrecMod.smx

    for c in Set.Intersect(csm.cnames, sm.cnames):
        revprop = csm.RevProps[c]
        sm.Delete(c)
        sm.NewReaction(c, csm.InvolvedWith(c), revprop)

        
def ReverseIsomerases(sm):
    """ make isomerases in sm reversible """

    for c in sm.cnames:
        invw = sm.InvolvedWith(c)
        if len(invw) == 2 and abs(invw.values()[0]) == abs(invw.values()[1]) == 1:
            sm.RevProps[c] = StoMat.t_Rever


def SubMetabolites(sm):
 
    reload(Substitutes)    
    subdic = Substitutes.MetNames

    for met in subdic.keys():
       if met in sm.rnames[:]:
            row = sm[met]
            newmet = subdic[met]
            print "subst met", met, newmet
            if newmet in sm.rnames:
                sm.AddRow(newmet, row)
            else:
                sm.NewRow(row,newmet)
          
            sm.DelRow(met)
     

def FixNads(sm): # TODO - re-write using proper sm.Add/DelReac()
 
    def AsNAD(metname):
        return metname.replace("(P)","")

    def AsNADP(metname):
        return metname.replace("(P)","P")

    def NADPRName(rname):
        return rname+'-(NADP)'

    def NADRName(rname):
        return rname+'-(NAD)'

    reacs = sm.InvolvedWith("NAD(P)H").keys()
    if len(reacs)>0:
        reacsd = {}
        revprops = {}

        for reac in reacs:
            reacsd[reac] = sm.InvolvedWith(reac)
            revprops[reac] = sm.RevProps[reac]
            sm.DelCol(reac)        
        
        sm.DelRow("NAD(P)H")
        sm.DelRow("NAD(P)")

        for reac in reacsd.keys():
            reacnad  = NADRName(reac)
            sm.NewCol(name=reacnad)
            reacnadp = NADPRName(reac)
            sm.NewCol(name=reacnadp)
            
            stod = reacsd[reac]
            for met in stod.keys():
                sm[AsNADP(met),reacnadp] = stod[met]
                sm[AsNAD(met),reacnad] = stod[met]
            
            sm.RevProps[reacnad] = revprops[reac]
            sm.RevProps[reacnadp] = revprops[reac]
            del sm.RevProps[reac]
     

def DelIsoforms(sm):

    for iso in sm.FindIsoforms():
        for reac in iso[1:]:
           sm.DelReac(reac)
    return

def GetQuoteDic(sm):

    rv = {}
    for name in sm.cnames+sm.rnames:
        if name[0] == '"':
            qname = name
            name = name[1:-1]
        else:
            qname = '"' + name + '"'
        rv[name] = qname
        rv[qname] = name
    return rv
    

 
def MakeCompNames(sm,comp,CompDic):

    cnames = CompDic.MapNames(sm.cnames,comp)
    NewOld = zip(cnames, sm.cnames)
    
    for new,old in NewOld:
        sm.RevProps[new] = sm.RevProps[old]
        del sm.RevProps[old]
        
    sm.cnames = cnames
    sm.rnames = CompDic.MapNames(sm.rnames,comp)

  




def MakeCore(m):

    Externs = m.sm.Externs

    orphans  = m.sm.OrphanMets()    
    while len(orphans) > 0:
        print "len(o)", len(orphans)
        reacd = {}
        for met in orphans:
            for reac in m.sm.InvolvedWith(met).keys():
                reacd[reac] =1

        m.DelReactions(reacd.keys())
        orphans  = m.sm.OrphanMets()

def BuildModel(Reactions,CompDic,db):
    reload(Substitutes)
    reload(Unwanted)
    Unwanted.Reacs =  Set.MakeSet(Unwanted.Reacs) # ensure duplicates are removed
    Unwanted.Mets = Set.MakeSet(Unwanted.Mets)
    CorrecMod.Hide()

    smdic = {}
    if len(CompDic)==0:
        CompDic[NoCompID] = NoCompID
        CompDic.DefComp = NoCompID
   
    for comp in CompDic:
        smdic[comp] = StoMat.StoMat(Conv=float)

    for reac in Reactions:
        if not db.Reaction.has_key(reac):
            print >>sys.stderr, "!!", reac, "not found !!"
        else:
            reacrec = db[reac]
            stod = dict(reacrec.CoeffDic)
            if reacrec.has_key(ReacDirKey):
                direc = reacrec[ReacDirKey][0]  
            else:
                direc = DefaultDirec
            direc = DirectionMap[direc]
            
            for comp in CompDic.Reac2Comps(reac):
                smdic[comp].NewReaction(reac, stod,direc)

    for comp in smdic:
        
        out = open(comp+".spy","w")
        print >>out, "\n\n#\n##\n### ", comp, "\n##\n#\n\n"
        sm = smdic[comp]
       
    
        RemovePipes(sm)
        RemoveBadMets(sm)
        
        RemoveBadReacs(sm)
        SubMetabolites(sm)
        
        AddCorrections(sm)  
        ReverseIsomerases(sm) 
        FixNads(sm)  
        DelIsoforms(sm)

        if comp != NoCompID:
            MakeCompNames(sm,comp,CompDic)
            
        qdic = GetQuoteDic(sm)
        sm.ToScrumPy(qdic, out)
        

def RebuildModel(m,Reactions, CompDic,db):

    reload(Substitutes)
    reload(Unwanted)

    m.Hide()

    BuildModel(Reactions, CompDic,db)
    m2 =  ScrumPy.Model(m.md.GetRootFilePath())
    del m
    return m2

    
