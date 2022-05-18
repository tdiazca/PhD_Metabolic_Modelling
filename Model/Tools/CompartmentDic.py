

EnQuote = lambda s: '"'+s+'"'
IsQuote = lambda s: s[0] == s[-1] == '"'
DeQuote = lambda s: s[1:-1] if IsQuote(s) else s

MakeSfx = lambda s: "_"+s[:4].title()
## not a class member so that user can redefine on the fly


class CompartmentDic(dict):

    def __init__(self, CompDic={}, DefaultCompartment="Cytosol"):

        self.update(CompDic)
        self.DefComp = DefaultCompartment     
        self.SfxMap = dict(zip(self.keys(), map(MakeSfx, self.keys())))


    def Reac2Comps(self, Reac):
        
        rv = []
        for comp in self:
            if Reac in self[comp]:
                rv.append(comp)
        if rv == []:
            rv.append(self.DefComp)

        return rv


    def MapNames(self, names,comp):

        sfx =  self.SfxMap[comp] if self.SfxMap.has_key(comp) else MakeSfx(comp)
        return map(lambda name: name+sfx, names)


    def AllReacs(self):
        """ list of all reacs in all compartments """

        rvd = {}
        for comp in self.values():
            for reac in comp:
                rvd[reac] = 1
                
        return rvd.keys()
                
    

    

    

    
