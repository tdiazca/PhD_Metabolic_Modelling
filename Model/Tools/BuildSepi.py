#
##
### Build the S. epi RP62A model
##
#

OrgName = 'Staph_epi_RP62A_20.1' # name of the database passed to PyoCyc.Organism
CompartModName = None
DBUpdateName = "ExtraCompounds.dat"

UpdateFrom = { # This dic is not being used, but kept for record
    "ACETOIN" : "CPD-255",
    "CPD-10490" : "CPD-10490", # 1 H too many in BioCyc RP62A
    "ADP-D-GLUCOSE" : "ADP-D-GLUCOSE",
    "ARABINOSE" : "ARABINOSE",
    "ISOVALERYL-COA" : "ISOVALERYL-COA",
    "3-METHYL-CROTONYL-COA" : "3-METHYL-CROTONYL-COA",
    "RIBOSE" : "CPD0-1110",
    "RIBOSE-5P" : "CPD-15895",
    "CPD-12297" : "CPD-12297",
    "D-Glucopyranuronate" : "CPD-12521",
    "DEOXY-RIBOSE-5P" : "DEOXY-RIBOSE-5P",
    "D-GLUCOSAMINE-6-P" : "CPD-13469",
    "N-ACETYL-D-GLUCOSAMINE-6-P" : "N-ACETYL-D-GLUCOSAMINE-6-P",
    #"D-SEDOHEPTULOSE-7-P" : "D-SEDOHEPTULOSE-7-P", unknown in db version 20.5
    "MENAQUINONE" : "CPD-9728",
    "O-SUCCINYL-L-HOMOSERINE" : "O-SUCCINYL-L-HOMOSERINE",
    "CPD-12297" : "CPD-12297",
    "CPD-12242" : "CPD-12242",
    "CPD-12243" : "CPD-12243",
    "CPD-1111"  : "CPD-1111",
    "CPD-19301" : "CPD-19301",
    "CPD-19302" : "CPD-19302",
    "CPD-19303" : "CPD-19303",
    "CPD-19305" : "CPD-19305",
    "CPD-19317" : "CPD-19317",
    "CPD-372"   : "CPD-372",
    "BUTANEDIOL" : "BUTANEDIOL", 
    "CPD-10490" : "CPD-10490",
    "CPD-12521" : "CPD-12521",
    "PALMITATE" : "PALMITATE",
    "S-ADENOSYLMETHIONINE" : "S-ADENOSYLMETHIONINE",
    "ISOBUTYRYL-COA" : "ISOBUTYRYL-COA",
    "CPD-12557" : "CPD-12557",
    "D-GLUCOSAMINE-6-P" : "D-GLUCOSAMINE-6-P",
    "Cytochromes-C-Oxidized" : "Cytochromes-C-Oxidized",
    "Cytochromes-C-Reduced" : "Cytochromes-C-Reduced",
    "ETF-Oxidized" : "ETF-Oxidized",
    "ETF-Reduced" : "ETF-Reduced",
    "ETR-Quinols" : "ETR-Quinols",
    "ETR-Quinones" : "ETR-Quinones",
    "Ox-Thioredoxin" : "Ox-Thioredoxin",
    "Red-Thioredoxin" : "Red-Thioredoxin",
    "Oxidized-NrdH-Proteins" : "Oxidized-NrdH-Proteins",
    "Reduced-NrdH-Proteins" : "Reduced-NrdH-Proteins",
    "Oxidized-ferredoxins" : "Oxidized-ferredoxins",
    "Reduced-ferredoxins" : "Reduced-ferredoxins",
    "Oxidized-flavodoxins" : "Oxidized-flavodoxins",
    "Reduced-flavodoxins" : "Reduced-flavodoxins",
    "Ox-Glutaredoxins" : "Ox-Glutaredoxins",
    "Red-Glutaredoxins" : "Red-Glutaredoxins",
    "GLY-tRNAs" : "GLY-tRNAs",
    "Charged-GLY-tRNAs" : "Charged-GLY-tRNAs",
    "ASN-tRNAs" : "ASN-tRNAs",
    "Charged-ASN-tRNAs" : "Charged-ASN-tRNAs",
    "L-aspartyl-tRNAAsn" : "L-aspartyl-tRNAAsn",
    "5-10-METHENYL-THF-GLU-N" : "5-10-METHENYL-THF-GLU-N",
    "METHYLENE-THF-GLU-N" : "METHYLENE-THF-GLU-N",
    "5-METHYL-THF-GLU-N" : "5-METHYL-THF-GLU-N",
    "THF-GLU-N" : "THF-GLU-N",
    "D-alanine-carrier-protein" : "D-alanine-carrier-protein",
    "D-Ala-DltC" : "D-Ala-DltC",
    "1-Phosphatidyl-2-O-D-Ala-Glycerol" : "1-Phosphatidyl-2-O-D-Ala-Glycerol",
    "Type-I-LTA" : "Type-I-LTA",
    "CDPdag" : "CDPdag",
    "PIA1" : "PIA1",
    "PIA2" : "PIA2",
    "glycogen" : "glycogen",
    "PALMITATE" : "PALMITATE",
    "L-PHOSPHATIDATE" : "L-PHOSPHATIDATE",
    "L-1-PHOSPHATIDYL-GLYCEROL" : "L-1-PHOSPHATIDYL-GLYCEROL",
    "L-1-PHOSPHATIDYL-GLYCEROL-P" : "L-1-PHOSPHATIDYL-GLYCEROL-P",
    "CDPDIACYLGLYCEROL" : "CDPDIACYLGLYCEROL",
    "CARDIOLIPIN" : "CARDIOLIPIN",
    "DIACYLGLYCEROL" : "DIACYLGLYCEROL",
    "D-Glucosyl-12-diacyl-glycerols" : "D-Glucosyl-12-diacyl-glycerols",
    "diacyl-3-O-glucl-1-6-gluc-sn-glycerol" : "diacyl-3-O-glucl-1-6-gluc-sn-glycerol",
    "Gro-P-Glc2-DAG" : "Gro-P-Glc2-DAG",
    "Gro-P-n-Gro-P-Glc2-DAG" : "Gro-P-n-Gro-P-Glc2-DAG",
    "GlcNAc-Gro-P-n-Gro-P-Glc2-DAG" : "GlcNAc-Gro-P-n-Gro-P-Glc2-DAG",
    "x_PROTON" : "x_PROTON",
    "CPD-19301" : "CPD-19301",
    "CPD-19311" : "CPD-19311",
    "CPD-19311_peri" : "CPD-19311_peri",
    "O-SUCCINYL-L-HOMOSERINE" : "O-SUCCINYL-L-HOMOSERINE",
    "x_AWork" : "x_AWork",
    "CPD-12297" : "CPD-12297",
    "O-SUCCINYL-L-HOMOSERINE" : "O-SUCCINYL-L-HOMOSERINE",
    "ADP-D-GLUCOSE" : "ADP-D-GLUCOSE",
    "CPD-372" : "CPD-372",
    "XYLOSE" : "XYLOSE",
    "CPD-19317" : "CPD-19317",
    "CPD-19311" : "CPD-19311",
    "CPD-1111" : "CPD-1111",
    "CPD-19793" : "CPD-19793",
    "Rbo-P-Teichoic-aureus-peptidoglycan" : "Rbo-P-Teichoic-aureus-peptidoglycan"
    }


import UsrBuildOrg

def Init():
    UsrBuildOrg.Init(OrgName,DBUpdateName)

Init()

from UsrBuildOrg import BuildModel, RebuildModel, ShowCorrec, HideCorrec, orgdb

