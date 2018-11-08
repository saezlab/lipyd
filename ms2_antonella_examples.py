

# BMP/PG
# BMP(36:2)
pbmp         = l.get_scan('BPIFB2',  'pos', 792.57)
# BMP(36:2)
nbmp         = l.get_scan('BPIFB2',  'neg', 773.5258)
# PG/BMP(38:3)
npg          = l.get_scan('GM2A',    'neg', 799.54)

# Ceramides: HexCer, HexCer-OH, Cer, Cer1P
# 
nhexcer      = l.get_scan('GLTP',    'neg', 744.5627)
nhexceroh    = l.get_scan('GLTP',    'neg', 760.557)
ndihexcer    = l.get_scan('GLTP',    'neg', 1004.689)
ncerp        = l.get_scan('GLTPD1',  'neg', 616.47)
# Cer(34:1)
ncer         = l.get_scan('SEC14L1', 'neg', 582.509)
ncer2        = l.get_scan('STARD11', 'neg', 582.51)

phexcer      = l.get_scan('GLTP',    'pos', 810.68)
phexceroh    = l.get_scan('GLTP',    'pos', 826.67)
pdihexceroh  = l.get_scan('GLTP',    'pos', 988.73)
pcerp        = l.get_scan('GLTPD1',  'pos', 728.59)
pcerp2       = l.get_scan('GLTPD1',  'pos', 590.45)
pcerp3       = l.get_scan('GLTPD1',  'pos', 702.58)
pcerp4       = l.get_scan('GLTPD1',  'pos', 618.430)
pcerp5       = l.get_scan('GLTPD1',  'pos', 616.415)
pcerp6       = l.get_scan('GLTPD1',  'pos', 640.409)

pcer         = l.get_scan('SEC14L1', 'pos', 538.52)
pcer2        = l.get_scan('STARD11', 'pos', 538.526)

# SM/PC
nsm          = l.get_scan('GLTPD1',  'neg', 745.55)
npc          = l.get_scan('BPI',     'neg', 804.57)

psm          = l.get_scan('GLTPD1',  'pos', 703.57)
ppc          = l.get_scan('BPI',     'pos', 760.58)

# PI
npi          = l.get_scan('BPI',     'neg', 861.55)

ppi          = l.get_scan('SEC14L2', 'pos', 906.607)

# PS
pps          = l.get_scan('ORP9',    'pos', 790.556)

nps          = l.get_scan('ORP9',    'neg', 788.544)
nps2         = l.get_scan('BPI',     'neg', 760.51)

# PE
npe          = l.get_scan('GM2A',    'neg', 714.507)

ppe          = l.get_scan('BPI',     'pos', 718.536)

# DAG/TAG
## DAG(32:1) => (16:0/16:1)
pdag         = l.get_scan('SEC14L2', 'pos', 584.52)
## TAG(48:3) => (16:1/16:1/16:1) or (16:1/14:1/18:1)
ptag         = l.get_scan('STARD11', 'pos', 818.724)

# Vitamin A
pva          = l.get_scan('RBP1',    'pos', 269.226)
pva2         = l.get_scan('RBP4',    'pos', 269.226)
