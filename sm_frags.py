
import os
from lipyd import ms2

sm_frags = (
    'PC/SM [N+3xCH3] (60.0808)',
    'PC/SM [Ch] (86.096)',
    'PC/SM [Ch+H2O] (104.107)',
    'PC/SM [P+Et] (124.9998)',
    'PC/SM [Ch-Et] (58.0651)',
    'NL PC/SM [P+Ch] (NL 183.066)',
    'NL SM [P+Ch] (NL 201.0766)',
    'NL SM [N+3xCH3] (77.0841)',
    'NL [H2O] (NL 18.0106)',
)

def check_sm_frags(scan):
    
    for frag in sm_frags:
        
        print(frag, scan.has_fragment(frag))


mgfdir = '/home/denes/archive/ltp/mgf_invitro'
mgfname = '161203_MLH_EMV_GLTPD1_pos_%s.mgf'
fractions = ('E07', 'E08', 'E09', 'E10')

mgfs = dict((fr, os.path.join(mgfdir, mgfname % fr)) for fr in fractions)

mzs = (
    813.6844,
    787.6689,
    705.5904,
    815.7003,
    717.5903,
    729.5914,
    759.6378,
    799.6701,
    787.6691,
)

rts = (
    (11.11, 12.03),
    (11.10, 11.88),
    (7.29, 7.95),
    (11.64, 13.49),
    (7.17, 8.07),
    (7.04, 7.86),
    (9.03, 10.40),
    (10.42, 11.54),
    (10.22, 12.10),
)


hdr = [
    'mz', 'rt', 'scan_id', 'deltaRT', 'frac', 'MS1', '[Sph-2xH2O+H]+'
] + list(sm_frags)

result = []

for mz, rt in zip(mzs, rts):
    
    fe = ms2.MS2Feature(mz, ionmode = 'pos', mgfs = mgfs, rt = rt)
    fe.build_scans()
    
    for sc, drt in zip(fe.scans, fe.deltart):
        
        ids = sc.identify()
        
        for ms1id, ms2ids in ids.items():
            
            if ms2ids and ms2ids[0].hg.main == 'SM':
                
                this_line = [
                    '%.07f' % mz,
                    '%.02f - %.02f' % rt,
                    '%u' % sc.scan_id,
                    '%.02f' % drt,
                    sc.sample,
                    ms1id,
                ]
                this_line.append(
                    sc.has_chain_fragment_type(
                        frag_type = 'Sph-2xH2O+H'
                    ).__str__()
                )
                
                for frag in sm_frags:
                    
                    this_line.append(sc.has_fragment(frag).__str__())
                
                result.append(this_line)

with open('GLTPD1_pos_invitro_SM.csv', 'w') as fp:
    
    _ = fp.write('\t'.join(hdr))
    _ = fp.write('\n')
    
    for line in result:
        
        _ = fp.write('\t'.join(line))
        _ = fp.write('\n')
