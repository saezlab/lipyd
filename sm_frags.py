
import os
from lipyd import ms2
from lipyd import lipproc

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
    'SM [Ch+P+C2NH] (225.0999)',
    'SM [Ch+P+C2NH+CO] (253.0948)',
)

def check_sm_frags(scan):
    
    for frag in sm_frags:
        
        print(frag, scan.has_fragment(frag))


mgfdir_invitro = '/home/denes/archive/ltp/mgf_invitro'
mgfname_invitro = '161203_MLH_EMV_GLTPD1_pos_%s.mgf'
fractions_invitro = ('E07', 'E08', 'E09', 'E10')

mgfdir_invivo = '/home/denes/archive/ltp/mgf_invivo'
mgfname_invivo = '150628_Popeye_MLH_AC_GLTPD1_pos_%s.mgf'
fractions_invivo = ('A11', 'A12')

mgfs_invitro = dict(
    (fr, os.path.join(mgfdir_invitro, mgfname_invitro % fr))
    for fr in fractions_invitro
)

mgfs_invivo = dict(
    (fr, os.path.join(mgfdir_invivo, mgfname_invivo % fr))
    for fr in fractions_invivo
)

mzs_invitro = (
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

rts_invitro = (
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

mzs_invivo = (
    703.575398430908,
    813.682274069939,
    787.666910609034,
    785.652436668606,
    701.560224499535,
    731.606735812837,
    811.669100112695,
    675.543961050273,
    759.637273205114,
    799.667855362805,
    815.700047967678,
    689.55967974059,
    705.590672360987,
    729.591161886638,
    717.554817590802,
)

rts_invivo = (
    (8.36, 9.02),
    (12.58, 13.34),
    (12.62, 13.18),
    (11.22, 12.36),
    (7.48, 8.29),
    (9.78, 10.37),
    (11.49, 12.12),
    (7.26, 7.81),
    (11.18, 11.78),
    (11.87, 12.77),
    (14.02, 14.62),
    (7.83, 8.37),
    (8.43, 9.54),
    (8.70, 9.54),
    (7.13, 7.72),
)


hdr = [
    'mz', 'rt', 'scan_id',
    'deltaRT', 'frac',
    'MS1',
    'score',
    '[Sph-2xH2O+H]+',
] + list(sm_frags)
hdr.append('chain_frags')

result = []

def collect_scan_lines(mzs, rts, mgfs):
    
    result = []
    
    for mz, rt in zip(mzs, rts):
        
        fe = ms2.MS2Feature(mz, ionmode = 'pos', mgfs = mgfs, rt = rt)
        fe.build_scans()
        
        for sc, drt in zip(fe.scans, fe.deltart):
            
            ids = sc.identify()
            
            smrec = list(
                sc.get_ms1_records(hg = 'SM', databases = {'lipyd.lipid'})
            )
            
            chainid = set()
            
            for sm in smrec:
                
                ccomb = sc.chain_combinations(sm[0])
                
                for cc in ccomb:
                    
                    chainid.add('%s[full=%s,frags=[%s,%s],i=[%.03f,%.03f]]' % (
                        sm[0].summary_str(),
                        lipproc.full_str(sm[0].hg, cc[0]),
                        cc[1].fragtype[0],
                        cc[1].fragtype[1],
                        cc[1].i[0] * 100,
                        cc[1].i[1] * 100,
                    ))
            
            for ms1id, ms2ids in ids.items():
                
                if ms2ids and ms2ids[0].hg.main == 'SM':
                    
                    this_line = [
                        '%.07f' % mz,
                        '%.02f - %.02f' % rt,
                        '%u' % sc.scan_id,
                        '%.02f' % drt,
                        sc.sample,
                        ms1id,
                        '%u' % ms2ids[0].score_pct,
                    ]
                    
                    isph = min(
                        sc.fragments_by_chain_type(frag_type = 'Sph-2xH2O+H'),
                        default = None
                    )
                    
                    irel_sph = (
                        '%.03f' % (sc.inorm[isph] * 100)
                        if isph is not None
                        else ''
                    )
                    
                    this_line.append(irel_sph)
                    
                    for frag in sm_frags:
                        
                        i = sc.fragment_by_name(frag)
                        reli = (
                            '%.03f' % (sc.inorm[i] * 100)
                            if i is not None
                            else ''
                        )
                        
                        this_line.append(reli)
                    
                    this_line.append(';'.join(chainid))
                    
                    result.append(this_line)
    
    return result

result.extend(collect_scan_lines(mzs_invitro, rts_invitro, mgfs_invitro))
result.extend(collect_scan_lines(mzs_invivo,  rts_invivo,  mgfs_invivo))

with open('GLTPD1_pos_SM.csv', 'w') as fp:
    
    _ = fp.write('\t'.join(hdr))
    _ = fp.write('\n')
    
    for line in result:
        
        _ = fp.write('\t'.join(line))
        _ = fp.write('\n')
