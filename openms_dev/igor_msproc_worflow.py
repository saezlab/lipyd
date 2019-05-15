#!/usr/bin/env python
# -*- coding: utf-8 -*-

# EMBL 2019
# Igor Bulanov

from lipyd.msproc_alternative import *

def main():
    
    a = Preprocessing()
    
    param_pp = {"signal_to_noise": 1.0}
    
    param_ff_com = {"noise_threshold_int": 9.0,
                "chrom_peak_snr": 2.0,
                "chrom_fwhm": 4.0}
    
    param_mtd = {"mass_error_ppm" : 19.0,
                "reestimate_mt_sd": "true",
                "quant_method": "area"}
    
    
    param_epd = {"enabled": "true"}
    
    
    param_ffm = {"mz_scoring_13C": "true",
                "report_convex_hulls": "true",
                "remove_single_traces": "true"}

    # We need to apply and update common parameters for each 3 processes in Feature Finding
    param_mtd.update(param_ff_com)
    param_epd.update(param_ff_com)
    param_ffm.update(param_ff_com)

    param_ma = {"max_num_peaks_considered": 990}


    a.peak_picking(src = "/home/denes/archive/ltp/STARD10_invivo_raw/mzml/pos/.+\.mzML$",
                                dst = "/home/denes/archive/ltp/STARD10_invivo_raw/mzml/pos/picked/",
                                suffix_dst_files = "_picked",
                                ext_dst_files = "mzML",
                                )


    a.feature_finding_metabo(
                                src = "/home/denes/archive/ltp/STARD10_invivo_raw/mzml/pos/picked/.+\.mzML$",
                                dst = "/home/denes/archive/ltp/STARD10_invivo_raw/mzml/pos/featured/",
                                suffix_dst_files = "_feature",
                                ext_dst_files = "featureXML",
                                )


    a.map_alignment(src = "/home/denes/archive/ltp/STARD10_invivo_raw/mzml/pos/featured/.+\.featureXML$",
                                dst = "/home/denes/archive/ltp/STARD10_invivo_raw/mzml/pos/aligned",
                                suffix_dst_files = "_aligned",
                                ext_dst_files = "featureXML",
                                reference_file = "/home/denes/archive/ltp/STARD10_invivo_raw/mzml/pos/featured/150310_Popeye_MLH_AC_STARD10_A09_pos_picked_feature.featureXML",
                                ) 


    a.setup_pp_params(**param_pp)

    a.setup_ff_mtd_params(**param_mtd)
    a.setup_ff_epd_params(**param_epd)
    a.setup_ff_ffm_params(**param_ffm)
    a.setup_ma_params()



    a.run()

    p = "/home/denes/archive/ltp/STARD10_invivo_raw/mzml/pos/aligned/.+\.featureXML$"
    #a.get_data( p )
    s = a.get_sampleset(p)
    
    return s, a
