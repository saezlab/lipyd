#!/usr/bin/env python

# Igor Bulanov
# Denes Turei
# turei.denes@gmail.com


import lipyd.msproc_alternative as msprc




preproc = msprc.Preprocessing()


param_pp = {"signal_to_noise": 1.0}

pp_src_dir = ""
pp_dst_dir = ""

preproc.peak_picking(
                src = pp_src_dir,
                dst = pp_dst_dir,
                suffix_dst_files = "_picked",
                ext_dst_files = "mzML")

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

#We need to apply and update common parameters for each 3 process in Feature Finding
param_mtd.update(param_ff_com)
param_epd.update(param_ff_com)
param_ffm.update(param_ff_com)


param_ma = {"max_num_peaks_considered": 990}

ff_src_dir = "" 
ff_dst_dir = ""

preproc.feature_finding_metabo(
                            src = ff_src_dir,
                            dst = ff_dst_dir,
                            suffix_dst_files = "_feature",
                            ext_dst_files = "featureXML"
                            )

ma_src_dir = ""
ma_dst_dir = ""
reference_file = ""

preproc.map_alignment(src = ma_src_dir,
                            dst = ma_dst_dir,
                            suffix_dst_files = "_aligned",
                            ext_dst_files = "featureXML",
                            reference_file = reference_file
                            )    


preproc.setup_pp_params(**param_pp)

preproc.setup_feature_finding_metabo( **param_ff_com )
preproc.setup_ff_mtd_params(**param_mtd)
preproc.setup_ff_epd_params(**param_epd)
preproc.setup_ff_ffm_params(**param_ffm)
preproc.setup_ma_params()


preproc.run()

"""
Converation to MGF format

mzml_file = ""
mgf_file = ""

convert = Convert2mgf(mzml_file = mzml_file,
                mgf_file = mgf_file)

convert.convert()
"""


src_set = ""

preproc.get_sampleset(src = src_set)
