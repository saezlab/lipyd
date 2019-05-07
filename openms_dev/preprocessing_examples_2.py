#!/usr/bin/env python

# Igor Bulanov
# Denes Turei
# turei.denes@gmail.com

preproc = Preprocessing()


param_pp = {"signal_to_noise": 1.0}

preproc.peak_picking(src = "/home/igor/Documents/Scripts/Data/Raw_data_STARD10/.+\.mzML$",
                            dst = "/home/igor/Documents/Scripts/Data/Picked_STARD10",
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
param_ma = {"max_num_peaks_considered": 990}

preproc.feature_finding_metabo(
                            src = "/home/igor/Documents/Scripts/Data/Picked_STARD10/.+\_picked\.mzML$",
                            dst = "/home/igor/Documents/Scripts/Data/feature_STARD10/",
                            suffix_dst_files = "_feature",
                            ext_dst_files = "featureXML"
                            )

preproc.map_alignment(src = "/home/igor/Documents/Scripts/Data/feature_STARD10/.+\.featureXML$",
                            dst = "/home/igor/Documents/Scripts/Data/Aligned_STARD10/",
                            suffix_dst_files = "_aligned",
                            ext_dst_files = "featureXML",
                            reference_file = "/home/igor/Documents/Scripts/Data/feature_STARD10/150310_Popeye_MLH_AC_STARD10_A09_pos_picked.featureXML"
                            )    

preproc.run()

"""
convert = Convert2mgf(mzml_file = "/home/igor/Documents/Scripts/Data/Picked_STARD10/150310_Popeye_MLH_AC_STARD10_A11_pos_picked.mzML",
                mgf_file = "/home/igor/Documents/Scripts/Data/Picked_STARD10/150310_Popeye_MLH_AC_STARD10_A11_pos_picked.mgf")

convert.convert()
"""
preproc.get_data()
p2.run()
