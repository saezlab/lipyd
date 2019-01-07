import pyopenms as pms

if __name__ == "__main__":
    
    raw_file_name = ""
    in_map = pms.MSExperiment()
    out_map = pms.MSExperiment()
    pms.MzMLFile().load(raw_file_name, in_map)
    
    pp = pms.PeakPickerCWT()
    param = pms.Param()
    param.setValue("peak_width", 0.1)
    pp.setParameters(param);
    pp.pickExperiment(in_map, out_map);
    pp.pickExperiment(in_map, out_map);
    out_map.updateRanges();
    
    pms.MzMLFile().store("first_result.mzML", out_map)
