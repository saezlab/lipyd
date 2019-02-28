import csv
import xlrd
from pyteomics.openms import featurexml

from lipyd import settings
from lipyd import sample
from lipyd import sampleattrs
from lipyd import lookup



if __name__ == "__main__":
    peaksfile = 'stard10_pos.csv'
    reader = sample.SampleReader(
        input_type = 'peaks',
        fname = peaksfile
    )

    samples = reader.get_sampleset(
        sampleset_args = {
            'sample_id_proc': sampleattrs.plate_sample_id_processor(),
        }
    )

    idx = samples.mzs_by_sample[
        :,samples.attrs.sample_id_to_index[('A', 10)]
    ].argsort()

    samples.sort_all(by = idx)

    mzs = samples.mzs_by_sample[:,samples.attrs.sample_id_to_index[('A', 10)]]
    rts = samples.rts[:,samples.attrs.sample_id_to_index[('A', 10)]]
    intens = samples.intensities[:,samples.attrs.sample_id_to_index[('A', 10)]]

    file_name = "feature_map.featureXML"
    out_iterator = featurexml.read(source  = file_name,
                                    read_schema=False)
    
    dim_name    = "dim"
    pos_name    = "position"
    intens_name = "intensity"
    intensity   = 0.0
    rt_dim_value = 0
    mz_dim_value = 1
    RT = 0.0
    MZ = 0.0
    feature_list = []
    
    for dic in out_iterator:
        for k,v in dic.items():
            if k == pos_name:
                dic_list = v
                if dic_list[0][dim_name] == rt_dim_value: #was 0;
                    RT = dic_list[0][pos_name]
                else:
                    MZ = dic_list[0][pos_name]
                    
                if dic_list[1][dim_name] == mz_dim_value: #was 1;
                    MZ = dic_list[1][pos_name]
                else:
                    RT = dic_list[1][pos_name]
            if k == intens_name:
                intensity = v
            else:
                continue
            # RT in seconds convert into minutes:
            feature_list.append( (float(MZ), float(RT)/60., float(intensity)) )

    feature_list.sort()
    
    intersection_dic = {}

    for e in feature_list:
        #closest_index = lookup.find(
        result_list = lookup.findall(
            samples.mzs_by_sample[
                :,samples.attrs.sample_id_to_index[('A', 10)]
            ], # all masses in the sample
            e[0],        # mass to search for
            t = 10 # tolerance in ppm
        )
        if len( result_list ):
            for closest_index in result_list:
                intersection_dic[(e[0], e[1], e[2])] = \
                    ((mzs[closest_index],
                    rts[closest_index],
                    intens[closest_index],
                    (e[0]/mzs[closest_index]-1)*1e6 )
                )
            
        else:
            intersection_dic[(e[0], e[1], e[2])] = (result_list, None)
      

    for k, v in intersection_dic.items():
        print(k, v)


    
    
    
