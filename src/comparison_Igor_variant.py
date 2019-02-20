# -- 2019/02/05--
# ver.5
# feature_map.featureXML
# stard10_pos.csv
# STARD10_top_features_MLH_final.xlsx
# 
import csv
import xlrd
from pyteomics.openms import featurexml

def is_number(str):
    try:
        float(str)
        return True
    except ValueError:
        return False

def get_point(the_list, mz, delta_mz):
    mz_min = mz - delta_mz/2.
    mz_max = mz + delta_mz/2.
    #rt_min = rt - delta_rt/2.
    #rt_max = rt + delta_rt/2.

    for e in the_list:
        if mz_min <= e[0] and e[0] <= mz_max: #and \
            #rt_min <= e[1] and e[1] <= rt_max:
            return e

    return None

if __name__ == "__main__":
    csv_file_name = "stard10_pos.csv"
    mz_col = 16     #stard10_a10_pos m/z,
    rt_col = 17     #stard10_a10_pos RT mean,
    intens_col = 18 #stard10_a10_pos Normalized Area,
    csv_list = []   # list of tuple: (m/z, RT, Normalized Area)

    with open(csv_file_name, "r") as f:
        reader = csv.reader(f, delimiter="," )
        next(f) # skip header line;
        for row in reader:
            #print( row[mz_col], row[rt_col], row[intens_col] )
            if is_number( row[mz_col] ) and \
                is_number( row[rt_col] ) and \
                is_number( row[intens_col] ):
                    csv_list.append( (float(row[mz_col]), float(row[rt_col]), float(row[intens_col])) )
            

    file_name = "STARD10_A10_feature_map.featureXML"
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
                #key= position val= [{'dim': 0, 'position': 332.445618411506}, {'dim': 1, 'position': 262.164012851789}]
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

    csv_list.sort()
    feature_list.sort()
    
    total_list = []
    delta_mz = .01
    #delta_rt = 1.
    
    for c in csv_list:
        point_from_feature = get_point( feature_list, 
                                    c[0],
                                    #c[1],
                                    delta_mz,
                                    #delta_rt 
                            )
        if point_from_feature:
            ppm = ((point_from_feature[0] / c[0]) - 1) * 1e6
            total_list.append( (c, point_from_feature, ppm) )                       
        else:
            total_list.append( (c, point_from_feature) )

    
    for e in total_list:
        if e[1] != None: 
            print("{}, {}, {}, {}, {}, {}, {}".format(e[0][0], e[0][1], e[0][2],
                e[1][0], e[1][1], e[1][2], e[2] ))
        else:
            print(e[0][0], e[0][1], e[0][2])
    

