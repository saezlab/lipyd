import pyopenms as oms

def extract_data(in_file_name, mz_search):
    xml_file = oms.FeatureXMLFile()
    fmap = oms.FeatureMap()
    xml_file.load(in_file_name, fmap)
    delta = 0.01
    #print( "FeatureMap size=", fmap.size() )
    for n in fmap:
        
        #if mz_search == n.getMZ():
        if abs(mz_search - n.getMZ()) < delta:
            """
            print( "mz=", n.getMZ(),
                    "rt=", n.getRT(),
                    " intensity=", n. getIntensity(),
                    "width=", n.getWidth(),
                    "charge=", n. getCharge() )
            """
            
            hull_list = n.getConvexHulls()  #getConvexHull() return ConvexHull2D;
            for hull in hull_list:
                hull_points = hull.getHullPoints()  # hull_points is numpy.ndarray;
                #print( "hull_points.size=", hull_points.size )
                for p in hull_points:
                    print( p[0], p[1] )
            
            subord_feature = n.getSubordinates()
            if subord_feature:
                #print("getSubordinates:")
                for f in subord_feature:
                    hull_list = f.getConvexHulls()
                    for hull in hull_list:
                        hull_points = hull.getHullPoints()  # hull_points is numpy.ndarray;
                        #print( "hull_points.size=", hull_points.size )
                        for p in hull_points:
                            print( p[0], p[1] )
        else:

            continue

if __name__ == "__main__":
    in_file_name = "/home/igor/Documents/Lipyd/src/lipyd_ms_preproc/150310_Popeye_MLH_AC_STARD10_A10_pos/STARD10_A10_pos_IterativePP.featureXML "
    extract_data(in_file_name, 704.52)
    
    """
    760.5822
    732.5533
    730.538
    784.5838
    746.6045
    772.6186
    718.5742
    704.5206
    """