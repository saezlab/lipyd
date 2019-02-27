import pyopenms as oms

def extract_data(in_file_name, mz_search):
    xml_file = oms.FeatureXMLFile()
    fmap = oms.FeatureMap()
    xml_file.load(in_file_name, fmap)
    delta = 0.01
    #print( "FeatureMap size=", fmap.size() )
    for n in fmap:
        
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
    extract_data()