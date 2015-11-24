simpleTypes = [int, float, str, unicode]

charTypes = [str, unicode]

def uniqList(seq):
    # Not order preserving
    # from http://www.peterbe.com/plog/uniqifiers-benchmark
    keys = {}
    for e in seq:
        try:
            keys[e] = 1
        except:
            print e
            print seq
            print keys
    return keys.keys()

def flatList(lst):
    return [it for sl in lst for it in sl]

def delEmpty(lst):
    return [i for i in lst if len(i) > 0]

def uniqOrdList(seq, idfun = None): 
   # Order preserving
   # from http://www.peterbe.com/plog/uniqifiers-benchmark
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

def addToList(lst, toadd):
    if isinstance(toadd, list):
        lst += toadd
    else:
        lst.append(toadd)
    if None in lst:
        lst.remove(None)
    return uniqList(lst)

def something(anything):
    return not (anything is None or \
        (type(anything) in [list, set, dict, str, unicode] \
            and len(anything) == 0))
