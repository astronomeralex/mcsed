# -*- coding: utf-8 -*-

# =========================================================================== #
import numpy as np

# =========================================================================== #
def fracsum(num, denom):
    """
    this functions returns sum(top/bottom) where top and bottom are 1d arrays.
    if bottom == 0 at a index, that index is skipped when creating the sum.
    
    its basically a smarter np.nansum
    """
    result = 0.0
    assert len(np.where(denom < 0)[0])==0 #denom should be positive or 0
    goodIDXs = np.where(denom != 0)[0]
    for i in range(0, len(goodIDXs)):
        idx = goodIDXs[i]
        result = result + num[idx]/denom[idx]
        
    """
    # for debugging
    print
    print num
    print
    print denom
    print
    print 'RESULT: ', result
    print "=========================================="
    """
    
    return result
    
# =========================================================================== #
def cartesian(arrays, out=None):
    '''
    Magic function that takes arrays and returns all possible combinations of their 
    entries.  This function is responsable for allowing me to go from parameterPointsRange
    to all the unique combinations that are pluged into python-fsps.
    
    example: input=[ [0,1,2]    would yield output=[ [0,7]
                     [7,8,9] ]                       [0,8]
                                                     [0,9]
                                                     [1,7]
                                                     [1,8]
                                                     [1,9]
                                                     [2,7]
                                                     [2,8]
                                                     [2,9] ]
    
    CODE FROM:  http://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays
    '''
    import numpy as np
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype
    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)
    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out


# =========================================================================== #    
def getMode(array, decimal_places):
    '''More intellegent way to get the mode of a 1d array'''
    ''' THIS MIGHT BE BROKEN '''
    # round data
    array = np.around(array, decimal_places)
    # get unique entries
    unique = np.unique(array)
    # create count array
    count = np.zeros(len(unique))
    # loop through unique entries, geting their number of occurence
    for i in range(0, len(unique)):
        count[i] = len(np.where(array == unique[i])[0])
    # modes are values where count = max(count), get them all
    modeIDXs = np.where(count == np.max(count))[0]
    return unique[modeIDXs], count[modeIDXs]
    
# =========================================================================== #    
