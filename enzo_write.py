"""
Let's write an enzo file, for restarting with 
TurbulenceSimulation
"""
import h5py
def dump_h5(array,filename):
    """presently this will not do multiple-component sets."""


    setname = filename.split("/")[-1]
    GridRank = len(array.shape)

    if GridRank >1:
        to_dump = array.swapaxes(0,GridRank-1)
    #to_dump = array


    #u'Component_Rank', array([1], dtype=int32))
    #u'Component_Size', array([2097152], dtype=int32))
    #u'Rank', array([3], dtype=int32))
    #u'Dimensions', array([64, 64, 64], dtype=int32))
    fptr = h5py.File(filename,'w')
    try:
        fptr.create_dataset(setname,data=to_dump)
        fptr[setname].attrs['Component_Rank'] = 1
        fptr[setname].attrs['Component_Size'] = to_dump.size
        fptr[setname].attrs['Rank'] = GridRank
        fptr[setname].attrs['Dimensions'] = array.shape #to_dump.shape


        #for a in attrs:
        #    fptr[setname].attrs[a]=attrs[a]
    except:
        raise
    finally:
        fptr.close()



if 0:
    setname = fn.split("/")[-1]
    fptr = h5py.File(fn,'r')
    print("ATTRS")
    for a in fptr[setname].attrs:
        print(a,fptr[setname].attrs[a])
    fptr.close()
