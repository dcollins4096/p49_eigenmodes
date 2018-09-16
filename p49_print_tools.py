import numpy as np

sfrm = "%12s"
if 'thresh' not in dir():
    thresh = 1e-11
def nz(arr):
    return  np.where(np.abs(arr) > thresh)
def nonzero(arr):
    return arr[ nz(arr)]
def pnz(arr):
    print(arr[ nz(arr)])
    print(nz(arr))
cf = '({0.real:5.2e} + {0.imag:5.2e}i)'
def ampl_str(arr):
    out = ""
    for v in arr:
        out += cf.format(v)
    return out
def maxis(inarr,dim):
    arr=np.zeros_like(inarr)
    ok =  inarr>1e-13
    arr[ok] = inarr[ok]
    out = np.log(np.max(np.abs(arr),axis=dim))

    #out = np.sum(np.abs(inarr),axis=dim)
    #out[ out==0] == -1
    return out
def print_nz_2(arr):
    nzb = nz(arr)
    nozob = nonzero(arr)
    for n,p in enumerate(zip(*nzb)):
        print(p, cf.format(nozob[n]))
