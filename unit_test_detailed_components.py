
class miscellaneous_things():
    def __init__(self):
        self.wave_list=['f-', 'a-','s-','c','f+','a+','s+']
        self.field_list = ['d','px','py','pz','hx','hy','hz','p']
        self.form_s = " %5s"
        self.form_f = " %5.2f"

t=miscellaneous_things()
if 1:
    print("Wave  Content hat(FSA)")
    for w in t.wave_list:
        x1 = read_stuff['waveset'].wave_content[w]
        x2 = ic_stuff['waveset'].wave_content[w]
        x= np.abs(x1-x2)
        val = x.sum()
        print("%8s %0.3e"%(w,val))
if 0:
    print("Real frame: XYZ")
    for f in t.field_list:
        x1=read_stuff['cubes'][f]
        x2=ic_stuff['cubes'][f]
        x= np.abs(x1-x2)
        #val = np.std(x)/np.abs(x1.max())
        val=x.sum()

        print("%8s %0.3e"%(f,val))

if 1:
    print("Wave frame: hat(ABC)")
    for f in t.field_list:
        x1=read_stuff['waveset'].wave_frame[f]
        x2=ic_stuff['waveset'].wave_frame[f]
        x= np.abs(x1-x2)
        val=x.sum()
        #val = np.std(x)/np.abs(x1.max())
        print("%8s %0.3e"%(f,val))

if 0:
    print("Other wave frame, hat(ABC) in hat_system.quan")
    for f in read_stuff['waveset'].hat_system.quan.keys():
        x1=read_stuff['waveset'].hat_system.quan[f]
        x2=ic_stuff['waveset'].hat_system.quan[f]
        x= np.abs(x1-x2)
        val=x.sum()

        print("%12s %0.3e"%(f,val))#/np.abs(x1.max())))

if 0:
    print("B0hat")
    for dim in [0,1,2]:
        x1=read_stuff['waveset'].B0hat[dim,...]
        x2=ic_stuff['waveset'].B0hat[dim,...]
        x= np.abs(x1-x2)
        print("%12s %0.3e"%(str(dim),np.std(x)))#/np.abs(x1.max())))
if 0:
    print("THING 1")

    x1=read_stuff['waveset'].hat_system.thing1
    x2=ic_stuff['waveset'].hat_system.thing1
    x= np.abs(x1-x2)
    var = x.sum()/x1.size
    print("%5s    %0.3e"%('thing1',var))#/np.abs(x1.max())))

if 0:
    print("hat(ABC) left vectors")
    for w in read_stuff['waveset'].hat_system.left.keys():
        for f in t.field_list:
            x1=read_stuff['waveset'].hat_system.left[w][f]
            x2=ic_stuff['waveset'].hat_system.left[w][f]
            x= np.abs(x1-x2)
            print("%5s %5s   %0.3e"%(w,f,np.sum(x)))#/np.abs(x1.max())))

if 0:
    print("Wave frame: Aunit")
    x1=read_stuff['waveset'].a_unit
    x2=ic_stuff['waveset'].a_unit
    x= np.abs(x1-x2)
    var=x.sum()
    #var = np.std(x)/np.abs(x1.max())))
    print("%8s %0.3e"%(f,var))
    print("Wave frame: Bunit")
    x1=read_stuff['waveset'].b_unit
    x2=ic_stuff['waveset'].b_unit
    x= np.abs(x1-x2)
    var=x.sum()
    #var = np.std(x)/np.abs(x1.max())))
    print("%8s %0.3e"%(f,var))
    print("Wave frame: C unit")
    x1=read_stuff['waveset'].c_unit
    x2=ic_stuff['waveset'].c_unit
    x= np.abs(x1-x2)
    var=x.sum()
    #var = np.std(x)/np.abs(x1.max())))
    print("%8s %0.3e"%(f,var))

