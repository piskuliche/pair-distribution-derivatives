import numpy as np


def Conv_Main(Iargs):
    header,sections=Read_Data(Iargs)
    Gen_Configs(sections,Iargs)

def Read_Data(Iargs):
    """This function takes a data file and saves the important bits, 
    but also clears the non-important bits."""
    def _Pull_Sections(lines):
        # Grab the data file header
        header = []
        for line in lines:
            if "Masses" in line:
                break
            header.append(line)
        nheader = len(header)
        lines=lines[nheader:]
        # Grab the sections
        flag=0
        sections={}
        sec = None
        for line in lines:
            cleanline = line.strip().split()
            if len(cleanline)==1 and cleanline != "\n":
                sec=cleanline[0]
                sections[sec]=[]
            if sec is not None and len(cleanline) > 1:
                sections[sec].append(line)
        return header,sections

    # Function Start
    with open(Iargs.datafile,'r') as f:
        lines=f.readlines()
        header,sections=_Pull_Sections(lines)
        print("Found the following keys:")
        for key in sections:
            print("-"+key+" has %d lines" % len(sections[key]))
    return header, sections

def Gen_Configs(sections, Iargs):
    def _Set_Map(sections):
        masslines=sections["Masses"]
        M={}
        for mass in masslines:
            mvals=mass.strip().split()
            M[int(mvals[0])]=float(mvals[1])
        atomlines=sections["Atoms"]
        mid={}
        typ={}
        chg={}
        for atom in atomlines:
            avals=atom.strip().split()
            mid[avals[0]]=int(avals[1])
            typ[avals[0]]=int(avals[2])
            chg[avals[0]]=float(avals[3])
        return M,mid

    print("Beginning Mapping Process")
    M,mid = _Set_Map(sections)
    natoms=len(mid)
    nmols=None
    with open(Iargs.dumpfile, 'r') as df:
        for i in range(Iargs.nconfigs):
            if i%100 == 0:
                print("Reached config %d of %d" % (i,Iargs.nconfigs))
            L,r,ts,mol,typid=Read_Config(df,natoms)
            if i == 0:
                Prep_Dir(mol)
            if i%skip == 0:
                comr=Calc_Com(r,mol,typid,M)
                nmols=len(np.unique(mol))
                dr = Get_Distances(comr,L)
                #idx,sdr=Sort_Distances(dr)
                Write_Groups(dr,ts,cut)
    Write_System(sysfile,nmols)
    return

def Read_Config(df,natoms):
    L=[0,0,0]
    r=[]
    molid=[]
    typval=[]
    df.readline()
    ts=int(df.readline().strip())
    df.readline()
    df.readline()
    df.readline()
    line = df.readline().strip().split()
    L[0] = -float(line[0])+float(line[1])
    line = df.readline().strip().split()
    L[1] = -float(line[0])+float(line[1])
    line = df.readline().strip().split()
    L[2] = -float(line[0])+float(line[1])
    df.readline()
    for l in range(natoms):
        line = df.readline().strip().split()
        mol,typ,x,y,z = int(line[1]),int(line[2]),float(line[3]),float(line[4]),float(line[5])
        r.append([x,y,z])
        molid.append(mol)
        typval.append(typ)
    return np.array(L),np.array(r),ts,molid,typval


def Get_Distances(r,L):
    """
    This is a pretty complicated function to do a simple thing.
    This code takes two vectors of size(m,3) and size(n,3)
    Then it uses numpy broadcasting to calculate ALL the pairwise distances and outputs them in a matrix, sized (m,n)
    Then I use it to calculate the distances.
    (described here: https://stackoverflow.com/questions/60039982/numpy-python-vectorize-distance-function-to-calculate-pairwise-distance-of-2-ma/60040269#60040269)
    """
    vecdr = r[:,np.newaxis,:2]-r[np.newaxis,:,:2]
    vecdr = vecdr - np.multiply(L[:2],np.round(np.divide(vecdr,L[:2])))
    dr = np.linalg.norm(vecdr,axis=-1)
    return dr

def Sort_Distances(dr):
    """
    This function takes an array (dr) and sorts it by the last axis.

    It returns:
    1) the idx map of the sort
    2) the sorted dr array (sdr)
    """
    natoms = np.shape(dr)[0]
    idx = np.argsort(dr,axis=-1,kind='quicksort')
    #print(idx)
    sdr=np.take_along_axis(dr, idx, axis=-1)
    return idx,sdr

def Prep_Dir(mol,kspace=False):
    import os
    if not os.path.exists("./include_dir"):
        os.makedirs("./include_dir")
    if not os.path.exists("./ener_dir"):
        os.makedirs("./ener_dir")
    if not os.path.exists("./calc_dir"):
        os.makedirs("./calc_dir")
    if not os.path.exists("./calc_dir/logs"):
        os.makedirs("./calc_dir/logs")
    nmols=len(np.unique(mol))
    for i in range(nmols):
        fi = open("./include_dir/include.groups-%d"%i,'w')
        fi.write("group solu molecule 0\n")
        fi.write("group close molecule 0\n")
        fi.write("group far molecule 0\n")
        fi.write("compute sske solu ke\n")
        fi.write("compute sspair solu group/group solu pair yes\n")
        if kspace==True: fi.write("compute sskspc solu group/group solu kspace yes\n")
        fi.write("compute ssbnd solu bond\n")
        fi.write("compute ssang solu angle\n")
        #fi.write("compute ssdih solu dihedral\n")
        #fi.write("compute ssimp solu improper\n")
        #fi.write("variable ssintra equal c_ssbnd+c_ssang+c_ssdih+c_ssimp\n")
        fi.write("variable ssintra equal c_ssbnd[1]+c_ssang[1]\n")
        fi.write("compute ccpair close group/group close pair yes\n")
        if kspace==True: fi.write("compute cckspc close group/group close kspace yes\n")
        fi.write("compute ccke close ke\n")
        fi.write("compute ccbnd close bond\n")
        fi.write("compute ccang close angle\n")
        #fi.write("compute ccdih close dihedral\n")
        #fi.write("compute ccimp close improper\n")
        fi.write("variable ccintra equal c_ccbnd[1]+c_ccang[1]\n")
        fi.write("compute ffpair far group/group far pair yes\n")
        if kspace==True: fi.write("compute ffkspc far group/group far kspace yes\n")
        fi.write("compute ffke far ke\n")
        fi.write("compute ffbnd far bond\n")
        fi.write("compute ffang far angle\n")
        #fi.write("compute ffdih far dihedral\n")
        #fi.write("compute ffimp far improper\n")
        #fi.write("variable ffintra equal c_ffbnd+c_ffang+c_ffdih+c_ffimp\n")
        fi.write("variable ffintra equal c_ffbnd[1]+c_ffang[1]\n")
        fi.write("compute sfpair solu group/group far pair yes\n")
        if kspace==True: fi.write("compute sfkspc solu group/group far kspace yes\n")
        fi.write("compute scpair solu group/group close pair yes\n")
        if kspace==True: fi.write("compute sckspc solu group/group close kspace yes\n")
        fi.write("compute cfpair close group/group far pair yes\n")
        if kspace==True: fi.write("compute cfkspc close group/group far kspace yes\n")
        fi.write("fix 1 all ave/time 1000 1 1000 c_sske c_ccke c_ffke file ../ener_dir/ke.compute%d\n"%i)
        fi.write("fix 2 all ave/time 1000 1 1000 c_sspair c_ccpair c_ffpair c_scpair c_sfpair c_cfpair file ../ener_dir/pair.compute%d\n"%i)
        if kspace==True: fi.write("fix 3 all ave/time 1000 1 1000 c_sskspc c_cckspc c_ffkspc c_sckspc c_sfkspc c_cfkspc file ../ener_dir/kspc.compute%d\n"%i)
        fi.write("fix 4 all ave/time 1000 1 1000 v_ssintra v_ccintra v_ffintra file ../ener_dir/intra.compute%d\n"%i)
        fi.close()






def Write_Groups(dr,ts,cut):
    nmols=np.shape(dr)[0]
    for i in range(nmols):
        fi = open("./include_dir/include.groups-%d"%i,'a')
        close=np.where(dr[0]<cut)[-1]
        fi.write("group lt molecule ")
        for j in close:
            fi.write("%d "% (j+1))
        fi.write("\n")
        fi.write("group solu molecule %d\n"%(i+1))
        fi.write("group close subtract lt solu\n")
        fi.write("group far subtract all lt solu\n")
        fi.write('print "beginning rerun %d"\n'%ts)
        fi.write("rerun ../dump.lammpsdump first %d last %d dump x y z\n" % (ts,ts))
        fi.write('print "finishing rerun"\n')
        fi.close()

def Calc_Com(r,mol,typ,M):
    natoms=len(mol)
    nmols=len(np.unique(mol))
    comr=np.zeros((nmols,3))
    Mtot=np.zeros(nmols)
    for i in range(natoms):
        mid=mol[i]-1
        atype=typ[i]
        for j in range(3):
            comr[mid][j] += r[i][j]*M[atype]
        Mtot[mid]+=M[atype]
    for i in range(nmols):
        comr[i]=comr[i]/Mtot[i]
    return comr 

def Gen_Configs(datafile,nconfigs,sections,sysfile,cut,skip):
    print("Beginning Mapping Process")
    M,mid = set_map(sections)
    natoms=len(mid)
    nmols=None
    with open(datafile, 'r') as df:
        for i in range(nconfigs):
            if i%100 == 0:
                print("Reached config %d of %d" % (i,nconfigs))
            L,r,ts,mol,typid=Read_Config(df,natoms)
            if i == 0:
                Prep_Dir(mol)
            if i%skip == 0:
                comr=Calc_Com(r,mol,typid,M)
                nmols=len(np.unique(mol))
                dr = Get_Distances(comr,L)
                #idx,sdr=Sort_Distances(dr)
                Write_Groups(dr,ts,cut)
    Write_System(sysfile,nmols)
    return

def Write_System(sysfile,nmols):
    print("Writing System Files")
    lines=None
    with open(sysfile,'r') as f:
        lines = f.readlines()
    for i in range(nmols):
        with open("calc_dir/"+sysfile+"-%d"%i,'w') as g:
            for line in lines:
                g.write(line)
            g.write("log logs/log.%d.out\n"%i)
            g.write("include ../include_dir/include.groups-%d\n"%i)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='This code creates a series of files for lammps simulations that quickly rerun and re-calculate energies')
    parser.add_argument('-data',    default="equil.data",       type=str,help='Name of data file [default equil.data]')
    parser.add_argument('-dump',    default="dump.lammpsdump",  type=str,help='Name of dump file [default dump.lammpsdump]')
    parser.add_argument('-nc',      default=10000,              type=int, help='Number of configurations [default 10000]')
    parser.add_argument('-sf',      default="system.in",        type=str, help='System input file [default system.in]')
    parser.add_argument('-cut',     default=30,                 type=float, help='Cutoff to use for defining what is close')
    parser.add_argument('-skip',    default=1,                  type=int, help='Frequency to use frames')
    Iargs = parser.parse_args()

    datafile = args.data
    dumpfile = args.dump
    nconfigs = args.nc
    sysfile  = args.sf
    cut      = args.cut
    skip     = args.skip
    
    Conv_Main(Iargs):
    header,sections=Read_Data(datafile,Iargs)
    Gen_Configs(dumpfile,nconfigs,sections,sysfile,cut,skip)
