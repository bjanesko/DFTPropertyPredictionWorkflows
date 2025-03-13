import sys, os , re
a2k=627.5095 # kcal/mol per Hartree
RTfac = -2625.5/627.5095/(2.303*0.008314*298)
Gproton=-270.296
SIdata=dict()

def getG(f):
    G=0
    for line in open(f,'r'):
        if re.search('hermal Free',line):
           fields=line.split()
           G = float(fields[-1])
    return(G)

def getGeom(f):
    nat=0
    for line in open(f,'r'):
      if re.search('NAtoms=',line):
        fields=line.split()
        nat=int(fields[1])
    geom=''
    counter=0
    for line in open(f,'r'):
        if (re.search('Input orientation:',line) or re.search('Standard orientation:',line)):
            counter=nat+5
            geom=''
        if(counter>0):
            geom=geom+line
        counter=counter-1
    return(geom)

def getGmin(name,prot):
    Gmin=100
    Emin=100
    filemin=''
    confnum=0
    listconf= ['01','02','03','04','05']
    prot=int(prot)
    readdir='./%s/'%(name)
    #print('Checking dir ',readdir)
    for d in os.listdir(readdir):
        d2='%s/%s'%(readdir,d)
        if os.path.isdir(d2):
            #print('Checking dir ',d2)
            iprot=int(re.sub('_.*','',d))
            if(iprot == prot):
                E=0
                G=0
                Eminthis=0
                Gminthis=0
                Sminthis=''
                for i in listconf:
                    f='./%s/%s/%s.out'%(name,d,i)
                    if(os.path.isfile(f)):
                        for line in open(f,'r'):
                            if re.search('SCF Don',line):
                                fields=line.split()
                                E = float(fields[4])
                            if re.search('hermal Free',line):
                                fields=line.split()
                                G = float(fields[-1])
                        #print('Directory %s %s file %s  gives free energy %.6f '%(name,d,i,G))
                        if(G<Gmin):
                            Gmin = G
                            Emin = E
                            filemin = d
                            confnum = i
                        if(G<Gminthis):
                            Gminthis = G
                            Eminthis = E
                            Sminthis=getGeom(f)
                            filemin = d
                if(Gminthis<0):
                    SIdata[d2]='Most stable total energy (Ha), Gibbs free energy (Ha), and geometry for protomer %s\n E: %12.6f \n G: %12.6f\nGeometry:\n%s\n'%(d2,Eminthis,Gminthis,Sminthis)


    #print('Best free energy %s file %s  %.6f '%(filemin,confnum,Gmin))
    return(Gmin,d)

# Find the range of protonation states
protstates=[]
for name in os.listdir('./'):
    if(os.path.isdir(name)):
        protstates.append(int(name.split('_')[0]))
protstates=list(set(protstates))


# Calculate the protonation free energy DeltaG for each protonation state
for base in protstates[1:]:
    (GX,f1) = getGmin('./',base-1)
    (GHX,f2) = getGmin('./',base)
    DG0 = a2k*(GX-GHX)+Gproton  # HX  --> X(-) + H(+)
    calc0=-1.0*DG0*RTfac
    print('%20s %2d %6.2f  '%(name,base,calc0))

print('## Computed geometries and energies of the most stable structure at each protonation state')
for key in SIdata.keys():
    print(SIdata[key])
