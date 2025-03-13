import itertools
import os ,sys , shutil
import argparse

# This script takes as input the XYZ coords of a 'fully protonated' 3D structure,
# plus the charge of the fully protonated form, a list of removable proton
# indices, and a number of explicit water molecules. It writes CREST
# metadynamics inputs for all possible protonation states

CLI=argparse.ArgumentParser()
CLI.add_argument('filename')
CLI.add_argument('--qfull',type=int)
CLI.add_argument('--nsolv',type=int)
CLI.add_argument('--prots',nargs="*",type=int)
args=CLI.parse_args()
filename=args.filename
qfull=int(args.qfull)
nsolv=int(args.nsolv)
theprots=args.prots

if(not(os.path.isfile(filename))):
    exit('No file %s'%(filename))
if(not(isinstance(theprots,list))):
    exit('This is not a list of protons: ',theprots)
if(nsolv<0):
    exit('Needs more than zero explicit waters, not ',nsolv)

nprot =len(theprots)
f=open(filename,"r")
data0=f.read()
f.close()
ats=data0.split('\n')
natall=len(ats)-1
print('Fully protonated molecule has %d atoms'%(natall))
tag=filename.replace('.xyz','')
newdir0="%s_%dw"%(tag,nsolv)
if not os.path.exists(newdir0):
   os.mkdir(newdir0)
newfile="%s/data.txt" %(newdir0)
f=open(newfile,"w")
f.write('Charge: %d\nExchangable protons: %s\nGeometry: %s\n'%(qfull,theprots,filename))

for ilost in range(nprot+1):
    iprot=nprot-ilost
    q=qfull- ilost
    nat = natall - ilost
    subsets = list(itertools.combinations(theprots,ilost))
    for iset in range(len(subsets)):
        subset = subsets[iset]
        print(subset)
        newdir = "%s/%d_%d" %(newdir0,iprot,iset+1)
        if not os.path.exists(newdir):
            os.mkdir(newdir)

        # Write the solute and solvent files and a dummy GJF for viewing
        newfile="%s/struc.xyz" %(newdir)
        newfile2="%s/struc.gjf" %(newdir)
        f=open(newfile,"w")
        f2=open(newfile2,"w")
        f.write("%d\n\n" % (nat))
        f2.write('#P hf 3-21g \n\nTEST\n\n%d 1\n'%(q))
        for iline in range(natall): # Skip the first two lines
            if not iline in subset:
                f.write("%s \n" % (ats[iline]))
                f2.write("%s \n" % (ats[iline]))
        f2.write('\n\n')
        f.close()
        f2.close()
        if(nsolv>0):
            newfile='%s/h2o.xyz'%(newdir)
            shutil.copyfile('h2o.xyz',newfile)

        # Write the batch script
        newfile="%s/crest.slurm" %(newdir)
        f=open(newfile,"w")
        f.write("""#!/bin/tcsh
#SBATCH --job-name=crest
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20GB
#SBATCH --time=48:00:00
#SBATCH --output=ensemble_%j.o
#SBATCH --error=ensemble_%j.e
pwd
hostname
whoami
date
""")
        st2='crest struc.xyz --gfn2 --gbsa h2o -T 6 --chrg %d   --noreftopo ' %(q)
        if(nsolv>0):
            st2='crest struc.xyz --qcg h2o.xyz --nsolv %d --ensemble --mdtime 50 --nofix --gfn2 --gbsa h2o -T 6 --chrg %d   --noreftopo ' %(nsolv,q)
        f.write(st2)
        f.write("\n\n")
        f.close()

