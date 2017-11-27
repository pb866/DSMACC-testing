#!/usr/local/anaconda/bin/python
import re,os,sys,glob

# Find all KPP input files
filedata = ''.join(tuple(open('./model.kpp')))
mech = re.findall(r'[^/\h]#INCLUDE(\s*.+)' ,filedata,re.IGNORECASE)

# Move KPP input files to save folder
for i in mech:
    i = re.sub(r'//.+$','',i)
    os.system('cp %s ./save/exec/%s/' %(i,sys.argv[1]))

# Move model.kpp and executable
os.system('cp ./model.kpp ./save/exec/%s/' %(sys.argv[1]))
os.system('cp ./model ./save/exec/%s/' %(sys.argv[1]))

# List all nc files
nc = glob.glob('*.nc')
print "Choose nc file to save:\n"
print "0  -  Save no nc files"
for i,f in enumerate(nc): print i+1 , ' - ', f.replace('.nc','')
selected_nc = filter(lambda x: len(x)>0 ,   raw_input('Enter Number(s): ').split(' '))
# Adjust file index
file_idx = []
for i in selected_nc:
    file_idx.append(int(i)-1)

# Rename and save first file
print "cp %s ./save/results/%s.nc"%(nc[file_idx[0]],sys.argv[1])
if file_idx[0] >= 0: os.system("cp %s ./save/results/%s.nc"%(nc[file_idx[0]],sys.argv[1]))
# Rename all other files using time step (or index if obsolete) and save
n = 0
c = 0
for i in range(1,len(file_idx)):
    f = nc[file_idx[i]]
    idx = f.rfind("_")
    if idx > 0 and f[:5] != "ropa_":
        ofile = sys.argv[1]+f[idx:]
    elif idx > 0 and f[:5] == "ropa_":
        if c==0:
            c += 1
            ofile = "ropa_"+sys.argv[1]+".nc"
        else:
            ofile = "ropa_"+sys.argv[1]+f[idx:]
    else:
        n += 1
        idx = f.rfind(".")
        ofile = "%s_%d%s"%(sys.argv[1],n,f[idx:])
    os.system("cp %s save/results/%s"%(f,ofile))

# Save DSMACC python library
if file_idx[0] >= 0: os.system('cp dsmacc.pyc ./save/results/%s.pyc'%(sys.argv[1]))
print sys.argv[1]+" saved."
