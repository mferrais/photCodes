import numpy as np, glob, os, errno
from shutil import copyfile

#--- set things here ----------------------------------------------------------
year = 2017

field = 348
field = '348_20171214'
field = '348_20171213'
field = '348_20171215'
field = '348_20171219'

night = '20171214B'
night = 'R56'
night = '20171214V'
night = '20171213R'
night = '20171215B'
night = '20171219R'

filt = 'V'
filt = 'I'
filt = 'B'
filt = 'R'

telescope = 'TS'
telescope = 'TN'

nbrImField = 62
nbrImField = 60

pattern = str(field[:3]) + '_' + str(year) + '-' + "*"

if telescope == 'TS':
    folder = "E:/trappist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/trappist/pho/"   # TS folder
if telescope == 'TN':
    folder = "E:/troppist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/troppist/pho/"   # TN folder
#-----------------------------------------------------------------------------

folder += str(field) + '/' + night + '/'
destFolder = folder

fnames = sorted(glob.glob(folder + pattern))
nIm = len(fnames)
if nIm == 0:
    raise ValueError('No files found ! Check files path: \n' + folder + '\n and pattern: \n' + pattern)
print(str(field) + ' : field')
print(night + ' : night')
print(str(nIm) + ' : nIm')
r = nIm % nbrImField
if r != 0:
    nbrFolder = int(nIm/nbrImField) + 1
else:
    nbrFolder = int(nIm/nbrImField)
print(str(nbrFolder) + ' : nbrFolder')

def createFolder(nFolder):
    if not os.path.exists(nFolder):
        try:
            os.makedirs(nFolder)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

createFolder(folder+filt)

i = 1
n = 0
oldFolderRacc = ''
newFolderRacc = ''
halfNIF = int(nbrImField / 2)
print r, nIm-(r/2)
for fname in fnames:
    if n % nbrImField == 0:
        newFolder = folder+filt+'/'+filt+str(i)
        createFolder(newFolder)
        if i < nbrFolder:
            oldFolderRacc = newFolderRacc
            newFolderRacc = folder+filt+'/'+filt+filt+str(i)+str(i+1)
            createFolder(newFolderRacc)
        print(i)
        i += 1
    copyfile(fname,newFolder+'/'+fname[-23:])
    if n > halfNIF and n < nIm-(r/2) and n % nbrImField > halfNIF:
        copyfile(fname,newFolderRacc+'/'+fname[-23:])
    elif n > halfNIF and n < nIm-(r/2) and n % nbrImField < halfNIF:
        if i == nbrFolder + 1:
            oldFolderRacc = newFolderRacc
        copyfile(fname,oldFolderRacc+'/'+fname[-23:])
    n += 1
