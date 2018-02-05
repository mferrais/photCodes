import numpy as np, glob, os, errno
from shutil import copyfile

"""
Sort photometry output files by filter in dedicated folder.
"""

#--- set things here ----------------------------------------------------------
year = 2017

field = 348

fieldExt = 'raccBVRI_obs1'
fieldExt = '_20171213'
fieldExt = '_20171214'
fieldExt = '_20171215'

night = '20171210'
night = '20171209'
night = '20171213'
night = '20171214'
night = '20171215'

telescope = 'TS'
telescope = 'TN'

pattern = str(field) + '_' + str(year) + '-' + "*"
if telescope == 'TS':
    folder = "E:/trappist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/trappist/pho/"   # TS folder
if telescope == 'TN':
    folder = "E:/troppist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/troppist/pho/"   # TN folder
#-----------------------------------------------------------------------------
folder += str(field) + fieldExt + '/' + night

filterIDtoName = {'12.0':'B','13.0':'V','14.0':'R','15.0':'I'}
filterIDbool = {'12.0':False,'13.0':False,'14.0':False,'15.0':False}

fnames = sorted(glob.glob(folder + '/' + pattern))
nIm = len(fnames)
if nIm == 0:
    raise ValueError('No files found ! Check files path: \n' + \
                     folder + '\n and pattern: \n' + pattern)
print(str(field) + ' : field')
print(night + ' : night')
print(str(nIm) + ' : nIm')

def checkFiltFolder(filtFolder):
    if not os.path.exists(filtFolder):
        try:
            os.makedirs(filtFolder)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

i = 1
for fname in fnames:
    try:
        IDs = np.loadtxt(fname)[:,-1]
    except ValueError:
        print('INDEF found')
    filt = filterIDtoName[str(IDs[0])]
    filtFolder = folder + filt
    if not filterIDbool[str(IDs[0])]:
        checkFiltFolder(filtFolder)
        filterIDbool[str(IDs[0])] = True
    copyfile(fname,folder+filt+'/'+fname[-23:])
    if i%10 == 0:
        print(str(i) + ' /  '+ str(nIm))
    i += 1
print('DONE')
