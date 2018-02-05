import numpy as np, glob, os.path as op

#---------------------------------------
year = 2017

field = 341
field = 342
field = 340
field = 344
field = 347
field = 343
field = 346
field = 345
field = 348

folder = "E:/trappist/pho/"                       # win TS folder
folder = "/media/marin/TRAPPIST3/trappist/pho/"   # TS folder
folder = "/media/marin/TRAPPIST3/troppist/pho/"   # TN folder
#---------------------------------------
pattern = '/' + str(field) + '_' + str(year) + '-' + "*"
folder += str(field)
destFolder = folder + 'racc/'
print(folder)
print(destFolder)
if not op.isdir(destFolder):
    raise ValueError('/'+str(field)+'racc' + ' folder does not exist')

formatData = ['%10.3f','%10.3f','%10.0f','%10.0f','%10.3f','%10.3f','%10.3f',
              '%10.3f','%10.3f','%10.3f','%10.10f','%10.6f','%10.0f']

fnames = sorted(glob.glob(folder + pattern))
nIm = len(fnames)
if nIm == 0:
    raise ValueError('No files found ! Check files path: \n' + \
                     folder + '\n and pattern: \n' + pattern)

i = 1
for fname in fnames:
    data = np.loadtxt(fname)
    data = np.hstack((data[:,0:10], data[:,-3:]))
    np.savetxt(destFolder + fname[len(fname)-23:], data, fmt=formatData )
    if i%10 == 0:
        print(str(i) + ' /  '+ str(nIm))
    i += 1
print('DONE')
