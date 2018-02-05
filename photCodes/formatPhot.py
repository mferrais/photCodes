import numpy as np, glob, os, sys
import matplotlib.pyplot as plt
import rtPhot as rtp

"""
Get together datas generated by the IRAF phot functions in
a unique file, and select only the star present on all images
"""

#--- set things here ----------------------------------------------------------
year = 2017

field = 341
field = 342
field = 340
field = 344
field = 347
field = 350
field = 343
field = 346
field = 345
field = 348
fieldExt = '_20171214/20171214V/V'
fieldExt = '_20171213/20171213R/R'
fieldExt = '_20171215/20171215B/B'
fieldExt = '_20171219/20171219R/R'

night = '20170621'
night = '20170618'
night = '20170628'
night = '20170709'
night = '20170711'
night = '20170603'
night = '20170602'
night = '20170930'
night = '20170921'
night = '20170922'
night = '20170923'
night = '20171108'
night = '20171110'
night = '20170920'
night = '20171118'
night = '20171107'
night = '20171202'
night = '20171127'
night = '20170903'
night = '20170903I'
night = '20171209B'
night = 'I45'
night = '20171223'
night = '20171207'
night = '20171218'
night = '20171209'
night = '20171214B'
night = 'R56'
night = 'V1'
night = str(sys.argv[1])

nightExt = '_20171215'
nightExt = ''
nightExt = '_20171219'

telescope = 'TS'
telescope = 'TN'

# only for racc file: h.runout1Crop
if night == '20170603' and field == 340:
    lastj = 314610
else :
    lastj = 0

pattern = str(field) + '_' + str(year) + '-' + "*"
if telescope == 'TS':
    folder = "E:/trappist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/trappist/pho/"   # TS folder
if telescope == 'TN':
    folder = "E:/troppist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/troppist/pho/"   # TN folder
#-----------------------------------------------------------------------------

folder += str(field) + fieldExt + '/' + night + '/'
destFolder = folder
outputFname = str(field) + nightExt + '_' +  night + '.dat'

fnames = sorted(glob.glob(folder + pattern))
nIm = len(fnames)
if nIm == 0:
    raise ValueError('No files found ! Check files path: \n' + folder + '\n and pattern: \n' + pattern)
print(str(field) + ' : field')
print(night + ' : night')
print(str(nIm) + ' : nIm')

time = []
airmass = []
listID = []
itime = []
msky = []
aperture1 = []
aperture2 = []
aperture3 = []
errAp1 = []
errAp2 = []
errAp3 = []
flux1 = []
flux2 = []
flux3 = []
area1 = []
area2 = []
area3 = []

# Fill list above with datas from all images for the night
badIm = 0
for fname in fnames:
    #print(fname)
    try:
        imData = np.loadtxt(fname)
    except ValueError:
        print('INDEF found')
    listID.append(imData[:,2])
    if imData[:,2][0] != 1.0:
        badIm += 1
        print(fname)
        print(imData[:,2][0])
        dfname = fname[0:-24-len(night)]  + fname[-23:]
        os.rename(fname, dfname)
    airmass.append(imData[0][-2])
    time.append(imData[0][-3])
    itime.append(imData[:,-4])
    msky.append(imData[:,-8])
    aperture1.append(imData[:,4])
    errAp1.append(imData[:,5])
    aperture2.append(imData[:,6])
    errAp2.append(imData[:,7])
    aperture3.append(imData[:,8])
    errAp3.append(imData[:,9])
    flux1.append(imData[:,10])
    flux2.append(imData[:,11])
    flux3.append(imData[:,12])
    area1.append(imData[:,13])
    area2.append(imData[:,14])
    area3.append(imData[:,15])
print('--')
print(str(badIm) + ' : badIm')

compStarsID = rtp.inAllIm(listID) # list of comp star ID + target ID
print(str(len(compStarsID)-1) + ' : len compStarsID')

#----- if using racc data -----
"""
raccName = 'h.runout1Crop'
raccName = '346_TN_20171207_2_14.arc'

'''
def refStarsID(refData):
    time = refData[:,3]
    print time.shape
    allID = refData[:,2]
    t = time[0]
    IDs = []
    j = 0
    for i in xrange(690,870):
        ID = []
        while (abs(time[j] - t) < 0.00001):
            ID.append(allID[j])
            j += 1
        IDs.append(ID)
        if 34000030.0 in ID : print 'yess'
        else : print 'nooo'
        t = time[j]
    return IDs
'''

def refStarsIDALLold(refData):
    time = refData[:,1]
    print(str(time.shape) + '  : time.shape')
    allID = refData[:,0] - (field*1e5)
    t = time[lastj]
    IDs = []
    j = lastj
    for i in xrange(nIm):
        ID = []
        while (j < nbrRefLine) and (abs(time[j] - t) < 0.00001):
            ID.append(allID[j])
            j += 1
        IDs.append(ID)
        if j < nbrRefLine:
            t = time[j]
    print(str(j) + ' : j')
    print(str(len(IDs)) + 'len IDs')
    return IDs


def refStarsIDALL(refData):
    time = refData[:,2]
    print(str(time.shape) + '  : time.shape')
    allID = refData[:,1] - (field*1e5)
    print allID[0], 'allid0'
    t = time[lastj]
    IDs = []
    j = lastj
    for i in xrange(nIm):
        ID = []
        while (j < nbrRefLine) and (abs(time[j] - t) < 0.00001):
            ID.append(allID[j])
            j += 1
        IDs.append(ID)
        if j < nbrRefLine:
            t = time[j]
    print(str(j) + ' : j')
    print(str(len(IDs)) + 'len IDs')
    return IDs


'''
refStarsData = np.loadtxt(folder[:-9]+'h.out1cro')
print refStarsData.shape
tmpRefIDlist = refStarsID(refStarsData)
'''
#refStarsData = np.loadtxt(folder[:-9]+'raccName)
refStarsData = np.loadtxt(folder+raccName)
nbrRefLine = refStarsData.shape[0]
tmpRefIDlist = refStarsIDALL(refStarsData)
refStarsID = rtp.inAllIm(tmpRefIDlist)


print('----------------------------------------------------------')
print(refStarsID)
print('----------------------------------------------------------')
print(compStarsID)
print('----------------------------------------------------------')

refAndCompStarsID = rtp.inAllIm([refStarsID,compStarsID])
#print(refAndCompStarsID)
print(str(len(refAndCompStarsID)) + ' : len of ref and compStarsID')

compStarsID = [1] + refAndCompStarsID # only IDs ok for racc and phot are kept

print compStarsID
"""


#  -----------------------------
nStars = len(compStarsID)
print(str(nStars-1) + ' : number of stars present on all images')
print('--')
#----- get phot data ----------------------------------------------------------
# create an array of nIm lines, with on each lines :
# time, airmass, expTime, mag, merr, flux, area, msky
# for aperture 1,2,3 and for the target and then all comparison stars
newData = np.empty((nIm, 3 + 12*nStars + nStars))
print(str(newData.shape) + ' : newData shape')
for i in xrange(nIm):
    newData[i][0], newData[i][1] = time[i], airmass[i]
    newData[i][2] = itime[i][0]
    IDs = listID[i]
    index = 0
    for j in xrange(len(IDs)): # for each line (each stars detected in image)
        if IDs[j] in compStarsID:
            newData[i][3+3*index] = aperture1[i][j]
            newData[i][4+3*index] = aperture2[i][j]
            newData[i][5+3*index] = aperture3[i][j]
            newData[i][3+3*(index+nStars)] = errAp1[i][j]
            newData[i][4+3*(index+nStars)] = errAp2[i][j]
            newData[i][5+3*(index+nStars)] = errAp3[i][j]
            newData[i][3+3*(index+(2*nStars))] = flux1[i][j]
            newData[i][4+3*(index+(2*nStars))] = flux2[i][j]
            newData[i][5+3*(index+(2*nStars))] = flux3[i][j]
            newData[i][3+3*(index+(3*nStars))] = area1[i][j]
            newData[i][4+3*(index+(3*nStars))] = area2[i][j]
            newData[i][5+3*(index+(3*nStars))] = area3[i][j]
            newData[i][3+index+12*nStars] = msky[i][j]
            index += 1

#--- saving on file -----------------------------------------------------------
print('*** Output file has been saved on folder : ' + destFolder)
print('*** With name: ' + outputFname)
np.savetxt(destFolder + outputFname, newData)
