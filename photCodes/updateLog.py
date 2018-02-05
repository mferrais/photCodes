import glob
from astropy.io import fits


"""
Update WCS Status on TRAP.log so those images who had their WCS corrected
can be used in photometry.
Remember to do ! sh TRAPall3.sh  after to update TRAPall3.log
"""

tele = 'TS'
tele = 'TN'

year = 2017
jd1 = 8085
jd2 = 8091
binn = 1  # image binning

if tele == 'TS':
    folder = "E:/trappist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/trappist/r"   # TS folder
    pattern = 'TRAP.' + str(year) + '-' + '*'
    histInd = 1
if tele == 'TN':
    folder = "E:/troppist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/troppist/r"   # TN folder
    pattern = 'TROP.' + str(year) + '-' + '*'
    histInd = 2

#------------------------------------------------------------------------------

folder += str(year) + '/' + str(jd1)+'_'+str(jd2)+'_B'+str(binn)+'R11'+ '/'
print(folder+pattern, 'fold + patt')
fnames = sorted(glob.glob(folder + pattern))
nIm = len(fnames)
print(str(nIm)+ ' images found')

logFile = folder + 'TRAP.log'
with open(logFile, 'r') as f:
    log = f.readlines()[1:]
print(str(len(log)) + ' lines on log file')
log = [x.strip().split() for x in log]

nbrCorrWCS = 0
nbrNoWCS = 0
print('-----')
for i in xrange(nIm):
    fname = fnames[i]
    hdulist = fits.open(fname)
    header = hdulist[0].header
    try:
        wcsOk = header['HISTORY'][histInd]
        if log[i][-7] == 'NO':
            log[i][-7], log[i][-6] = 'WCS', 'OK'
            nbrCorrWCS += 1
    except IndexError:
        print('NO WCS : ' + fname[54:])
        print('-----')
        nbrNoWCS += 1
    hdulist.close()
    if i%1000 == 0:
        print(str(i)+' images checked over '+str(nIm))

log = [' '.join(line) for line in log]
print(len(log))
print('No WCS to Ok : ' + str(nbrCorrWCS))
print('No WCS : ' + str(nbrNoWCS))

newLogName = folder+'TRAP.logTmp'
newLog = open(newLogName, 'w')
textList = map(lambda x: x+"\n", log)
textList.insert(0, '\n')
newLog.writelines(textList)
newLog.close()
print('updated log saved as TRAP.logTMP')
