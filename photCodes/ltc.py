import numpy as np, matplotlib.pyplot as plt, rtPhot as rtp

"""
Compute and apply light time correction to individual files
"""

ap = 3
ap = 2
ap = 1

fileOut = False
fileOut = True

telescope = 'TS'

field = 343

night = '20170903'

distFname = 'brol1_' + str(field)
distFname = 'brol1_' + str(field) + '_' + night

if telescope == 'TS':
    folder = "E:/trappist/pho/"                    # win TS folder
    folder = "/media/marin/TRAPPIST3/trappist/pho/"   # TS folder
if telescope == 'TN':
    folder = "E:/troppist/pho/"                       #clean win TS folder
    folder = "/media/marin/TRAPPIST3/troppist/pho/"   # TN folder
folder += str(field) +'/'

fname = 'AIJ/03122AIJ.txtcorr'
data = np.loadtxt(folder + fname)
time = data[:,0]
astDiffFlux = data[:,1]
astFluxErr = data[:,2]

#astDiffFlux /= np.mean(astDiffFlux)


#-----   Light travel time correction   ---------------------------------------
c_au = (299792458. * 24*3600)/149597870700. # speed of light [au/JD]
if len(distFname) == 18:
    distFolder = folder
else:
    distFolder = folder[:-2]

distFolder = '/media/marin/TRAPPIST3/trappist/pho/343/' + night +'/'

brol1 = np.loadtxt(distFolder + distFname)
topoRange = brol1[:,3] #target topocentric range, [au]
topoTime = brol1[:,0] # time for topo range
# light travel time correction, in JD :
lightTimeCorrVal = [d/c_au for d in topoRange]
lightTimeCorrFull = []
for t in time:
    deltaTime = [abs(tt - t) for tt in topoTime]
    minDt = min(deltaTime)
        if minDt > 1 and abs(minDt-2450000) > 1:
            raise ValueError('time diff too big : '+str(minDt))
    index = deltaTime.index(minDt)
    lightTimeCorrFull.append(lightTimeCorrVal[index])
corrTime = time - lightTimeCorrFull # corrected time
print(time[0],corrTime[0])

u = np.arange(len(time))
plt.plot(u, time, '.', label='time')
plt.plot(u, corrTime, '.', label='corrTime')
plt.legend()
plt.show()

#---- saving on file ---------------------------------------------------------
if fileOut:
    print('* File has been saved on folder : ' + folder)
    print('With name: ' + fname)
    data = np.dstack((corrTime, astDiffFlux, astFluxErr))[0]
    np.savetxt(folder + fname, data, fmt=['%10.6f','%.6e','%.2e'])
