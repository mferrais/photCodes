import numpy as np, matplotlib.pyplot as plt, rtPhot as rtp

#---------------------------------------
ap = 3
ap = 2
ap = 1


field = 342
field = 341
field = 340
field = 344
field = 347
field = 343
field = 346
field = 345
field = 348

fieldExt = '_20171214/20171214V/V'
fieldExt = '_20171214/20171214R/R'
fieldExt = '_20171213/20171213R/R'
fieldExt = '_20171215/20171215B/B'
fieldExt = '_20171219/20171219R/R'

obs = 'obs1And2'
obs = 'obs1V'
obs = 'obs2B'
obs = 'obs2V'
obs = 'obs2R'
obs = 'obs2I'
obs = 'obs2'
obs = 'obs3'
obs = '20171214V'
obs = '20171214R'
obs = '20171215B'
obs = '20171219R'

telescope = 'TS'
telescope = 'TN'

fileOut = False
fileOut = True


if field == 340: # 596 Scheila
    nights = ['20170602', '20170603']
    Prot = 15.844
elif field == 342: # 476 Hedwig
    nights = ['20170618', '20170628', '20170709', '20170711']
    Prot = 27.33
elif field == 343: # 476 Hedwig
    nights = ['1','2','2bis','3','4','5','6','7']
    nights = ['1','2','3','4']
    nights = ['20170903I']
    Prot = 2.358
elif field == 345: # 89 Julia
    Prot = 11.387
    if telescope == 'TS' and obs == 'obs1':
        nights = ['20170920', '20170921', '20170922', '20170930']
    elif telescope == 'TN' and obs == 'obs1':
        nights = ['20171127', '20171202']
    elif telescope == 'TN' and obs == 'obs2':
        nights = ['20171223']
    elif telescope == 'TN' and obs == 'obs3':
        nights = ['20171207','20171223']
elif field == 346: # 31 Euprhosine
    if telescope == 'TN' and obs == 'obs1':
        nights = ['20171107']
    elif telescope == 'TN' and obs == 'obs2':
        nights = ['20171127','20171202']
    elif telescope == 'TN' and obs == 'obs1And2':
        nights = ['20171107','20171127','20171202']
    if telescope == 'TN' and obs == 'obs3':
        nights = ['20171207']
    Prot = 5.53
elif field == 344: # 24 Themis
    if telescope == 'TS' and obs == 'obs1':
        nights = ['20170920', '20170921', '20170922', '20170923'] # TS
    elif telescope == 'TS' and obs == 'obs2':
        nights = ['20171107', '20171108']
    elif telescope == 'TN' and obs == 'obs1':
        nights = ['20171108']
    Prot = 8.374
elif field == 347: # 20 Massalia
    if telescope == 'TS' and obs == 'obs1':
         nights = ['20171110']
    elif telescope == 'TS' and obs == 'obs2':
         nights = ['20171118']
    elif telescope == 'TN' and obs == 'obs1':
         nights = ['20171110']
    elif telescope == 'TN' and obs == 'obs2':
         nights = ['20171118']
    Prot = 8.098
elif field == 348: # 3200 phaethon
    if telescope == 'TN' and obs == 'obs2B':
         nights = ['B1','B2','B3','B4','B5','B6']
         nightsRacc = ['B12','B23','B34','B45','B56']
    elif telescope == 'TN' and obs == 'obs2R':
         nights = ['R1','R2','R3','R4','R5']
         nightsRacc = ['R12','R23','R34','R45']
    elif telescope == 'TN' and obs == 'obs2V':
         nights = ['V1','V2','V3','V4','V5']
         nightsRacc = ['V12','V23','V34','V45']
    elif telescope == 'TN' and obs == 'obs2I':
         nights = ['I1','I2','I3','I4','I5']
         nightsRacc = ['I12','I23','I34','I45']
    elif telescope == 'TN' and obs == '20171213R':
         nights = ['R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12','R13','R14']
         nightsRacc = ['RR12','RR23','RR34','RR45','RR56','RR67','RR78','RR89','RR910','RR1011','RR1112','RR1213','RR1314']
    elif telescope == 'TN' and obs == '20171214R':
         nights = ['R1','R2','R3','R4','R5','R6']
         nightsRacc = ['R12','R23','R34','R45','R56']
    elif telescope == 'TN' and obs == '20171214V':
         nights = ['V1','V2','V3','V4','V5','V6','V7','V8']
         nightsRacc = ['V12','V23','V34','V45','V56','V67','V78']
    elif telescope == 'TN' and obs == '20171215I':
         nights = ['I1','I2','I3','I4','I5','I6','I7']
         nightsRacc = ['II12','II23','II34','II45','II56','II67']
    elif telescope == 'TN' and obs == '20171215B':
         nights = ['B1','B2','B3','B4','B5','B6','B7']
         nightsRacc = ['BB12','BB23','BB34','BB45','BB56','BB67']
    elif telescope == 'TN' and obs == '20171219R':
         nights = ['R1','R2','R3','R4','R5']
         nightsRacc = ['RR12','RR23','RR34','RR45']
    Prot = 3.603957


yAxisLim = (0.90, 1.10)
yAxisLim = (0.5, 1.5)
ylim = True
ylim = False

stdID = 1752  # 342 ?
stdID = 165   # 344 TS obs1

binData = True
binData = False
binW = 2       # bin width in minutes

matching = False
matching = True

Correction = True
Correction = False

if telescope == 'TS':
    folder = "E:/trappist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/trappist/pho/"   # TS folder
if telescope == 'TN':
    folder = "E:/troppist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/troppist/pho/"   # TN folder
#------------------------------------------------------------------------------
folder += str(field) + fieldExt
plt.style.use('bmh')


getDataRes = rtp.getData(nights, folder, field, ap)
data, datas = getDataRes[0], getDataRes[1]
print(data.shape)
lList = [datai.shape[0] for datai in datas] # list of night's len
print(lList)

######### CORR ##################
if Correction:
    print('CORRECTION')
    #data[:,4][18:420] -= 0.02
    #data[:,4][:605] += 0.039
    #data[:,4][320:] -= 0.04
    #data[:,4][384:] -= 0.04

#################################

time = data[:, 0]
airmass = data[:, 1]
astDiffMag = data[:, 2]
astMagErr = data[:, 3]
astDiffFlux = data[:, 4]
astFluxErr = data[:, 5]
nIm = len(data)


astDiffFlux /= np.mean(astDiffFlux)
#astDiffFlux -= np.mean(astDiffFlux)


#----- field matching ---------------------------------------------------------

if matching:
    getDataResRacc = rtp.getData(nightsRacc, folder, field, ap)
    datasRacc = getDataResRacc[1]
    lListRacc = [datai.shape[0] for datai in datasRacc] # list of night's len
    print(lListRacc, 'lList datasRacc')
    l0 = 0
    for i in range(len(nights)-1):
        print('---'+str(i+1))
        f1Shift, f2Shift = rtp.bridgeMatching(datas[i],datas[i+1],datasRacc[i])
        l1, l2 = len(datas[i]), len(datas[i+1])
        print('---'+str(f1Shift)+' '+str(f2Shift))
        astDiffFlux[l0+l1:l0+l1+l2] += (-f1Shift + f2Shift)
        datas[i+1][:,4] += (-f1Shift + f2Shift)
        datasRacc[i][:,4]  += -f1Shift
        l0 += l1

    astDiffFlux /= np.mean(astDiffFlux)

#----Phase computation --------------------------------------------------------

t0 = time[0]
Porb_JD = Prot/24.
print(Porb_JD, 'PorbJD')

phase = [((t - t0)/Porb_JD) % 1. for t in time]
oldPhase = phase
oldlList = lList

#------ Data Binning ----------------------------------------------------------

print(lList, 'lList')
if binData:
    oldAstFluxErr = astFluxErr
    airmass = rtp.binning(airmass, time, binW, lList)[0]
    astDiffMag = rtp.binning(astDiffMag, time, binW, lList)[0]
    astMagErr = rtp.binning(astMagErr, time, binW, lList)[0]
    astDiffFlux = rtp.binning(astDiffFlux, time, binW, lList)[0]
    astFluxErr = rtp.binning(astFluxErr, time, binW, lList)[0]
    phase = rtp.binning(phase, time, binW, lList)[0]
    resBin = rtp.binWeighted(time, time, oldAstFluxErr, binW, lList)
    time, lList = resBin[0], resBin[1]

#----- get data from racc script ----------------------------------------------
"""
refStarsData = np.loadtxt(folder + '/h.runout1Crop')
astData = rtp.getAstData(refStarsData, field, 1)

'''
ratioCorr = 0.953844443843
i1, i2 = lList[-1], lList[-2]
astData[3][-i1:] /= np.mean(astData[3][-i1:])
astData[3][-i1:] *= (np.mean(astData[3][-i1-i2:-i1]) * ratioCorr)
#astData[3][-i1:] *= (np.mean(astData[3][-i1-i2:-i1]) * ratioCorr)
'''

'''
n12 = np.loadtxt(folder + '/h.runout1Crop12')
n23 = np.loadtxt(folder + '/h.runout1Crop23')
n34 = np.loadtxt(folder + '/h.runout1Crop34')
ast12 = rtp.getAstData(refStarsData, field, 1)
ast23 = rtp.getAstData(refStarsData, field, 1)
ast12 = rtp.getAstData(refStarsData, field, 1)
'''
refAstMag = astData[1]
refAstFlux = astData[3]
# choose flux or mag:
astPhot = refAstFlux
astPhot = refAstMag

print(len(astPhot))

astPhot /= np.mean(astPhot)
plt.plot(astData[0], astPhot, '.')
plt.show()

t0T = astData[0][0]
phaseT = [((t - t0T)/Porb_JD) % 1. for t in astData[0]]
plt.plot(phaseT, astPhot, '.')
plt.xlim((0,1))
plt.show()
n1 = astPhot[:619] ###############
n2 = astPhot[619:] ###############
mn1, mn2 = np.mean(n1), np.mean(n2)
print(mn1, mn2, mn2 / mn1)

astData[3] /= np.mean(astPhot[:lList[0]])# normalization to first night mean


#----- plot std star LC, phase curve of phot and racc data --------------------
# std star check
stdStarData = rtp.getAstData(refStarsData, field, stdID)
stdStarFlux = stdStarData[3]
if len(stdStarFlux) ==0:
    print('Not a valid std ID - std ID not found')
else :
    stdStarFlux /= np.mean(stdStarFlux)
    stdStarTime = stdStarData[0]
    print(len(stdStarFlux), 'len std flux')
    plt.plot(stdStarTime, stdStarFlux, '.', label='stdStarFlux')
    plt.ylim((0.95,1.05))
    plt.legend()
    plt.show()


# Plot flux or mag vs rotational phase  RACC
l0 = 0
for i in range(len(nights)):
    li = oldlList[i]
    plt.plot(oldPhase[l0:li+l0], astPhot[l0:li+l0], '.', label=nights[i])
    l0 += li
plt.xlabel('Rotational Phase')
plt.ylabel('Differential Flux')
plt.title('phase LC of adjusted data RACC')
plt.legend()
plt.xlim((0,1))
if ylim:
    plt.ylim(yAxisLim)
plt.show()

# Plot flux or mag vs rotational phase  PHOT
l0 = 0
for i in range(len(nights)):
    li = lList[i]
    #plt.plot(phase[l0:li+l0], astDiffMag[l0:li+l0], '.', label=nights[i]) # mag
    plt.plot(phase[l0:li+l0], astDiffFlux[l0:li+l0], '.', label=nights[i])  # flux
    l0 += li
#plt.gca().invert_yaxis() # invert mag axis
plt.xlabel('Rotational Phase')
#plt.ylabel('$\Delta$ Rmag')
plt.ylabel('Differential Flux')
plt.title('phase LC of non adjusted data PHOT')
plt.legend()
plt.xlim((0,1))
if ylim:
    plt.ylim(yAxisLim)
plt.show()

# ----- Night shifting --------------------------------------------------------
nightsRacc = []
nightsPhot = []
l0 = 0
oldl0 = 0
for i in range(len(nights)):
    li = lList[i]
    oldli = oldlList[i]
    nightsRacc.append(astPhot[oldl0:oldl0+oldli])
    nightsPhot.append(astDiffFlux[l0:l0+li])
    l0 += li
    oldl0 += oldli

meanNightsRacc = [np.mean(night) for night in nightsRacc]
meanNightsPhot = [np.mean(night) for night in nightsPhot]
ratiosRacc = [mnR / meanNightsRacc[0] for mnR in meanNightsRacc]
print(ratiosRacc, 'ratiosRacc')

colors = ['C0.','C1.','C2.','C3.','C4.','C5.','C6.','C7.']

l0 = 0
for i in range(len(nights)):
    li = oldlList[i]
    u = np.arange(l0,l0+li)
    plt.plot(u, nightsRacc[i], colors[i], label=nights[i])
    plt.plot(u, [meanNightsRacc[i] for e in u], colors[i]+'-')
    l0 += li
plt.title('racc')
plt.legend()
if ylim:
    plt.ylim(yAxisLim)
plt.show()

l0 = 0
for i in range(len(nights)):
    li = lList[i]
    u = np.arange(l0,l0+li)
    plt.plot(u, nightsPhot[i], colors[i], label=nights[i])
    plt.plot(u, [meanNightsPhot[i] for e in u], colors[i]+'-')
    l0 += li
plt.title('before shift')
plt.legend()
if ylim:
    plt.ylim(yAxisLim)
plt.show()

# shift all night except first one
l0 = lList[0]
for i in range(len(nights)-1):
    li = lList[i+1]
    data[:,4][l0:l0+li] *= ratiosRacc[i+1]
    l0 += li

# update of phot nights mean and nightsPhot
l0 = 0
meanNightsPhot = []
for i in range(len(nights)):
    li = lList[i]
    meanNightsPhot.append(np.mean(data[:,4][l0:l0+li]))
    nightsPhot.append(astDiffFlux[l0:l0+li])
    l0 += li

l0 = 0
for i in range(len(nights)):
    li = lList[i]
    u = np.arange(l0,l0+li)
    plt.plot(u, nightsPhot[i], colors[i], label=nights[i])
    plt.plot(u, [meanNightsPhot[i] for e in u], colors[i]+'-')
    l0 += li
plt.title('after shift')
plt.legend()
if ylim:
    plt.ylim(yAxisLim)
plt.show()

ratiosPhot = [mnP / meanNightsPhot[0] for mnP in meanNightsPhot]
print(ratiosPhot, 'ratiosPhot')

# renormalization of the whole curve to 1
data[:,4] /= np.mean(data[:,4])
"""
#-----    saving on file:   ---------------------------------------------------
if fileOut and binData:
    print('Warning : data are binned, files not saved !')
fnameOut = '/' + str(field) + '_diff_ap' + str(ap) +'_' +telescope+obs + '.dat'
if fileOut and not binData:
    print('* Full output file has been saved on folder : ' + folder)
    print('With name: ' + fnameOut[1:-4]+'_full.dat')
    #np.savetxt(folder + fnameOut, data)
    np.savetxt(folder + fnameOut[:-4]+'_full.dat', data,
               fmt=['%10.6f','%10.4f','%.6e','%.2e','%.6e','%.2e'])
    print('* Output file has been saved on folder : ' + folder)
    print('With name: ' + fnameOut[1:])
    dataFlux = np.dstack((data[:,0],data[:,4],data[:,5]))[0]
    np.savetxt(folder + fnameOut, dataFlux,
               fmt=['%10.6f','%.6e','%.2e'])

# raw file:
fnameRaw = '/' + str(field) + '_raw_ap' + str(ap) + '_'+telescope+ obs + '.dat'
dataRaw = rtp.getRawData(nights, folder, field, ap)
if fileOut and not binData:
    print('* Raw file has been saved on folder : ' + folder)
    print('With name: ' + fnameRaw[1:])
    np.savetxt(folder + fnameRaw, dataRaw)
#-----  Graphs  ---------------------------------------------------------------

# Plot flux or mag in a sequential way
l0 = 0
for i in range(len(nights)):
    li = lList[i]
    plt.plot(np.arange(li)+l0, astDiffFlux[l0:li+l0], '.', label=nights[i])
    l0 += li
plt.xlabel('n')
plt.ylabel('Differential Flux')
plt.title('Sequential LC')
plt.legend()
plt.show()

# Plot flux error in a sequential way
l0 = 0
for i in range(len(nights)):
    li = lList[i]
    plt.plot(np.arange(li)+l0, astFluxErr[l0:li+l0], '.', label=nights[i])
    l0 += li
plt.xlabel('n')
plt.ylabel('error')
plt.title('Error on flux')
plt.legend()
plt.show()


# Plot flux or mag vs time in JD
l0 = 0
for i in range(len(nights)):
    li = lList[i]
    plt.plot(time[l0:li+l0], astDiffFlux[l0:li+l0], '.', label=nights[i])
    l0 += li
if matching:
    for i in range(len(nightsRacc)):
        plt.plot(datasRacc[i][:,0], datasRacc[i][:,4], 'k,', label=nightsRacc[i])
plt.xlabel('JD time')
plt.ylabel('Differential Flux')
plt.title('JD time LC')
plt.legend()
plt.show()

# Plot flux or mag vs rotational phase
l0 = 0
for i in range(len(nights)):
    li = lList[i]
    #plt.plot(phase[l0:li+l0], astDiffMag[l0:li+l0], '.', label=nights[i]) # mag
    plt.plot(phase[l0:li+l0], astDiffFlux[l0:li+l0], '.', label=nights[i])  # flux
    l0 += li
#plt.gca().invert_yaxis() # invert mag axis
plt.xlabel('Rotational Phase')
#plt.ylabel('$\Delta$ Rmag')
plt.ylabel('Differential Flux')
plt.title('phase LC of shifted data')
plt.legend()
plt.xlim((0,1))
if ylim:
    plt.ylim(yAxisLim)
plt.show()
