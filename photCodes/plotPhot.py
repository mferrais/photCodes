import numpy as np, matplotlib.pyplot as plt, rtPhot as rtp, sys

"""
Do differential photometry. Input files have time, airmass, raw mag, raw flux.
Light time corrected time, differential flux and or mag and errors are
computed.
"""

#--- set things here ----------------------------------------------------------
#choice of aperture = 1, 2 or 3
ap = 2
ap = 3
ap = 1

fileOut = False
fileOut = True

showCstars = True
showCstars = False
showCstars10 = True
showCstars10 = False

graph = True
graph = False

field = 341
field = 340
field = 342
field = 344
field = 343
field = 347
field = 346
field = 345
field = 348
fieldExt = '_20171214'
fieldExt = '_20171214/20171214V/V'
fieldExt = '_20171213/20171213R/R'
fieldExt = '_20171215/20171215B/B'
fieldExt = '_20171219/20171219R/R'

telescope = 'TS'
telescope = 'TN'

night = '20171127'
night = '20171118'
night = '20170903R'
night = '4'
night = '20171223'
night = '20171207'
night = '20171202'
night = 'R56'
night = str(sys.argv[1])

nightExt = '_20171213'
nightExt = '_20171215'
nightExt = '_20171219'

#------------------------
# to do : files with those bad ID
# ID of stars that are variable:
if night == '20170709' and field == 342:
    badID = [249,309,310,321,331,363,369,456,639,1141]            # 00476
elif night == '20170711' and field == 342:
    badID = [441,475,590,717,1138,1148,1171,1374,1428,1462,1589,1629,1710,
             1729,1799]                                           # 00476
elif night == '20170618' and field == 342:
    badID = [42,44,149,170,174,177,230,252,461,534,543,666,611,656,705,735,
             745,757,777,830,840,913,945,981,992,1008,1011,1013,1020,1123,
             1145,1197,1200,1336,1361,1383]                       # 00476
elif night == '20170628' and field == 342:
    badID = [577,793]                                             # 00476
elif night == '20170602' and field == 340:
    badID = [3,34,35]                                             # 00596
elif night == '20170603' and field == 340:
     badID = [4,6,26,52,55,56,59,84,91,92,99,111,128,133,141,150] # 00596
elif night == '20170920' and field == 344:
    badID = [47,68]
elif night == '20170921' and field == 345:
    badID = [23,24,26,48]
elif night == '20170930' and field == 345:
    badID = [12,15,17]
elif night == '20170922' and field == 344:
    badID = [48,53,61,66]
elif night == '20170923' and field == 344:
    badID = [1,5,34,39,45,47,73]
elif night == '20171107' and field == 344:
    badID = [12]
elif night == '20171127' and field == 346:
    badID = [12,40,66,83,91,96]
elif night == '20171207' and field == 346:
    badID = np.arange(40)
    badID = []
elif night == '20171110' and field == 347:
    if telescope == 'TS':
         badID = [137,380]
         badID = [94,101,232,282,355,373]
         badID =[94,230,280,318]
    else:
        badID = []
elif night == '20171118' and field == 347:
    if telescope == 'TS':
         badID = [28,83,150,190,193,267,342,353,425]
    elif telescope == 'TN':
        badID = [11,13]
    else:
        badID = []
else:
    badID = []
    #badID = np.concatenate((np.arange(4,8),np.arange(14,18),np.arange(31,39)))
    #badID = np.concatenate((np.arange(4,9) ,np.arange(13,31)))
    #badID = np.arange(13,18)
    print badID
    print('No badID. Check field, telescope and night?')
#------------------------

distFname = 'brol1_' + str(field) + '_' + night
distFname = 'brol1_' + str(field)

#**** telescope characteristics ***
D = 60         # telescope diameter in cm
if telescope == 'TS':
    folder = "E:/trappist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/trappist/pho/"   # TS folder
    h = 2315   # altitude in m
    RON = 12.   # ReadOutNoise in el343
    gain = 1.1 # Gain in electron per ADU
elif telescope == 'TN':    ### to change !!! ###
    folder = "E:/troppist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/troppist/pho/"   # TN folder
    h = 2751
    RON = 12
    gain = 1.1

#-----------------------------------------------------------------------------
# load data from folder+fname as an array. Contain a line for each images with:
# otime, airmass, itime, m1, m2, m3, merr1, merr2, merr3, f1, f2, f3,
# area1, area2, area3, msky
# for target and each comp stars
folder += str(field) + fieldExt + '/' + night +'/'
fname = str(field) + fieldExt + '_' +  night + '.dat'
fname = str(field) + nightExt + '_' +  night + '.dat'

# load data from folder + fname as as array. Contain in each lines: time,
# airmass and mag for target and all comparison stars on the image
data = np.loadtxt(folder + fname)

nCstar = (len(data[0])-6) / 13 # number of comparison stars
print(len(data[0])-6, (len(data[0])-6) / 13., nCstar, 'test nCstar')
time = data[:, 0]
airmass = data[:, 1]
itime = data[:, 2]

astMsky = data[:, 3 + 12*(nCstar+1)] * gain
astArea = data[:, 2 + ap + 9*(nCstar+1)]

# target raw mag:
astRawMag = data[:, 2 + ap]
astMagErr = data[:, 2 + ap + 3*(nCstar+1)]
# target raw flux:
astRawFlux = data[:, 2 + ap + 6*(nCstar+1)]

#-----   Error computation of each measurements   -----------------------------
scintillation = [0.09 * D**(-2/3.) * x**(1.75) * (2*it)**(-1/2.) *
                 np.exp(-h/8000.) for x, it in zip(airmass, itime)]
astNoise = [Ft + (Ft * sc)**2 + npix * (bg+RON**2) for Ft,sc,npix,bg
            in zip(astRawFlux, scintillation, astArea, astMsky)]

cStarsSumFluxArea = rtp.getSumCfluxArea(data, ap, nCstar, badID)
cStarsSumFlux = cStarsSumFluxArea[0] # sum of comp star flux
cStarsSumArea = cStarsSumFluxArea[1] # sum of comp star areas
# mean of all comp star background :
cStarsMsky = rtp.getMeanCstarMsky(data, nCstar, badID, gain)

cStarsNoise = [Fc + (Fc * sc)**2 + ncNpix * (bg+RON**2) for Fc,sc,ncNpix,bg
               in zip(cStarsSumFlux, scintillation, cStarsSumArea, cStarsMsky)]
astFluxErr = [np.sqrt((Nt/(Ft**2)) + (Nc/(Fc**2))) for Nt,Ft,Nc,Fc
              in zip(astNoise, astRawFlux, cStarsNoise, cStarsSumFlux)]

print("%.4f" % np.mean(astFluxErr), ' : mean flux error')

nIm = len(data)
print('nbr of images = ' + str(nIm))
print('nbr of comp stars = ', nCstar, ', nbr of bad stars = ', len(badID))

#-----   Light travel time correction   ---------------------------------------
if len(distFname) == 10+len(night):
    distFolder = folder
    print distFolder
else:
    distFolder = folder[:-len(night)-1]
corrTime = rtp.lightTimeCorrection(distFolder+distFname, time)
'''
u = np.arange(len(time))
plt.plot(u, time, '.', label='time')
plt.plot(u, corrTime, '.', label='corrTime')
plt.legend()
plt.show()
'''
#------ differential mag and fluxes -------------------------------------------
#----- mag -----
# mean of all comp star mag for each images:
zeroMag = 15

meanCmag = rtp.getMeanCmag(data, ap, nCstar, badID)
astDiffMag = astRawMag - meanCmag  # differential mag of target

#offset = zeroMag - np.mean(meanCmag)
#astDiffMag -= offset

#astDiffMag -= np.mean(astDiffMag)

#----- flux -----

astDiffFlux = [float(Ft) / Fc for Ft, Fc in zip(astRawFlux, cStarsSumFlux)]
astDiffFlux /= np.mean(astDiffFlux)           #   normalization to 1

"""
zeroMag = 25
zeroMagOffset = 15
zeroFlux = pow(10, -0.4*(zeroMag-25))

meancStarsSumFlux = cStarsSumFlux / nCstar
astDiffFlux = astRawFlux / meancStarsSumFlux

astDiffMag = -2.5 * np.log10(astDiffFlux) + zeroMag
meanCstarsMags = -2.5 * np.log10(np.mean(meancStarsSumFlux)) + zeroMag
offset = zeroMagOffset - meanCstarsMags
astDiffMag += offset
astDiffFlux = pow(10, -0.4*(astDiffMag-zeroMag))
"""
#astDiffFlux = astDiffMag
#print offset, np.mean(astDiffMag)


"""
astDiffFlux = np.asarray(astDiffFlux)
astDiffFlux /= float(nCstar)
cStarMeanFluxes = cStarsSumFlux / nCstar
cStarMeanFluxAll = np.mean(cStarMeanFluxes)
print cStarMeanFluxAll
f0 = 50000.0
fluxCoeff = f0 / cStarMeanFluxAll
print fluxCoeff
astDiffFlux *= fluxCoeff
"""
"""
print cStarMeanFluxAll * nCstar
astDiffFlux = np.asarray(astDiffFlux)
astDiffFlux *= nCstar
"""
"""
astDiffFlux = np.asarray(astDiffFlux)
cStarMeanSumFlux = np.mean(cStarsSumFlux)
astDiffFlux *= cStarMeanSumFlux
"""
#----- saving on file ---------------------------------------------------------


dataOut = np.dstack((corrTime, airmass, astDiffMag, astMagErr,
                     astDiffFlux, astFluxErr))[0]
fnameOut = str(field) + '_' +  night + '_diff_ap' + str(ap) + '.dat'

# Write corrTime, airmass, diff mag of target and mag error,
# diff flux of targer and error on flux,  on a new file
if fileOut:
    print('* Output file has been saved on folder : ' + folder)
    print('With name: ' + fnameOut)
    np.savetxt(folder + fnameOut, dataOut)

#write file for unmatched night ready to go
dataOut = np.dstack((corrTime, astDiffFlux, astFluxErr))[0]
fnameOut = str(field) + '_' +  night + '_ap' + str(ap) + '.dat'
if fileOut:
    print('* File for unmatched data has been saved on folder : ' + folder)
    print('With name: ' + fnameOut)
    np.savetxt(folder + fnameOut, dataOut, fmt=['%10.6f','%.6e','%.2e'])

# write uncorrected time and raw flux
dataRaw = np.dstack((time, astRawMag))[0]
fnameRaw = str(field) + '_' +  night + '_raw_ap' + str(ap) + '.dat'
if fileOut:
    print('* Raw file has been saved on folder : ' + folder)
    print('With name: ' + fnameRaw)
    np.savetxt(folder + fnameRaw, dataRaw)

#--- Graphs ------------------------------------------------------------------

if showCstars10 :
    rtp.seeTenByTenCstars(data, ap, nCstar, corrTime, meanCmag)
if showCstars :
    rtp.seeCstars(data, ap, nCstar, corrTime, meanCmag)

'''
for i in xrange(nCstar + 1):
    error = data[:, (nCstar+1)*3 + 1 + ap + i*3]
    plt.plot(corrTime, error, '.r', label='error')
    plt.title(str(i+1))
    plt.legend()
    plt.show()
'''
'''
plt.plot(corrTime, airmass, 'k.')
plt.plot(corrTime, data[:, 4 + ap + 0*3]-np.mean(data[:, 4 + ap + 0*3]), 'r.')
plt.show()
'''

if graph:
    plt.plot(time, astFluxErr, 'r.', label='astFluxErr')
    plt.legend()
    plt.show()

    plt.plot(corrTime, astMagErr, 'r.', label='astMagErr')
    plt.legend()
    plt.show()

    plt.plot(corrTime, astRawMag, 'r.', label='astRawMag')
    plt.plot(corrTime, meanCmag, 'b.', label='meanCstarsMag')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.show()
    '''
    plt.plot(corrTime, astDiffMag, 'k.', label='diffMag')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.show()
    '''
    u = np.arange(len(astDiffFlux))+1
    cStarMeanFluxes = cStarsSumFlux / nCstar
    cStarMeanFluxes /= np.mean(cStarMeanFluxes)
    plt.plot(u, astDiffFlux, 'k.', label='astDiffFlux')
    plt.plot(u, cStarMeanFluxes, 'r.')
    plt.xlabel('n')
    plt.ylabel('differential flux')
    plt.legend()
    plt.show()

    plt.plot(corrTime, astDiffFlux, 'k.', label='astDiffFlux')
    plt.errorbar(corrTime, astDiffFlux, yerr=astFluxErr, fmt='k.',
                 ecolor='gray')
    plt.xlabel('time')
    plt.ylabel('differential flux')
    plt.legend()
    plt.show()
