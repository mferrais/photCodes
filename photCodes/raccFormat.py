import numpy as np, rtPhot as rtp
#---------------------------------------
field = 346
field = 345
field = 348

fieldExt = 'n2'
fieldExt = 'n4'
fieldExt = 'n1'
fieldExt = '_20171219'
fieldExt = ''
fieldExt = '_20171214'
fieldExt = '_20171215'
fieldExt = '_20171213'
fieldExt = '_20171218'

telescope = 'TS'
telescope = 'TN'

night = '20171219'
night = '20171207'
night = '20171202'
night = '20171214'
night = '20171215'
night = '20171213'
night = '20171218'

ap = 2
ap = 3
ap = 18
ap = 13

filt = 13
filt = 15
filt = 12
filt = 14

#--------------------------------
#fname = 'h.runout1Crop'
#fname = str(field)+'_'+'racc_ap'+str(ap)+'_'+telescope+obs+'_'+night+'.arc'
fname = str(field)+'_'+telescope+'_'+night+'_'+str(ap)+'_'+str(filt)+'.arc'

#fnameOut = str(field)+'_'+'racc_ap'+str(ap)+'_'+telescope+obs+'_'+night
fnameOut = str(field)+'_'+telescope+'_'+night+'_'+str(ap)+'_'+str(filt)

idObject = 1

fileOut = True

if telescope == 'TS':
    folder = "E:/trappist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/trappist/pho/"   # TS folder
if telescope == 'TN':
    folder = "E:/troppist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/troppist/pho/"   # TN folder
#------------------------------------------------------------------------------
folder += str(field)+fieldExt+'/'
distFname = folder + '/' + night + '/' + 'brol1_' + str(field) + '_' + night
distFname = folder + '/' + 'brol1_' + str(field) + '_' + night
print('folder : '+folder)
print('fname : '+fname)

refStarsData = np.loadtxt(folder + fname)
# get time, mag, magerr, flux, fluxerr, filtID of the object
astData = rtp.getAstDataRacc(refStarsData, field, idObject)


astData[0] += 2450000
astData[0] = rtp.lightTimeCorrection(distFname,astData[0]) # time is LTC

#astDataMag = np.dstack((astData[0],astData[1],astData[2]))[0]
#astDataFlux = np.dstack((astData[0],astData[3],astData[4]))[0]
#astData = np.dstack((astData[0],astData[1],astData[2],astData[3],astData[4]))[0]
astDataMag = np.dstack((astData[0],astData[1],astData[2],astData[5]))[0]
astData[3] /= np.mean(astData[3])
astDataFlux = np.dstack((astData[0],astData[3],astData[4],astData[5]))[0]
print astData.shape

if fileOut:
    np.savetxt(folder + fnameOut+'.mag', astDataMag,
    fmt=['%10.6f','%.6e','%.2e','%10.0f'])
    print('With name: ' + fnameOut+'.mag')
    print('*** MAG file has been saved on folder : ' + folder)

    np.savetxt(folder + fnameOut+'.flux', astDataFlux,
    fmt=['%10.6f','%.6e','%.2e','%10.0f'])
    print('*** FLUX file has been saved on folder : ' + folder)
    print('With name: ' + fnameOut+'.flux')
