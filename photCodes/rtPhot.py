import numpy as np, matplotlib.pyplot as plt


#------ formatPhot ------------------------------------------------------------
def inAllIm(listID):
    """search the stars detected on all images, via their IDs"""
    inAll = listID[0] # we start with the IDs of the first image
    newInAll = []
    for IDs in listID:
        for ID in IDs:
            if ID in inAll:
                newInAll.append(ID)
        inAll = newInAll
        newInAll = []
    return inAll

#------ plotPhot --------------------------------------------------------------

def getMeanCmag(data, ap, nCstar, badID):
    """ Compute for each images the mean magnitude off all the comparison
        stars, except for those whose ID is in the list badID (variable stars)
    """
    meanCmag = []
    for line in data:
        magLine, index = [], 1
        for i in xrange(5 + ap, 5 + ap + 3*nCstar, 3):
            if index not in badID:
                magLine.append(line[i])
            index += 1
        meanCmag.append(np.mean(magLine))
    return np.array(meanCmag)

def getSumCfluxArea(data, ap, nCstar, badID):
    """ Compute for each images the flux sum of all the comparison
        stars, except for those whose ID is in the list badID (variable stars)
    """
    sumCflux = []
    sumCarea = []
    for line in data:
        fluxLine, areaLine, index = [], [], 1
        i1 = 5 + ap + 6*(nCstar+1)
        for i in xrange(i1, i1 + 3*nCstar, 3):
            if index not in badID:
                fluxLine.append(line[i])
                areaLine.append(line[i + 3*(nCstar+1)])
            index += 1
        sumCflux.append(np.sum(fluxLine))
        sumCarea.append(np.sum(areaLine))
    return np.array(sumCflux), np.array(sumCarea)

def getMeanCstarMsky(data, nCstar, badID, gain):
    """ Compute for each images the mean background level for all comp stars
    """
    meanMsky = []
    lenLine = len(data[0])
    for line in data:
        mskyLine, index = [], 1
        for i in xrange(4 + 12*(nCstar+1), lenLine):
            if index not in badID:
                mskyLine.append(line[i]*gain)
            index += 1
        meanMsky.append(np.mean(mskyLine))
    return np.array(meanMsky)


def seeCstars(data, ap, nCstar, corrTime, meanCmag):
    """ Plot raw and diff mag of each comp stars one by one"""
    for i in xrange(nCstar):
        cStarRawMag = data[:, 5 + ap + i*3]
        cStarDiffMag = cStarRawMag - meanCmag
        plt.plot(corrTime, cStarRawMag - np.mean(cStarRawMag), '.r', label='rawMag')
        plt.plot(corrTime, cStarDiffMag - np.mean(cStarDiffMag) + 0.2,
        '.k', label='diffMag')
        plt.title(str(i+1))
        plt.gca().invert_yaxis()
        plt.legend()
        plt.show()

def seeTenByTenCstars(data, ap, nCstar, corrTime, meanCmag):
    """ Plot raw and diff mag of each comp stars 10 by 10"""
    for i in xrange(0, nCstar, 10):
        for j in xrange(10):
            cStarRawMag = data[:, 5 + ap + (i+j)*3]
            cStarDiffMag = cStarRawMag - meanCmag
            plt.plot(corrTime, cStarRawMag - np.mean(cStarRawMag) + j, '.r')
            plt.plot(corrTime, cStarDiffMag - np.mean(cStarDiffMag) + j + 0.5,
            '.k')
        plt.title(str(i))
        plt.show()

#----- concatPhot + graphPhot -------------------------------------------------

def binning(listToBin, time, binT, nightLen):
    """ bin listToBin binT by binT """
    tmpList, listBinned, newNightLen = [], [], []
    tmpt0 = time[0]
    for i in xrange(len(listToBin)):
        c1, c3 = len(tmpList) == 0, i not in nightLen
        if (c1 or abs(tmpt0-time[i]) < binT/1440.) and c3:
            tmpList.append(listToBin[i])
        else:
            listBinned.append(np.mean(tmpList))
            tmpList = [listToBin[i]]
            tmpt0 = time[i]
        if i in nightLen:
            newNightLen.append(len(listBinned)-sum(newNightLen))
    newNightLen.append(len(listBinned)-sum(newNightLen))
    return listBinned, newNightLen

def binWeighted(listToBin, time, fluxErr, binT, nightLen):
    """ Bin datas in listToBin together in a time span given by binT,
        separated in several nights given by nightLen and using a
        weighted average using the flux errors
    """
    tmpList, tmpErr, listBinned, newNightLen = [], [], [], []
    tmpt0 = time[0]
    for i in xrange(1,len(listToBin)):
        c1, c3 = len(tmpList) == 0, i not in nightLen
        if (c1 or abs(tmpt0-time[i]) < binT/1440.) and c3:
            tmpList.append(listToBin[i])
            tmpErr.append(fluxErr[i])
        else:
            #listBinned.append(np.average(tmpList, weights=tmpErr))
            tmpErr = [1./x**2 for x in tmpErr]
            listBinned.append(np.average(tmpList, weights=tmpErr))
            tmpList = [listToBin[i]]
            tmpt0 = time[i]
            tmpErr = [fluxErr[i]]
        if i in nightLen:
            newNightLen.append(len(listBinned)-sum(newNightLen))
    newNightLen.append(len(listBinned)-sum(newNightLen))
    return np.asarray(listBinned), newNightLen

#----- concatPhot -------------------------------------------------------------

def getRawData(nightL, folder, field, ap):
    ''' stack togethere  raw datas from different night in the array data
    '''
    folderL = [folder + '/' + night + '/' for night in nightL]
    fnameL = [str(field) + '_' +  night + '_raw_ap' + str(ap) + '.dat'
              for night in nightL]
    datas = [np.loadtxt(folder+fname) for folder,fname in zip(folderL,fnameL)]
    data = np.vstack(datas)
    return data

def getData(nightL, folder, field, ap):
    ''' stack togethere datas from different night in the array data
        In datas the nights are kept separated
    '''
    folderL = [folder + '/' + night + '/' for night in nightL]
    fnameL = [str(field) + '_' +  night + '_diff_ap' + str(ap) + '.dat'
              for night in nightL]
    datas = [np.loadtxt(folder+fname) for folder,fname in zip(folderL,fnameL)]
    """ #normalization
    for i in range(len(nights)):
        datas[i][:,2] /= np.mean(datas[i][:,2])
    """
    data = np.vstack(datas)
    return data, datas

def bridgeMatching(f1,f2,fracc):
    """Match 2 consecutive field of the same obs night using a third field
       overlaping both
    """
    f1Time, f2Time, fraccTime = f1[:,0], f2[:,0], fracc[:,0]
    f1Flux, f2Flux, fraccFlux = f1[:,4], f2[:,4], fracc[:,4]
    l1 = len(f1Flux)
    jd1racc, jd2racc = fraccTime[0], fraccTime[-1]
    jdf1, n1 = f1Time[0], 0
    while jdf1 != jd1racc:
        n1 += 1
        jdf1 = f1Time[n1]
    jdf2, n2 = f2Time[0], 0
    while jdf2 != jd2racc:
        n2 += 1
        jdf2 = f2Time[n2]

    f1Fluxcut, fraccFluxcut1 = f1Flux[n1:],fraccFlux[0:l1-n1]
    meanf1Flux = np.mean(f1Flux)
    meanf2Flux = np.mean(f2Flux)
    f1fraccDist = []
    d = 0.02
    while len(f1fraccDist) < 5 and d <= 1:
        f1fraccDist = [yracc-y1 for y1,yracc in zip(f1Fluxcut, fraccFluxcut1)
                       if abs(y1-meanf1Flux)<d]
        d += 0.02
    print('d1='+str(d-0.02))
    f1Shift = np.mean(f1fraccDist[1:])
    print(str(len(f1fraccDist))+' : len f1fraccDist')

    f2Fluxcut, fraccFluxcut2 = f2Flux[:n2], fraccFlux[l1-n1:]
    f2fraccDist = []
    d = 0.02
    while len(f2fraccDist) < 5 and d <= 1:
        f2fraccDist = [yracc-y2 for y2,yracc in zip(f2Fluxcut, fraccFluxcut2)
                       if abs(y2-meanf2Flux)<d]
        d += 0.02
    print('d2='+str(d-0.02))
    f2Shift = np.mean(f2fraccDist[1:])
    print(str(len(f2fraccDist))+' : len f2fraccDist')
    return f1Shift, f2Shift

def meanShift(f1, f2, time1, time2):
    """ Compute mean shift between f1 and f2
    """
    shifts, m1, m2 = [], np.mean(f1), np.mean(f2)
    sign, i = 1, 0
    if len(f1) > len(f2): # we want len f1 < f2, same for time1 and time2
        f1, f2 = f2, f1
        time1, time2 = time2, time1
        sign = -1
    while i < len(f1):
        t1 = time1[i]
        deltaTime = [abs(t2-t1) for t2 in time2]
        i2 = deltaTime.index(min(deltaTime))
        if abs(f1[i]-m1) < 0.1 and abs(f2[i2]-m2) < 0.1:
            shifts.append(f1[i] - f2[i2])
        i += 1
    return np.mean(shifts)


# --- if racc Data used ---

def getAstData(refData, field, ID):
    ''' Get adjusted data from the raccord script
        refData :
        id = id number from iraf scripts - (field * 1e5)
    '''
    idAll = refData[:,0] - (field * 1e5)
    timeAll = refData[:,1]
    magAll = refData[:,2]
    errAll = refData[:,3]
    time = []
    mag = []
    merr = []
    for i in xrange(len(idAll)):
        if idAll[i] ==  float(ID):
            time.append(timeAll[i])
            mag.append(magAll[i])
            merr.append(errAll[i])
    time, mag, merr = np.array(time), np.array(mag), np.array(merr)
    flux, ferr = pow(10, -mag/2.5), 10*merr
    return np.array((time, mag, merr, flux, ferr))

def getAstDataRacc(refData, field, ID):
    ''' Get adjusted data from the raccord script
        refData :
        id = id number from iraf scripts - (field * 1e5)
    '''
    idFiltAll = refData[:,0]
    idAll = refData[:,1] - (field * 1e5)
    timeAll = refData[:,2]
    magAll = refData[:,4]
    errAll = refData[:,5]
    time = []
    mag = []
    merr = []
    filtId = []
    for i in xrange(len(idAll)):
        if idAll[i] ==  float(ID):
            time.append(timeAll[i])
            mag.append(magAll[i])
            merr.append(errAll[i])
            filtId.append(idFiltAll[i])
    time, mag, merr, filtId = np.array(time), np.array(mag), \
                              np.array(merr), np.array(filtId)
    flux, ferr = pow(10, -mag/2.5), 10*merr
    return np.array((time, mag, merr, flux, ferr, filtId))

#---- multiple places ---------------------------------------------------------

def lightTimeCorrection(fname, time):
    "Compute light time correction using asteroid topocentric distance"
    c_au = (299792458. * 24*3600)/149597870700. # speed of light [au/JD]
    brol1 = np.loadtxt(fname)
    topoRange = brol1[:,3] #target topocentric range, [au]
    topoTime = brol1[:,0] # time for topo range
    # light travel time correction, in JD :
    lightTimeCorrVal = [d/c_au for d in topoRange]
    lightTimeCorrFull = []
    for t in time:
        deltaTime = [abs(tt - t) for tt in topoTime]
        if deltaTime[0] > 1:
            deltaTime = [abs(x-2450000) for x in deltaTime]
        minDt = min(deltaTime)
        if minDt > 0.1:
            print deltaTime
            raise ValueError('time diff too big : '+str(minDt) +\
                             ', t: '+ str(t) + ', file: ' + fname )
        index = deltaTime.index(minDt)
        lightTimeCorrFull.append(lightTimeCorrVal[index])
    ltcTime = time - lightTimeCorrFull # corrected time
    return ltcTime
