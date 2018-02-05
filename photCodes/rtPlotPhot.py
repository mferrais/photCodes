#--- Functions ----------------------------------------------------------------

def getMeanCmag(data, nCstar, badID):
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


def getSumCfluxOLDWAY(data, nCstar, badID):
    """ Compute for each images the (flux/expTime) sum off all the comparison
        stars, except for those whose ID is in the list badID (variable stars)
    """
    sumCflux = []
    for line in data:
        fluxLine, index = [], 1
        for i in xrange(5 + ap, 5 + ap + 3*nCstar, 3):
            if index not in badID:
                fluxLine.append(pow(10, 0.4*(25-line[i])))
            index += 1
        sumCflux.append(np.sum(fluxLine))
    return np.array(sumCflux)

def getSumCflux(data, nCstar, badID):
    """ Compute for each images the flux sum of all the comparison
        stars, except for those whose ID is in the list badID (variable stars)
    """
    sumCflux = []
    for line in data:
        fluxLine, index = [], 1
        i1 = 5 + ap + 6*(nCstar+1)
        for i in xrange(i1, i1 + 3*nCstar, 3):
            if index not in badID:
                fluxLine.append(line[i])
            index += 1
        sumCflux.append(np.sum(fluxLine))
    return np.array(sumCflux)

def getSumCfluxArea(data, nCstar, badID):
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

def getMeanCstarMsky(data, nCstar, badID):
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


def compMagAndFluxMean(data, nCstar, badID):
    """ Compute for each images the mean magnitude/flux sum off all the
        comparison stars, except for those whose ID is in the list badID
        (variable stars)
    """
    meanCmag = []
    sumCflux = []
    for line in data:
        magLine, fluxLine = [], []
        index = 1
        for i in xrange(5 + ap, 5 + ap + 3*nCstar, 3):
            if index not in badID:
                mag = line[i]
                magLine.append(mag)
                fluxLine.append(pow(10, 0.4*(25-mag)))
            index += 1
        meanCmag.append(np.mean(magLine))
        sumCflux.append(np.sum(fluxLine))
    return np.array(meanCmag), np.array(sumCflux)


def seeCstars(data, nCstar, badID):
    """ Plot raw and diff mag of each comp stars one by one"""
    for i in xrange(nCstar):
        cStarRawMag = data[:, 5 + ap + i*3]
        cStarDiffMag = cStarRawMag - meanCmag
        plt.plot(corrTime, cStarRawMag - np.mean(cStarRawMag), '.r', label='rawMag')
        plt.plot(corrTime, cStarDiffMag - np.mean(cStarDiffMag) + 1,
        '.k', label='diffMag')
        plt.title(str(i+1))
        plt.gca().invert_yaxis()
        plt.legend()
        plt.show()

def seeTenByTenCstars():
    """ Plot raw and diff mag of each comp stars 10 by 10"""
    for i in xrange(0, nCstar, 10):
        for j in xrange(10):
            cStarRawMag = data[:, 5 + ap + (i+j)*3]
            cStarDiffMag = cStarRawMag - meanCmag
            plt.plot(corrTime, cStarRawMag - np.mean(cStarRawMag) + j, '.r')
            plt.plot(corrTime, cStarDiffMag - np.mean(cStarDiffMag) + j + 0.5,
            '.k')
        plt.title(str(i))
        p*lt.legend()
        plt.show()


def meanDiffCstarAll():
    cStarDiffMagSumAll = np.zeros_like(data[:, 5 + ap])
    for i in xrange(nCstar):
        cStarDiffMagSumAll += data[:, 5 + ap + i*3] - meanCmag
    cStarDiffMagSumAll /= nCstar
    plt.plot(corrTime, cStarDiffMagSumAll - np.mean(cStarDiffMagSumAll), '.b')
    tStar = data[:, 5 + ap + 3] - meanCmag  #random star in comp star set
    plt.plot(corrTime, tStar - np.mean(tStar), '.g')
    plt.plot(corrTime, meanCmag - np.mean(meanCmag), '.r')

def meanDiffCstarGood():
    cStarDiffMagSum = np.zeros_like(data[:, 5 + ap])
    for i in xrange(nCstar):
        if i+1 not in badID:
            cStarDiffMagSum += data[:, 5 + ap + i*3] - meanCmag
    cStarDiffMagSum /= (nCstar-len(badID))
    plt.plot(corrTime, cStarDiffMagSum - np.mean(cStarDiffMagSum), '.k')
    plt.show()
