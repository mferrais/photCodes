import numpy as np, matplotlib.pyplot as plt, rtPhot as rtp
from scipy.stats import sigmaclip
from scipy.interpolate import UnivariateSpline
from astroML.filters import savitzky_golay
from astroML.time_series import MultiTermFit

#---------------------------------------
ap = 3
ap = 18
ap = 1
ap = 13
ap = 2

field = 341
field = 342
field = 344
field = 347
field = 340
field = 343
field = 346
field = 348
field = 345

fieldExt = 'n2'
fieldExt = 'n1'
fieldExt = '_20171214'
fieldExt = '_20171215'
fieldExt = '_20171219'
fieldExt = '_20171213'
fieldExt = ''

telescope = 'TS'
telescope = 'TN'
telescope = 'TNandTS'

obs = 'jm'
obs = 'obs2'
obs = 'aij'
obs = 'obs1'
obs = 'obs1B'
obs = 'obs1V'
obs = 'obs1R'
obs = 'obs1I'
obs = 'obs2BVRI'
obs = 'obs1BVRI'
obs = 'obs1Bracc'
obs = 'obs2Bracc'
obs = 'obs2Iracc'
obs = '18'
obs = 'BVRIracc'
obs = 'obs1And2'
obs = '20171214VR'
obs = '20171215IB'
obs = 'obs123'
obs = 'racc'
obs = 'all'

filt = 12
filt = 15
filt = 13
filt = 'BVRI'
filt = 'IB'
filt = 14

night = '20171218'
night = '20171207'
night = '20171202'
night = '20171214'
night = '20171215'
night = '20171219'
night = '20171213'

clean = True
clean = False

fileOut = True
fileOut = False


if field == 340:
    Filter = 'R'
    title = '596 Scheila'
    subtitle = '3 - 4 Jun 2017, TRAPPIST South'
    ProtPos = 0.61, 1.035
    tShift = 0.005
    bShift = 2*tShift
    nights = ['20170602', '20170603', 'old']
    nights = ['20170602TS', '20170603TS']
    lList = [690,870] # uncleaned
    lList = [690,753] # last cleaned
    lList = [673,1445-673] # jm cleaned
    lList = [690,1426-690] # last cleaned ap3
    lList = [1344,3072-1344] # w596_weighted
    lList = [673,1534-673] # last cleaned ap2
    lList = [690,1443-690] # last cleaned ap2
    Prot = 16.0
    Prot = 15.84400,730,783
    Prot = 15.848
    Prot = 15.793
elif field == 342:
    Filter = 'R'
    title = '476 Hedwig'
    subtitle = '19 - 12 Jul 2017, TRAPPIST South'
    ProtPos = 0.25, 1.065
    tShift = 0.005
    bShift = 0.01
    nights = ['20170618', '20170628', '20170709', '20170711']
    Prot = 27.33
    Prot = 27.336
    lList = [720,46,522,161]
elif field == 343:
    Filter = 'V'
    title = '3122 Florence'
    ProtPos = 0.6, 1.09
    tShift = 0.01
    bShift = 2*tShift
    Prot = 2.3685
    Prot = 2.3581
    if telescope == 'TS':
        subtitle = '4 Sep 2017, TRAPPIST South'
        nights = ['20170618','20170618','20170618','20170618','20170618','20170618']
        lList = [69, 74, 67, 43, 50, 60]
    elif telescope == 'TN':
        subtitle = '4 Sep 2017, TRAPPIST North'
        nights = ['20170603TN','20170603TN','20170603TN','20170603TN','20170603TN','20170603TN','20170603TN','20170603TN']
        nights = ['20170903TN']
        lList = [78, 55, 50, 67, 100, 100, 73, 60]
        lList = [580]
    elif telescope == 'TNandTS':
        subtitle = '4 Sep 2017, TRAPPIST North and TRAPPIST South'
        nights = ['20170903TN','20170903TS','20170903TS','20170903TS','20170903TS','20170903TS','20170903TS']
        nights = ['20170903TN','20170903TS','20170903TS','20170903TS','20170903TS']
        nights = ['20170903TN','20170903TS']
        nights = ['V20170903TN','V20170903TS','I20170903TN']
        lList = [583, 56, 57, 55, 100]
        lList = [583, 268]
        lList = [583, 268, 18]
elif field == 344:
    Filter = 'R'
    title = '24 Themis'
    subtitle = '8 - 9 Nov 2017, TRAPPIST South and 9 Nov 2017, TRAPPIST North'
    subtitle = '21-24 Sep + 8-9 Nov 2017, TRAPPIST South and 9 Nov 2017, TRAPPIST North'
    ProtPos = 0.2, 1.065
    tShift = 0.005
    bShift = 2*tShift
    nights = ['20171107', '20171108']
    nights = ['20171108']
    nights = ['20171107TS', '20171108TS', '20171108TN']
    nights = ['20170920TS','20170921TS', '20170922TS', '20170923TS']
    nights = ['20170920TS', '20170921TS', '20170922TS', '20170923TS',
              '20171107TS', '20171108TS', '20171108TN']
    MFtimes = [2458066.4485532409] # 344 20171109 TN
    Prot = 8.374
    lList = [707,699]
    lList = [76]
    lList = [707,699,76]
    lList = [700,698,62]
    lList = [255,315,308,362]
    lList = [255,315,308,362,707,699,664]
    lList = [255,315,308,362,701,699,651] # cleaned 1
    lList = [255,315,308,362,698,698,628] # cleaned 2
elif field == 345:
    Filter = ' '
    title = '89 Julia'
    ProtPos = 0.2, 0.94
    tShift = 0.01
    bShift = 2*tShift
    Prot = 11.38834
    Prot = 11.387
    Prot = 11.3844
    if obs == 'obs1':
        MFTimes = [0.6]
        subtitle = '20 - 30 Sep 2017, TRAPPIST South'
        nights = ['20170920','20170921','20170922','20170930']
        lList = [553,537,619,726]
    elif obs == 'all':
        subtitle = '20 - 30 Sep 2017, TRAPPIST South ; 27 Nov - 24 Dec 2017 TRAPPIST North'
        nights = ['V20170920TS','V20170921TS','V20170922TS','V20170930TS','R20171127TN','R20171202TN','R20171207TN','R20171223TN']
        lList = [553,537,619,726,522,478,604,418]
    if obs == 'racc' and night == '20171202':
        subtitle = '2 Dec TRAPPIST North'
        nights = ['20171202TN']
        lList = [478]
elif field == 346:
    Filter = 'R'
    title = '31 Euphrosyne'
    ProtPos = 0.6, 1.04
    tShift = 0.005
    bShift = tShift*2
    if obs == 'obs1And2':
        #MFTimes = [0.6]
        subtitle = '07 Nov - 02 Dec 2017, TRAPPIST North'
        nights = ['20171107TN','20171127TN','20171202TN']
        lList = [374-15,261,365]
    if obs == 'obs123':
        subtitle = '07 Nov - 08 Dec 2017, TRAPPIST North'
        nights = ['20171107TN','20171127TN','20171202TN','20171207TN']
        lList = [359,261,365,864]
    if obs == 'racc' and night == '20171207':
        subtitle = '7 Dec TRAPPIST North'
        nights = ['20171207TN']
        lList = [1000]
        lList = [864]
    Prot = 5.53
    Prot = 5.57991
    Prot = 5.529597
    Prot = 5.52874
    Prot = 5.5312
elif field == 347:
    Filter = 'R'
    title = '20 Massalia'
    ProtPos = 0.4, 1.3
    tShift = 0.01
    tShift = 0.05
    bShift = 2*tShift
    if obs == 'all':
        MFTimes = [2458068.636574,2458068.810945,2458076.787568]
        subtitle = '11 + 18 Nov 2017, TRAPPIST South and TRAPPIST North'
        nights = ['20171110TN','20171110TS','20171118TN','20171118TS']
        lList = [700,729,700,783]
        lList = [600,729,700,783]
    Prot = 8.098
elif field == 348:
    Filter = ''
    title = '3200 Phaethon'
    ProtPos = 0.4, 1.075
    ProtPos = 0.4, 1.2
    tShift = 0.01
    bShift = 2*tShift
    if obs == 'obs1R':
        nights = ['R1','R2','R3','R4','R5']
        lList = [46, 50, 48, 50, 48]
    elif obs == 'obs1V':
        nights = ['V1','V2','V3','V4','V5','V6']
        lList = [41, 40, 40, 45, 35, 35]
    elif obs == 'obs1B':
        nights = ['B1','B2','B3','B4','B5','B6']
        lList = [22, 24, 27, 24, 24, 21]
    elif obs == 'obs1I':
        nights = ['I1','I2','I3','I4','I5','I6']
        lList = [41, 40, 40, 40, 40, 35]
    elif obs == 'obs1BVRI' or obs == 'BVRIracc':
        subtitle = '9 Dec 2017, TRAPPIST North'
        nights = ['B','V','R','I']
        lList = [142,236,242,236]
    elif obs == 'obs2BVRI':
        subtitle = '10 Dec 2017, TRAPPIST North'
        nights = ['B','V','R','I']
        lList = [122,205,225,210]
    elif obs == 'obs1Bracc':
        subtitle = '9 Dec 2017, TRAPPIST North'
        nights = ['B']
        lList = [142]
    elif obs == 'obs2Iracc':
        subtitle = '10 Dec 2017, TRAPPIST North'
        nights = ['B']
        nights = ['V']
        nights = ['R']
        nights = ['I']
        lList = [122]  # B
        lList = [205]  # V
        lList = [225]  # R
        lList = [210]  # I
    elif obs == 'racc' and night == '20171218Not':
        subtitle = '18 Dec 2017, TRAPPIST North'
        nights = ['R']
        lList = [294]  # R
    elif obs == 'racc' and night == '20171218':
        Filter  = ' R '
        subtitle = '18 and 19 Dec 2017, TRAPPIST North'
        nights = ['20171218TN','20171219TN']
        lList = [368,299]  # R
    elif obs == 'racc' and night == '20171219':
        Filter  = ' R '
        subtitle = '19 Dec 2017, TRAPPIST North'
        nights = ['20171219TN']
        lList = [299]  # R
    elif obs == 'racc' and night == '20171213':
        Filter = 'R'
        subtitle = '13 Dec 2017, TRAPPIST North'
        nights = ['20171213TN',]
        lList = [1353]  # R
        lList = [1346]  # R
    elif obs == 'racc' and night == '20171214':
        subtitle = '14 Dec 2017, TRAPPIST North'
        nights = ['B20171214TN']
        nights = ['I20171214TN']
        nights = ['V20171214TN']
        nights = ['B','V','R','I']
        lList = [24]  # B - I
        lList = [460]  # V
        lList = [360]  # R
        lList = [24,460,360,24]  # BVRI
    elif obs == '20171214VR':
        Filter = ''
        subtitle = '14 Dec 2017, TRAPPIST North'
        nights = ['R20171214TN','V20171214TN']
        lList = [369,496]
    elif obs == '20171213R':
        Filter = 'R'
        subtitle = '13 Dec 2017, TRAPPIST North'
        nights = ['20171213TN']
        lList = [1350]
    elif obs == 'racc' and night == '20171215':
        subtitle = '15 Dec 2017, TRAPPIST North'
        nights = ['I20171215TN','B20171215TN']
        lList = [444,368]  # IB
    elif obs == '20171215IB' and night == '20171215':
        subtitle = '15 Dec 2017, TRAPPIST North'
        nights = ['I20171215TN','B20171215TN']
        lList = [447,405]  # IB
    Prot = 3.603957

wsSG = 601

ffOrder = 7

#binning
binData = True
binData = False
binW = 2 # bin width in minutes

nStart = 1
gLine = True
gLine = False

correction = True
correction = False

yAxisLim = (0.5, 1.5)
yAxisLim = (0.85, 1.15)
yAxisLim = (0.9, 1.10)
yAxisLim = (12.2, 12.7)
setylim = True
setylim = False

setToMag = True
setToMag = False

showMF = True
showMF = False

if telescope == 'TS':
    folder = "E:/trappist/pho/"                    # win TS folder
    folder = "/media/marin/TRAPPIST3/trappist/pho/"   # TS folder
if telescope == 'TN':
    folder = "E:/troppist/pho/"                       #clean win TS folder
    folder = "/media/marin/TRAPPIST3/troppist/pho/"   # TN fo
if telescope == 'TNandTS':
    folder = "/media/marin/TRAPPIST3/northAndSouth/"
#----------------
folder += str(field) + fieldExt +'/'
plt.style.use('bmh')

#fname = '/' + str(field) + '_diff_ap' + str(ap) + 'Damit.txt'

if clean:   #cleaned data
    fname = str(field) + '_diff_ap' + str(ap) + '_clean.txt'
#elif obs == 'obs1':
#    fname = '/'+str(field)+'_diff_ap'+str(ap)+'_'+telescope+obs+'.txt'
elif obs == 'all':  # uncleaned data
    fname = str(field) + '_diff_ap' + str(ap) + '.dat'
elif obs == 'aij':
    fname = 'AIJ/03122AIJ.txt'
elif obs == 'BVRIracc':
    fname = str(field)+'_'+telescope+'_'+night+'_'+str(ap)+'_BVRI'+'.mag'
elif obs == 'racc':
    fname = str(field)+'_'+telescope+'_'+night+'_'+str(ap)+'_'+str(filt)
    if setToMag:
        fname += '.mag'
    else:
        fname += '.flux'
else:
    #fname = '/'+str(field)+'_diff_ap'+str(ap)+'_'+telescope+obs+'.txt'
    fname = str(field)+'_diff_ap'+str(ap)+'_'+telescope+obs+'.dat'

data = np.loadtxt(folder + fname)
time = data[:,0]
astDiffFlux = data[:,1]
astFluxErr = data[:,2]

if obs == 'jm':
    setToMag = True
    folder = '/home/marin/Downloads/Fwd_Scheila_LC/'
    fname = 'w596_weighted.txt'
    data = np.loadtxt(folder + fname)
    time = data[:,0]
    astDiffFlux = data[:,1]

nPts = len(astDiffFlux)
print(fname)
print(str(nPts) + ' nPts')

#----- detrend ----------------------------------------------------------------
"""
# fit a straight line to correct trend due to phase angle changes
l0 = 0
for i in range(len(nights)):
    li = lList[i]
    flux, t= astDiffFlux[l0:l0+li], time[l0:l0+li]

    spl = UnivariateSpline(t, flux,k=3)
    detrendFlux = flux / spl(t)
    plt.plot(t,flux,'C0.', t,spl(t),'k--', t,detrendFlux,'C1.')
    plt.show()
    '''
    pfit = np.poly1d(np.polyfit(t, flux, 2))
    detrendFlux = flux / pfit(t)
    plt.plot(t,flux,'.', t,pfit(t),'k--', t,detrendFlux,'.')
    plt.show()
    '''
    astDiffFlux[l0:l0+li] = detrendFlux
    l0 +=li
    #astDiffFlux[l0:l0+li] -= np.mean(astDiffFlux[l0:l0+li])
astDiffFlux = data[:,1]
"""
"""
# use the overlaping part in JD of two LC to match their fluxes
l0 = lList[0]
f1, time1 = astDiffFlux[0:l0], time[0:l0]
for i in range(1, len(nights)):
    li = lList[i]
    f2, time2 = astDiffFlux[l0:l0+li], time[l0:l0+li]
    shift = rtp.meanShift(f1, f2, time1, time2)
    astDiffFlux[l0:l0+li] += shift
    l0 += li
    print(str(shift) + ' : shift ' + nights[i])

astDiffFlux /= np.mean(astDiffFlux)
"""


#----- Correction -------------------------------------------------------------
i = 0
ind = []
for l in lList:
    i += l
    ind.append(i)

print(ind, 'ind')
if correction:
    if field == 343:
        #astDiffFlux = -2.5 * np.log10(astDiffFlux)
        print(field,obs)
        astDiffFlux[ind[1]:ind[2]] -= 0.065
    elif field == 344 and obs == 'all':
        print(field,obs)
        astDiffFlux[:ind[0]] -= 0.0035
        astDiffFlux[ind[3]:ind[4]] += 0.015
        astDiffFlux[ind[4]:ind[5]] += 0.010
        astDiffFlux[2685+250:] += 0.01
    elif field == 345:
        print(field,obs)
        astDiffFlux[:ind[0]] -= 0.0
        astDiffFlux[ind[0]:ind[1]] -= 0.0
        astDiffFlux[ind[1]:ind[2]] -= 0.0
        astDiffFlux[ind[2]:ind[3]] -= 0.0
        astDiffFlux[ind[3]:ind[4]] -= 0.0
        astDiffFlux[ind[4]:ind[5]] -= 0.0
        astDiffFlux[ind[5]:ind[6]] -= 0.0
        astDiffFlux[ind[6]:ind[7]] -= 0.0
    elif field == 346:
        print(field,obs)
    elif field == 347:
        print(field,obs)
        astDiffFlux[ind[0]:ind[1]] += 0.016
        astDiffFlux[ind[1]:ind[2]] -= 0.019
        astDiffFlux[ind[2]:ind[3]] -= 0.010
    elif field == 348 and obs == 'obs1B':
        print(field,obs)

    elif field == 348 and obs == 'racc' and not setToMag and night == '20171214':
        print(field,obs)

    elif field == 348 and obs == 'racc' and setToMag and night == '20171214':
        print(field,obs)
        astDiffFlux[:ind[0]] += 13 - np.mean(astDiffFlux[:ind[0]])
        astDiffFlux[ind[0]:ind[1]] += 13 - np.mean(astDiffFlux[ind[0]:ind[1]])
        astDiffFlux[ind[1]:ind[2]] += 13 - np.mean(astDiffFlux[ind[1]:ind[2]])
        astDiffFlux[ind[2]:ind[3]] += 13 - np.mean(astDiffFlux[ind[2]:ind[3]])
    elif field == 348 and obs == '20171214VR':
        print(field,obs)
        astDiffFlux[ind[0]:ind[1]] += 0.04
    elif field == 348 and obs == 'racc' and setToMag and night == '20171215':
        print(field,obs)
        astDiffFlux[ind[0]:ind[1]] += 0.0
    elif field == 348 and obs == '20171215IB' and night == '20171215':
        print(field,obs)
        astDiffFlux[ind[0]:ind[1]] -= 0.0

    if setToMag:
        #astDiffFlux -= np.mean(astDiffFlux)
        bob = 'bob'
    else:
        astDiffFlux /= np.mean(astDiffFlux)

#----- data smoothing ---------------------------------------------------------
"""
span = 10
for i in xrange(nPts):
    pts = astDiffFlux[i]
    mean = np.mean(astDiffFlux[i-span:i+span])
fluxSmooth = savitzky_golay(astDiffFlux, window_size=201, order=4)
plt.plot(time,astDiffFlux,'C1.', time,fluxSmooth,'k.')
plt.show()
"""

#------ Data Binning ----------------------------------------------------------
# index list
i = 0
ind = []
for l in lList:
    i += l
    ind.append(i)

print(lList, 'lList')
if binData:
    oldTime, oldAstDiffFlux,  = time, astDiffFlux
    oldAstFluxErr = astFluxErr
    astDiffFlux = rtp.binning(astDiffFlux, time, binW, ind)[0]
    astFluxErr = rtp.binning(astFluxErr, time, binW, ind)[0]
    resBin = rtp.binWeighted(time, time, oldAstFluxErr, binW, ind)
    time, lList = resBin[0], resBin[1]
print(len(time), ' new nPts')
print(lList, 'new lList')

# phase
nind = [0] + ind
t0 = time[nind[nStart-1]]
if field == 347:
    t0 = 2458068.423936
if field == 348:
    t0 = 2458097.280053
ProtJD = Prot/24.
phase = [((t - t0)/ProtJD) % 1. for t in time]

d = 60
if binData:
    nind = [0] + ind
    t0 = oldTime[nind[nStart-1]]
    if field == 347:
        t0 = 2458068.423936
    ProtJD = Prot/24.
    oldPhase = [((t - t0)/ProtJD) % 1. for t in oldTime]
    fluxByPhase = [f for _,f in sorted(zip(oldPhase,oldAstDiffFlux))]
    errByPhase = [f for _,f in sorted(zip(oldPhase,oldAstFluxErr))]
    timeByPhase = [f for _,f in sorted(zip(oldPhase,oldTime))]
    sortedPhase = sorted(oldPhase)
    sortedPhaseBin = sorted(rtp.binning(oldPhase, oldTime, binW, ind)[0])
    fluxByPhaseBin = rtp.binWeighted(fluxByPhase,oldTime,oldAstFluxErr,binW,ind)[0]
    weights = [ 1./(err**2) for err in errByPhase]
    pfit = np.poly1d(np.polyfit(sortedPhase, fluxByPhase, deg=d, w=weights))
    fitCurve = pfit(sortedPhase)
    strCurve = fluxByPhase / fitCurve
else:
    fluxByPhase = [f for _,f in sorted(zip(phase,astDiffFlux))]
    errByPhase = [f for _,f in sorted(zip(phase,astFluxErr))]
    sortedPhase = sorted(phase)
    sortedPhaseBin = sorted(rtp.binning(phase, time, binW, ind)[0])
    fluxByPhaseBin = rtp.binWeighted(fluxByPhase, time, astFluxErr,binW, ind)[0]
    weights = [ 1./(err**2) for err in errByPhase]
    pfit = np.poly1d(np.polyfit(sortedPhase, fluxByPhase, deg=d, w=weights))
    fitCurve = pfit(sortedPhase)
    strCurve = fluxByPhase / fitCurve

fluxByPhaseSG = savitzky_golay(fluxByPhase, window_size=wsSG, order=4)

#----- Fourier Fit ------------------------------------------------------------
omega = 2 * np.pi / Prot * 24.
mtf = MultiTermFit(omega, ffOrder)
mtf.fit(time, astDiffFlux, astFluxErr)
phaseFit, fluxFit, phasedTime = mtf.predict(1000, return_phased_times=True)

phaseFit -= phasedTime[0]
phasedTime -= phasedTime[0]
phasedTime = [1+t if t < 0 else t for t in phasedTime]
phaseFit = [1+t if t < 0 else t for t in phaseFit]

plt.errorbar(phasedTime, astDiffFlux, astFluxErr, fmt='.k', ecolor='gray',
             lw=1, ms=4, capsize=1.5)
plt.plot(phaseFit, fluxFit, 'b.', lw=10, label='order = '+str(ffOrder))

plt.xlabel('Rotational Phase')
plt.ylabel('Relative ' + Filter + ' Flux')
plt.suptitle(title, y=0.945, fontsize=18)
plt.title(subtitle, fontsize=10)
if setToMag:
    plt.gca().invert_yaxis()
    plt.ylabel('Differential ' + Filter + ' Mag ')
if setylim:
    plt.ylim(yAxisLim)
plt.xlim(-0.01, 1.01)
plt.tick_params(direction='in', right='on',top='on', bottom='on')
plt.text(ProtPos[0], ProtPos[1], 'Prot = '+str(Prot)+' h', fontsize=12)
plt.text(ProtPos[0], ProtPos[1]-tShift, 'Zero time = '+str('%.6f' % t0)+' JD',
         fontsize=12)
plt.minorticks_on()
plt.legend(loc='best')
plt.show()

#---- saving on file ----------------------------------------------------------
if fileOut and binData:
    print('Warning : data are binned !')
if fileOut and not binData:
    print('* File has been saved on folder : ' + folder)
    print('With name: ' + fname)
    data = np.dstack((time, astDiffFlux, astFluxErr))[0]
    np.savetxt(folder + fname, data, fmt=['%10.6f','%.6e','%.2e'])
    """
    ind = [0] + ind
    TScount = 0
    TNcount = 0
    for night in nights:
        count = TScount + TNcount
        if night[-1] == 'S':
            TScount += 1
            fnameN = fname[:-4] + '_' + night[-02:]+'obs'+str(TScount)+'_FM.dat'
        else:
            TNcount += 1
            fnameN = fname[:-4] + '_' + night[-2:]+'obs'+str(TNcount)+'_FM.dat'
        print('* Night file has been saved on folder : ' + folder)
        print('With name: ' + fnameN)
        dataN = data[ind[count]:ind[count+1]]
        magN = -2.5 * np.log10(dataN[:,1])
        magErrN = -2.5 * np.log10(dataN[:,2])
        dataN = np.dstack((dataN[:,0], dataN[:,1], dataN[:,2], magN, magErrN))[0]
        np.savetxt(folder + fnameN, dataN, fmt=['%10.6f','%.6e','%.2e','%.6e','%.2e'])
    """
#------- Graphs ---------------------------------------------------------------
#astDiffFlux = -2.5 * np.log10(astDiffFlux)

"""
u = np.arange(len(sortedPhaseBin))
plt.plot(u, sortedPhaseBin, '.')
plt.show()
"""
months = {'01':'Jan','02':'Feb','03':'Mar','04':'Apr','05':'May','06':'Jun',
          '07':'Jul','08':'Aug','09':'Sep','10':'Oct','11':'Nov','12':'Dec'}
colors = ['C0.','C1.','C2v','C3*','C4.','C5.','C6.','C7.']
colors = ['C0p','C1.','C2^','C3*','C4h','C5v','C6s','C7x']
colors = ['C0p','C1,','C2,','C3*','C4h','C5v','C6s','C7x']
colors = ['C0.','C1.','C2.','C3.','C4.','C5.','C6.','C7.','C8.','C9.']

colors = ['#0082c8','#3cb44b','#f58231','#800000'] # BVRI colors
mark = ['p','.','^','*','h','v','s','x']

colors = ['b','g','r','#800000'] # BVRI colors
colors = ['r','g'] # RV colors
colors = ['#800000','b'] # IB colors
colors = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9']
markDot = ['.' for c in colors]
mark = ['h','s','v','p']
mark = ['h','o','v','p','s']
mark = markDot

colorsLine = ['C0-','C1-','C2-','C3-','C4-','C5-','C6-','C7-']
labels = []
for i in range(len(nights)):
        if len(nights[0]) == 10:
            label = nights[i][-2:] + ' ' + nights[i][0:4] + ' '
            label += months[str(nights[i][4:6])]+' '+str(int(nights[i][-4:-2]))
            #label += months[str(nights[i][4:6])]+' '+str(int(nights[i][-4:-2]))
            labels.append(label)
        elif len(nights[0]) == 11:
            label = nights[i][-2:] + ' ' + nights[i][1:5] + ' '
            label += months[str(nights[i][5:7])]+' '
            label += str(int(nights[i][-4:-2]))+'  ('+nights[i][0] + ')'
            labels.append(label)
        else:
            labels.append(nights[i])

# Plot the whole phased curve
lab = ''.join([n + ', ' for n in labels])
#plt.plot(sortedPhase, fluxByPhase, colors[1], label=lab)
plt.plot(sortedPhaseBin, fluxByPhaseBin, 'C1.', label=lab)
#plt.plot(sortedPhaseBin, fluxByPhaseBin, colors[1])
plt.xlabel('Rotational Phase')
plt.ylabel('Relative ' + Filter + ' Flux')
plt.suptitle(title, y=0.945, fontsize=18)
plt.title(subtitle, fontsize=10)
if setToMag:
    plt.gca().invert_yaxis()
    plt.ylabel('Differential ' + Filter + ' Mag ')
if setylim:
    plt.ylim(yAxisLim)
plt.xlim(-0.01, 1.01)
plt.tick_params(direction='in', right='on',top='on', bottom='on')
plt.text(ProtPos[0], ProtPos[1], 'Prot = '+str(Prot)+' h', fontsize=12)
plt.text(ProtPos[0], ProtPos[1]-tShift, 'Zero time = '+str('%.6f' % t0)+' JD',
         fontsize=12)
plt.text(ProtPos[0],ProtPos[1]-bShift,'Bin = '+str(binW)+' min', fontsize=12)
plt.plot(sortedPhase, fluxByPhaseSG, 'k.', label='Savitzky_Golay')
plt.legend(loc='best')
plt.show()

"""
# test fit
plt.plot(sortedPhase, fluxByPhase, 'r.', label='flux')
plt.plot(sortedPhase, fitCurve, 'k-', label='fitted curve')
plt.plot(sortedPhase, strCurve, 'b.', label='redressed curve')
plt.xlabel('Rotational Phase')
plt.ylabel('Relative ' + Filter + ' Flux')
plt.suptitle(title, y=0.945, fontsize=18)
plt.title(subtitle, fontsize=10)
plt.legend(loc='best')
if setToMag:
    plt.gca().invert_yaxis()
    plt.ylabel('Differential ' + Filter + ' Mag ')
if setylim:
    plt.ylim(yAxisLim)
plt.xlim(-0.01, 1.01)
plt.tick_params(direction='in', right='on',top='on', bottom='on')
plt.text(ProtPos[0], ProtPos[1], 'Prot = '+str(Prot)+' h', fontsize=12)
plt.text(ProtPos[0], ProtPos[1]-tShift, 'Zero time = '+str('%.6f' % t0)+' JD',
         fontsize=12)
plt.text(ProtPos[0],ProtPos[1]-bShift,'Bin = '+str(binW)+' min', fontsize=12)
plt.show()
"""

# Plot flux or mag in a sequential way
l0 = 0
for i in range(len(nights)):
    li = lList[i]
    plt.plot(np.arange(li)+l0+1, astDiffFlux[l0:li+l0], '.', label=labels[i])
    l0 += li
plt.xlabel('n')
plt.ylabel('Relative ' + Filter + ' Flux')
if setToMag:
    plt.gca().invert_yaxis()
    plt.ylabel('Differential ' + Filter + ' Mag ')
plt.title('Sequential plot')
plt.legend()
plt.show()

"""
# Plot error in a sequential way
l0 = 0
for i in range(len(nights)):
    li = lList[i]
    plt.plot(np.arange(li)+l0+1, astFluxErr[l0:li+l0], '.', label=labels[i])
    l0 += li
plt.xlabel('n')
plt.ylabel('error')
plt.title('Error sequential plot')
plt.legend()
plt.show()
"""

# Plot flux or mag vs rotational phase
l0 = 0
if gLine:
    colors = colorsLine
for i in range(len(nights)):
    li = lList[i]
    plt.plot(phase[l0:li+l0], astDiffFlux[l0:li+l0], color=colors[i],
             marker=mark[i], linestyle='None' ,label=labels[i])
    #plt.plot(phase[l0:li+l0], astDiffFlux[l0:li+l0], colors[i])
    #plt.errorbar(phase[l0:li+l0], astDiffFlux[l0:li+l0]
    #              , yerr=astFluxErr[l0:li+l0], fmt=colors[i])
    l0 += li
plt.xlabel('Rotational Phase')
plt.ylabel('Relative ' + Filter + ' Flux')
plt.suptitle(title, y=0.945, fontsize=18)
plt.title(subtitle, fontsize=10)
if setylim:
    plt.ylim(yAxisLim)
if setToMag:
    plt.gca().invert_yaxis()
    plt.ylabel('Differential ' + Filter + ' Mag ')
plt.xlim(-0.01, 1.01)
plt.tick_params(direction='in', right='on',top='on', bottom='on')
plt.text(ProtPos[0], ProtPos[1], 'Prot = '+str(Prot)+' h', fontsize=12)
plt.text(ProtPos[0], ProtPos[1]-tShift, 'Zero time = '+str('%.6f' % t0)+' JD',
         fontsize=12)
if binData:
    plt.text(ProtPos[0],ProtPos[1]-bShift,'Bin = '+str(binW)+' min',
             fontsize=12)
if showMF:
    MFphases = [(t - t0)/ProtJD % 1. for t in MFTimes]
    for i in range(len(MFphases)):
        x1, x2 = MFphases[i], MFphases[i]
        plt.plot((x1,x2),(min(astDiffFlux),max(astDiffFlux)),colorsLine[i])
'''
d = 2458132
u = np.arange(d+0.3,d+0.4,0.001)
phaseT = [((t - t0)/ProtJD) % 1. for t in u]
plt.plot(phaseT, [1 for x in phaseT], 'k.')
'''
plt.plot(sortedPhase, fluxByPhaseSG, 'k.', label='Savitzky_Golay')
plt.legend(loc='best')
plt.show()

### Plot flux or mag vs JD + Error ############################################
"""
plt.figure(1)
plt.subplot(211)
#ax1 = plt.subplot(211)
l0 = 0
for i in range(len(nights)):
    li = lList[i]
    plt.plot(time[l0:li+l0]-2450000, astDiffFlux[l0:li+l0], color=colors[i],
             marker=markDot[i], linestyle='None' ,label=labels[i])
    #plt.plot(time[l0:li+l0]-2450000, astDiffFlux[l0:li+l0],'.')
    l0 += li
#plt.xlabel('Julian Day  - 2450000')
plt.ylabel('Relative ' + Filter + ' Flux')
plt.suptitle(title, y=0.945, fontsize=18)
plt.title(subtitle, fontsize=10)
plt.legend(loc='best')
if setylim:
    plt.ylim(yAxisLim)
if setToMag:
    plt.gca().invert_yaxis()
    plt.ylabel('Differential ' + Filter + ' Mag ')
#plt.setp(ax1.get_xticklabels(), visible=False)
plt.subplot(212)
#ax2 = plt.subplot(212, sharex=ax1)
l0 = 0
for i in range(len(nights)):
    li = lList[i]
    plt.plot(time[l0:li+l0]-2450000, astFluxErr[l0:li+l0], color=colors[i],
             marker=markDot[i], linestyle='None' ,label=labels[i])
    #plt.plot(time[l0:li+l0]-2450000, astFluxErr[l0:li+l0], '.')
    l0 += li
plt.xlabel('Julian Day - 2450000')
plt.ylabel('Error')
#plt.setp(ax2.get_xticklabels())
#plt.subplots_adjust(left=0.2, wspace=0.8, top=0.8)
plt.legend(loc='best')
"""
###############
f, (ax1, ax2) = plt.subplots(2, 1, sharex='all', gridspec_kw = {'height_ratios':[3, 1]})
l0 = 0
for i in range(len(nights)):
    li = lList[i]
    ax1.scatter(time[l0:li+l0]-2450000, astDiffFlux[l0:li+l0], color=colors[i],
             marker=mark[i], linestyle='None' ,label=labels[i])
    #ax1.plot(time[l0:li+l0]-2450000, astDiffFlux[l0:li+l0],'.')
    l0 += li
ax1.set_ylabel('Relative ' + Filter + ' Flux')
plt.suptitle(title, y=0.935, fontsize=18)
ax1.set_title(subtitle, fontsize=10)
ax1.legend(loc='best')
ax2.set_ylabel('Flux Error')
if setylim:
    ax1.set_ylim(yAxisLim)
if setToMag:
    ax1.invert_yaxis()
    ax1.set_ylabel('Differential ' + Filter + ' Mag ')
    ax2.set_ylabel('Mag Error')
# error
l0 = 0
for i in range(len(nights)):
    li = lList[i]
    ax2.plot(time[l0:li+l0]-2450000, astFluxErr[l0:li+l0], color=colors[i],
             marker=mark[i], linestyle='None' ,label=labels[i])
    #ax2.plot(time[l0:li+l0]-2450000, astFluxErr[l0:li+l0], '.')
    l0 += li
ax2.set_xlabel('Julian Day - 2450000')
ax1.minorticks_on()
#ax2.legend(loc='best')
plt.show()


"""
s = 2

plt.figure(1)
plt.subplot(221)
l0 = 0
for i in range(0,len(nights)-s):
    li = lList[i]
    plt.plot(time[l0:li+l0]-2450000, astDiffFlux[l0:li+l0],'.',label=labels[i])
    l0 += li
plt.xlabel('Julian Day  - 2450000')
plt.ylabel('Relative ' + Filter + ' Flux')
plt.title(title+' - 11 Novembre')
plt.legend(loc='best')
if setylim:
    plt.ylim(yAxisLim)
if setToMag:
    plt.gca().invert_yaxis()
    plt.ylabel('Differential ' + Filter + ' Mag ')
plt.subplot(222)
l0 = sum(lList[:s])
for i in range(s,len(nights)):
    li = lList[i]
    plt.plot(time[l0:li+l0]-2450000, astDiffFlux[l0:li+l0],'.',label=labels[i])
    l0 += li
plt.xlabel('Julian Day  - 2450000')
plt.ylabel('Relative ' + Filter + ' Flux')
plt.title(title+' - 19 Novembre')
plt.legend(loc='best')
if setylim:
    plt.ylim(yAxisLim)
if mag:
    plt.gca().invert_yaxis()
    plt.ylabel('Differential ' + Filter + ' Mag ')
plt.subplot(223)
l0 = 0
for i in range(0,len(nights)-s):
    li = lList[i]
    plt.plot(time[l0:li+l0]-2450000, astFluxErr[l0:li+l0], '.',label=labels[i])
    l0 += li
plt.xlabel('Julian Day - 2450000')
plt.ylabel('Error')
plt.legend(loc='best')

plt.subplot(224)
l0 = sum(lList[:s])
for i in range(s,len(nights)):
    li = lList[i]
    plt.plot(time[l0:li+l0]-2450000, astFluxErr[l0:li+l0], '.',label=labels[i])
    l0 += li
plt.xlabel('Julian Day - 2450000')
plt.ylabel('Error')
plt.legend(loc='best')
plt.show()
"""
