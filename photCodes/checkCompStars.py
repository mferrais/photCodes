import numpy as np, rtPhot as rtp, matplotlib.pyplot as plt
#---------------------------------------
field = 346
field = 345
field = 348

fieldExt = '_20171214'
fieldExt = '_20171215'
fieldExt = '_20171213'

night = '20171219'
night = '20171207'
night = '20171202'
night = '20171214'
night = '20171215'
night = '20171213'

telescope = 'TS'
telescope = 'TN'

ap = 2
ap = 3
ap = 18
ap = 13

filt = 13
filt = 15
filt = 12
filt = 14

fileSize = 'big'
fileSize = 'small'

m = 10000

fname = str(field)+'_'+telescope+'_'+night+'_'+str(ap)+'_'+str(filt)+'.hout1'
fname = str(field)+'_'+telescope+'_'+night+'_'+str(ap)+'_'+str(filt)+'.arc'

if telescope == 'TS':
    folder = "E:/trappist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/trappist/pho/"   # TS folder
if telescope == 'TN':
    folder = "E:/troppist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/troppist/pho/"   # TN folder
#------------------------------------------------------------------------------
folder += str(field)+fieldExt+'/'
print('folder : '+folder)
print('fname : '+fname)

arcData = np.loadtxt(folder + fname)

allID = arcData[:,1] - field*1e5
print(len(allID))
IDs = []
for ID in allID:
    if ID not in IDs:
        IDs.append(ID)
print(len(IDs))
print(IDs[0] == 1.0)
if IDs[0] == 1.0:
    IDs = IDs[1:]
print(len(IDs))

def seeTenByTenCstars(datas, m):
    ''' Plot raw and diff mag of each comp stars 10 by 10'''
    n = len(datas)
    for i in xrange(0, n, m):
        for j in xrange(m):
            if i+j < n:
                data = datas[i+j]
                mag = data[1]
                #mag -= (np.mean(mag)-(j+1))
                mag -= (np.mean(mag)-(j/20.))
                plt.plot(data[0], mag, '.-')
        plt.title(str(i))
        plt.show()

def methodSmallFile(m):
    print('small file method')
    datas = [rtp.getAstDataRacc(arcData, field, ID)[:2] for ID in IDs]
    print(len(datas))
    seeTenByTenCstars(datas, m)

def methodBigFile(m):
    print('big file method')
    n = len(IDs)
    for i in xrange(0, n, m):
        for j in xrange(m):
            if i+j >= n:
                break
            data = rtp.getAstDataRacc(arcData, field, IDs[i+j])[:2]
            mag = data[1] - (np.mean(data[1])-(j+1))
            plt.plot(data[0], mag, '.-')
        plt.title(str(i))
        plt.show()

if fileSize == 'small':
    methodSmallFile(m)
elif fileSize == 'big':
    methodBigFile(m)
