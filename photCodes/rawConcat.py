import numpy as np

#---------------------------------------
field = 341
field = 340
field = 342
field = 344
field = 347

telescope = 'TS'
telescope = 'TN'

obs = 'obs1'
obs = 'obs2'

if telescope == 'TS':
    folder = "E:/trappist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/trappist/pho/"   # TS folder
if telescope == 'TN':
    folder = "E:/troppist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/troppist/pho/"   # TN fo
if telescope == 'TNandTS':
    folder = "/media/marin/TRAPPIST3/northAndSouth/"
#---------------------------------------
folder += str(field)


apList = [1,2,3]
fnameL = ['/' + str(field) + '_raw_ap' + str(ap) + '_' + telescope+obs + '.txt'
          for ap in apList]
datas = [np.loadtxt(folder+fname) for fname in fnameL]
dataRaw = np.dstack((datas[0][:,0], datas[0][:,1],
                     datas[1][:,1], datas[2][:,1]))[0]

fnameRaw = '/' + str(field) + '_raw' + '_' + telescope+obs  + '.txt'

print('Raw file has been saved on folder : ' + folder)
print('With name: ' + fnameRaw[1:])
np.savetxt(folder + fnameRaw, dataRaw)
