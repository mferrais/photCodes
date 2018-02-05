import numpy as np

field = 341
field = 340
field = 342
field = 344
field = 345
field = 347
field = 346

telescope = 'TS'
telescope = 'TN'

night = '20171127'
night = '20171202'

cut1 = 0
cut2 = 47145
cut1 = 47144
cut2 = 119110

#------------------------------------------------------------------------------
if telescope == 'TS':
    folder = "E:/trappist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/trappist/pho/"   # TS folder
if telescope == 'TN':
    folder = "E:/troppist/pho/"                       # win TS folder
    folder = "/media/marin/TRAPPIST3/troppist/pho/"   # TN folder
folder += str(field)
fname = '/h.runout1Crop'

data = np.loadtxt(folder+fname)[cut1:cut2-1]
print(data.shape)

fname += '_' + str(field) + '_' + str(night)
np.savetxt(folder + fname, data, fmt=['%10.0f','%10.4f','%10.4f','%10.4f'])

print('File saved with name' + fname )
