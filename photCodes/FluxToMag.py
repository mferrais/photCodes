import numpy as np

data = np.loadtxt('/media/marin/TRAPPIST3/troppist/pho/348_20171214/348_TN_20171214_13_13.mag')

mag = data[:,1]
flux  =  pow(10, -mag/2.5)
flux /= np.mean(flux)
data[:,1] = flux
np.savetxt('/media/marin/TRAPPIST3/troppist/pho/348_20171214/348_TN_20171214_13_13.fluxTest',
            data, fmt=['%10.6f','%.6e','%.2e','%10.0f'])
