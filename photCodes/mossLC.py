import numpy as np, matplotlib.pyplot as plt, glob

#-----------------
def plotLC(fname, shift=0, Format = 'k.'):
    """ Plot LC """
    data = np.loadtxt(fname)
    data[:,1] += shift
    plt.errorbar(data[:,0], data[:,1], data[:,2], fmt=Format, ecolor='gray',
    lw=1, ms=4, capsize=1.5, label=fname[-14:-4])

#-----------------
folder = '/media/marin/TRAPPIST3/troppist/pho/348Moss/'

"""
night = '20171213'
night = '20171214b2'

plotLC(folder+night+'.txt')
plt.legend()
plt.gca().invert_yaxis()
plt.show()
"""

pattern = '201712*'
fnames = sorted(glob.glob(folder + pattern))

col=['C0.','C1.','C2.','C3.','C4.','C5.','C6.','C7.','C8.','C9.',
     'C0.','C1.','C2.','C3.','C4.','C5.','C6.','C7.','C8.','C9.']
for i in range(len(fnames)):
    fname = fnames[i]
    plotLC(fname, Format=col[i])
    plt.legend()
    plt.gca().invert_yaxis()
plt.show()

#-------------------

night1 = '20171213'
night2 = '20171213b'

plotLC(folder+night1+'.txt')
plotLC(folder+night2+'.txt', -0.12)
plt.gca().invert_yaxis()
plt.legend()
plt.show()
