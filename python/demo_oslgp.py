################################################################
## Python port of https://github.com/markusheinonen/adaptivegp
## (commit ed9ac09, Oct 18 2019) by @dereksoeder
##
## demo_oslgp.py
## Last updated July 25, 2024
################################################################

from scipy.io import loadmat
from adaptivegp import nsgp, plotnsgp, plotnsgpsamples, MLfigure, plt_show

data = loadmat('../data/datasets.mat')

# fit fully ns-gp
lsogp,_,_,_,_,_ = nsgp(data['Dlso']['x'].item(), data['Dlso']['y'].item(), 'lso', 'grad')
MLfigure()
plotnsgp(lsogp,True)



# fit with HMC and plot?(takes a while)
lgp,samples,_,_,_,_ = nsgp(data['Dl']['x'].item(), data['Dl']['y'].item(), 'l', 'hmc')
MLfigure()
plotnsgpsamples(lgp,samples,True)

plt_show()
