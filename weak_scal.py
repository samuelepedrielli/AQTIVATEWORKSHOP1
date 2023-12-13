import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
mpl.rc_context({'text.usetex' : True, 'font.family' : "Helvetica"})

X_train = [2,3,4,5,6,7,8,9,10] #number of processors, equal to the number of simulations
y_train = [] #put here the measured time
plt.scatter(X_train, y_train)
plt.plot(X_train, y_train)
plt.xlabel("$p$")
plt.ylabel("$T_{p}(N_{stat}=p)$")
_ = plt.title("Weak scalability: increase $p$ and $N_{stat}$ proportionally")
plt.savefig('/Users/samuelepedrielli/Desktop/Wscalab.png', dpi=300)
