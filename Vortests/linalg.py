import numpy as np
from time import *
import scipy.linalg

x = np.random.rand(1,1)
t1 = clock()
y = scipy.linalg.lu(x)
t2 = clock()-t1
print(t2)