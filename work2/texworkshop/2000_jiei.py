import numpy as np
import matplotlib.pyplot as plt
import copy

n=np.linspace(0,1000)
y = 2*(1-np.cos(n/1001*np.math.pi))

plt.plot(n,y)
plt.show()