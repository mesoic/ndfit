import ndfit
import numpy as np
import matplotlib.pyplot as plt

x = list(np.linspace(-10,10,101))
y = [i*i*i for i in x]
z = ndfit.derivative(y,x)
print z

plt.plot(x,z)
plt.plot(x,y)
plt.show()
