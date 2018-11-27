import numpy as np

x = np.arange(-5, 5, 0.1)
y = np.linspace(0,1,10)
xx, yy = np.meshgrid(x, y)
print(xx)