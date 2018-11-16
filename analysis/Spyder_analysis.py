import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# make sure labels are large enough in report
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (10,8)
mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 'x-large'
mpl.rcParams['lines.linewidth'] = 2.5
mpl.rcParams['lines.markersize'] = 10.0


#temp_range = np.hstack([np.arange(0.5,2.,0.5), np.arange(2.,2.5,0.05), np.arange(2.5,5,0.5)])

L = 10

pos = np.loadtxt("../results/position.dat")
vel = np.loadtxt("../results/speed.dat")

#pos = np.random.rand(1000,100)
#vel = np.random.rand(1000,100)

fig = plt.figure()
#plt.legend()
ax = fig.add_subplot(111)
line, = ax.plot([], [], 'bo', ms=5)
ax.set_xlim(0, 1)
ax.set_ylim(-1,1)

def make_frame(t):
    X = pos[(t%1000),:]
    V = vel[(t%1000),:]
    line.set_data(X, V)
    return line,

ap = animation.FuncAnimation(fig, make_frame, interval=1,blit=True)
plt.show()
