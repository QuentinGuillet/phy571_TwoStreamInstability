import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
'''
# make sure labels are large enough in report
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (10,8)
mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 'x-large'
mpl.rcParams['lines.linewidth'] = 2.5
mpl.rcParams['lines.markersize'] = 10.0


#temp_range = np.hstack([np.arange(0.5,2.,0.5), np.arange(2.,2.5,0.05), np.arange(2.5,5,0.5)])
'''
L = 10

pos = np.loadtxt("../results/position.dat")
vel = np.loadtxt("../results/speed.dat")

Writer = animation.writers['ffmpeg']
writer = Writer(fps=30, bitrate=1800)
#pos = np.random.rand(1000,100)
#vel = np.random.rand(1000,100)

fig = plt.figure()
#plt.legend()
ax = plt.axes(xlim=(0, 1000), ylim=(0,1000))
line, = ax.plot([], [], 'bo', ms=5)

def make_frame(t):
    X = t
    V = t
    line.set_data(X, V)
    return line,

ap = animation.FuncAnimation(fig, make_frame, frames = 1000,interval=33,blit=True)
#plt.show()
ap.save("ouioui.mp4",writer=writer)
