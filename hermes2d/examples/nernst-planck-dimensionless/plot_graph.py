# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
pylab.title("Error convergence")
pylab.xlabel("Physical time (s)")
pylab.ylabel("Error [%]")
data = numpy.loadtxt("time_error.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, "-s", label="error (est)")

legend()

# initialize new window
pylab.figure()

pylab.title("Problem size")
pylab.xlabel("Physical time (s)")
pylab.ylabel("Number of DOF")
data = numpy.loadtxt("time_dof.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, "-s", label="dof")

pylab.figure()

pylab.title("Time")
pylab.xlabel("Physical time (s)")
pylab.ylabel("CPU time")
data = numpy.loadtxt("time_cpu.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, "-s", label="cpu")
legend()

# finalize
show()
