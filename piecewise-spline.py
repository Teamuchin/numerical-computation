from matplotlib import pyplot as plt
import numpy as np

sp1x = np.linspace(0, 10, 100)
sp1y = -(1/165)*(sp1x-1)+(1/165)*(sp1x-1)**(3)

sp2x = np.linspace(0, 10, 100)
sp2y = (167/165)*(sp2x-2)+(1/55)*(sp2x-2)**(2)-(5/165)*(sp2x-2)**(3)

sp3x = np.linspace(0, 10, 100)
sp3y = (158/165)*(sp3x-3)-(4/55)*(sp3x-3)**(2)+(2/495)*(sp3x-3)**(3)+1
 
fig = plt.figure(figsize = (10, 10))
plt.plot(sp1x, sp1y,label = "spline 1")
plt.plot(sp2x, sp2y,label = "spline 2")
plt.plot(sp3x, sp3y,label = "spline 3")
plt.xticks(np.arange(0, 10, 1))
plt.yticks(np.arange(-5, 6, 1))
plt.grid()
plt.legend()
plt.show()
