import numpy as np
from matplotlib import pyplot as plt

#putting (x,y) values to ax+b=y equation to create matrixes

m1x = np.matrix([[-2,1],[0,1],[3,1],[5,1]])
m1y = np.matrix([[3],[1],[0],[2]])
m1xt = np.matrix.getT(m1x)

#xtranspose*x*solution = xtranspose*y

least_sq = np.linalg.solve(np.dot(m1xt,m1x),np.dot(m1xt,m1y))
print(least_sq[0,0])

x1 = np.linspace(-5, 10, 1000)
y1 = least_sq[0,0]*x1+least_sq[1,0]
plt.figure()
plt.plot(x1, y1)

print("\n")
#putting (x,y) values to ax^2+bx+c=y equation to create matrixes

m2x = np.matrix([[4,-2,1],[0,0,1],[9,3,1],[25,5,1]])
m2y = np.matrix([[3],[1],[0],[2]])
m2xt = np.matrix.getT(m2x)

#xtranspose*x*solution = xtranspose*y

least_sq2 = np.linalg.solve(np.dot(m2xt,m2x),np.dot(m2xt,m2y))
print(least_sq2)

x2 = np.linspace(-5, 10, 1000)
y2 = least_sq2[0,0]*(x2**2)+least_sq2[1,0]*x2+least_sq2[2,0]
plt.figure()
plt.plot(x2, y2)

print("\n")

#root mean squared error for m1

print(np.linalg.norm(np.subtract(m1y,np.dot(m1x,least_sq))))

print("\n")

#root mean squared error for m2

print(np.linalg.norm(np.subtract(m2y,np.dot(m2x,least_sq2))))

print("\n")

#reduced qr factorizaton for m1

#first values gathered for required calculations

m3 = m1x

m3_1 = np.matrix([[-2],[0],[3],[5]])

m3_2 = np.matrix([[1],[1],[1],[1]])

m3_2t = np.matrix.getT(m3_2)

m3_u1 = m3_1

m3_u1t = np.matrix.getT(m3_u1)

#u2 = a2-((transpose of a2*u1)/(transpose of u1*u1))

m3_u2 = np.subtract(m3_2,np.dot((np.dot(m3_2t,m3_u1)/np.dot(m3_u1t,m3_u1))[0,0],m3_u1))

#q1 = u1/length of vector u1              q2 = u2/length of u2

m3_q1 = np.divide(m3_u1,np.linalg.norm(m3_u1))

m3_q2 = np.divide(m3_u2,np.linalg.norm(m3_u2))

#merge q1 and q2 into one q

m3_q = np.matrix.getT(np.concatenate((np.matrix.getT(m3_q1),np.matrix.getT(m3_q2))))

#r = transpose of q*a

m3_r = np.dot(np.matrix.getT(m3_q),m3)

print(m3_q)

print(m3_r)

#solving x in R*x=transpose of q*y

qr_lsq1 = np.linalg.solve(m3_r, np.dot(np.matrix.getT(m3_q),m1y))

x3 = np.linspace(-5, 10, 1000)
y3 = qr_lsq1[0,0]*x3+qr_lsq1[1,0]
plt.figure()
plt.plot(x3, y3)



m4 = m2x

m4_1 = np.matrix([[4],[0],[9],[25]])

m4_2 = np.matrix([[-2],[0],[3],[5]])

m4_3 = np.matrix([[1],[1],[1],[1]])

m4_2t = np.matrix.getT(m4_2)

m4_3t = np.matrix.getT(m4_3)

m4_u1 = m4_1

m4_q1 = np.divide(m4_u1,np.linalg.norm(m4_u1))

m4_q1t = np.matrix.getT(m4_q1)

m4_u2 = np.subtract(m4_2,np.dot((np.dot(m4_2t,m4_u1)/np.dot(m4_q1t,m4_u1))[0,0],m4_u1))

m4_q2 = np.divide(m4_u2,np.linalg.norm(m4_u2))

m4_q2t = np.matrix.getT(m4_q2)

m4_u31 = np.dot((np.dot(m4_3t,m4_u1)/np.dot(m4_q1t,m4_u1))[0,0],m4_u1)

m4_u32 = np.dot((np.dot(m4_3t,m4_u2)/np.dot(m4_q2t,m4_u2))[0,0],m4_u2)

m4_u3 = np.subtract(m4_2,np.add(m4_u32,m4_u31))

m4_q3 = np.divide(m4_u3,np.linalg.norm(m4_u3))

m4_q = np.matrix.getT(np.concatenate(((np.concatenate((np.matrix.getT(m4_q1),np.matrix.getT(m4_q2)))),np.matrix.getT(m4_q3))))

m4_r = np.dot(np.matrix.getT(m4_q),m4)

print(m4_q)

print(m4_r)

qr_lsq2 = np.linalg.solve(m4_r, np.dot(np.matrix.getT(m4_q),m2y))

x4 = np.linspace(-5, 5, 1000)
y4 = qr_lsq2[0,0]*(x4**2)+qr_lsq2[1,0]*x4+qr_lsq2[2,0]
plt.figure()
plt.plot(x4, y4)

print("\n")



plt.show()