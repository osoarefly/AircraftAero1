from airfoil_discretisation import *
import solver as sl
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

numPan = 50  # must be an even integer
alpha = 0  # deg
x, z = gen_airfoil('0008', numPan)
numPan, x, z, phi, S, x_control, z_control = discretization(x, z, numPan)
AN, AT = sl.matrix_building(numPan, x,z,x_control, z_control, S, phi)
rhs = sl.rhs_vec(alpha, phi)
Cp, gamma = sl.solver(AN, AT, rhs, phi, 0)

fig= plt.figure()
title = r'NACA 0008 airfoil; $\alpha = $' + str(alpha) + r'$^{\circ}$'
plt.suptitle(title)
plt.xlabel('x')
plt.ylabel(r'$-C_{p}$')
plt.plot(x_control, -Cp)
plt.show()





# t1, t2 = sl.unit_panel_induction((x[0],0,z[0]), (x[1],0,z[1]), (x_control[-1],0,z_control[-1]), 2)
# mats = sl.matrix_building(numPan, x, z, phi, S, x_control, z_control)
# rhs = sl.rhs_vec(0, 1, phi)

# gamma = np.linalg.solve(mats[0], rhs)
# nvels = np.matmul(mats[0], gamma)  # + np.sin(phi)
# tvels = np.matmul(mats[1], gamma) + np.cos(phi)
# tvels1, tvels2 = np.split(tvels, 2)

# plt.plot(1 - tvels1**2)
# plt.plot(np.flip(1 - tvels2**2))
# plt.show()