from airfoil_discretisation import *
import solver as sl
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')

numPan = 50  # must be an even integer
foil = '0012'


'''initialise geometry and solver matrices'''
x, z = gen_airfoil(foil, numPan)
numPan, x, z, phi, S, x_control, z_control = discretization(x, z, numPan)
AN, AT = sl.matrix_building(numPan, x,z,x_control, z_control, S, phi)


'''Lift polar'''
alphas = np.arange(-17, 17.5, 0.5)
cls = np.zeros_like(alphas)
for i, a in enumerate(alphas):
    rhs = sl.rhs_vec(a, phi)
    Cp, gamma = sl.solver(AN, AT, rhs, phi, 0)
    cls[i] = sl.lift_calculator(Cp, x_control)

plt.suptitle('NACA 0012 Lift Curve')
plt.xlabel(r'$\alpha$ $[]$')
plt.ylabel(r'$C_{l}$')
plt.plot(alphas, cls)

# Read polar for 0012
exp_data = np.genfromtxt('0012.abbottdata.cl.dat', skip_header=6)
plt.scatter(exp_data[:,0],exp_data[:,1], c='b')



'''Plotting pressure ditributions at given α'''
# alpha = 4  # deg
# rhs = sl.rhs_vec(alpha, phi)
# Cp, gamma = sl.solver(AN, AT, rhs, phi, 0)
# Cl = sl.lift_calculator(Cp,x_control)
# fig= plt.figure()
# title = r'NACA 0008 airfoil; $\alpha = $' + str(alpha) + r'$^{\circ}$'
# plt.suptitle(title)
# plt.xlabel('x')
# plt.ylabel(r'$-C_{p}$')
# plt.plot(x_control, -Cp)
# plt.show()
