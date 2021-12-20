from airfoil_discretisation import *
import solver as sl


numPan = 100 #must be an even integer
x,y = gen_airfoil('0008', int(numPan/2))
numPan, x, y, phi, S, x_control, y_control = discretization(x,y,0, numPan-2)
mats = sl.matrix_building(numPan, x, y, phi, S, x_control, y_control)
rhs = sl.rhs_vec(0,1,phi)

gamma = np.linalg.solve(mats[0],rhs)
tvels = np.matmul(mats[1],gamma) + np.cos(phi)
plt.plot(1-tvels**2)