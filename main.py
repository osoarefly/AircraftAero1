from airfoil_discretisation import *

numpan = 100
x,y = gen_airfoil('2412', numpan)

AoA = 0

theta = discretization(x,y, AoA, 2*numpan-1)

print(theta)

