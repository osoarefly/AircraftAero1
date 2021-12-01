from airfoil_discretisation import *

x,y = load_airfoil('0008')

AoA = 0

theta = discretization(x,y, AoA)

print(theta)

