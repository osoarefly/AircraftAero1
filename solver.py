import numpy as np
import matplotlib.pyplot as plt
import math as m

def matrix_building(numPan, X_panel, Y_panel, theta, S, x_control, y_control):

    # Initialize N (normal) and T (tangential) system matrices
    N = np.zeros([numPan, numPan])
    T = np.zeros([numPan, numPan])

    # cycle through panels
    for i in range(numPan-1):

        for j in range(numPan):

            A = - (x_control[i] - X_panel[j]) * np.cos(theta[i]) - (y_control[i] - Y_panel[j]) * np.sin(theta[j])
            B = (x_control[i] - X_panel[j])**2 - (y_control[i] - Y_panel[j])**2
            C = np.sin(theta[i] - theta[j])
            D = np.cos(theta[i] - theta[j])
            E = (x_control[i] - X_panel[j]) * np.sin(theta[j]) - (y_control[i] - Y_panel[j]) * np.cos(theta[j])
            F = np.log(1 + (S[j]**2 * A * S[j])/B )
            G = np.arctan((E * S[j])/(B + A * S[j]))
            P = (x_control[i] - X_panel[j]) * np.sin(theta[i] - 2*theta[j]) + (y_control[i] - Y_panel[j]) * np.cos(theta[i] - 2*theta[j])
            Q = (x_control[i] - X_panel[j]) * np.cos(theta[i] - 2*theta[j]) + (y_control[i] - Y_panel[j]) * np.cos(theta[i] - 2*theta[j])

            Cn2 = D + 0.5*Q*F/S[j] - (A*C + D*E)*G/S[j]
            Cn1 = 0.5*D*F + C*G - Cn2
            if j==numPan-1:
                N[i,j] = N[i,j] + Cn1
                N[i,0] = N[i,0] + Cn2
            else:
                N[i,j] = N[i,j] + Cn1
                N[i,j+1] = N[i,j+1] + Cn2



    N[-1,0] = 1 #kutta condition
    N[-1,-1] = 1 #kutta condition

    return N, T

def rhs_vec(a, Vinf, theta):
    rhs = Vinf*np.sin(theta-a)
    np.append(rhs,0.)
    return rhs

# if __name__ == "__main__":



# def induced_velocity(xcontrol, ycontrol, ):
