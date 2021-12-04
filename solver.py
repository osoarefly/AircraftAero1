import numpy as np
import matplotlib.pyplot as plt
import math as m

def matrix_building(numPan, ):

    # Initialize N (normal) and T (tangential) system matricesw
    N = np.zeros([numPan, numPan])
    T = np.zeros([numPan, numPan])

    # cycle through panels
    for i in range(numPan):
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



