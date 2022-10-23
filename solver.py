import numpy as np

def matrix_building(numPan, x, z, xc, zc, S, phi):
    '''
    Function to build an incidence matrix. Code adapted from Fortran Code presented by Chow & Keuthe in
    'Foundations of Aerodynamics: Bases of Aerodynamic Design', 5th Ed.; Page 161-163.

    :param numPan: Number of panels
    :param x: Array of panel boundary x coordinates.
    :param z: Array of panel boundary z coordinates.
    :param xc: Array of control point x coordinates.
    :param zc: Array of control point z coordinates.
    :param S: Array of panel lengths.
    :param phi: Array of panel orientation angles.

    :return AN: Normal velocity unit induction matrix
    :return AT: Tangential velocity unit induction matrix
    '''

    sin = np.sin(phi)
    cos = np.cos(phi)

    #Define empty matrices
    CN1 = np.zeros((numPan,numPan)) #Matrix mapping induced normal vel of boudary points j on each control point
    CN2 = np.zeros_like(CN1) #Matrix mapping induced normal vel of boudary points j+1 on each control point
    CT1 = np.zeros_like(CN1) #Matrix mapping induced tangential vel of boudary points j on each control point
    CT2 = np.zeros_like(CN1) #Matrix mapping induced tangential vel of boudary points j+1 on each control point

    for i in range(numPan):
        for j in range(numPan):
            if i==j: #if control point and boundary point lie on the same panel.
                CN1[i, j] = -1.
                CN2[i, j] = 1.
                CT1[i, j] = 0.5*np.pi
                CT2[i, j] = 0.5*np.pi
            else:
                #The following vars (A,B,C,D,E,F,G,P,Q) are used only to make the distributed vortex induction equation
                #simpler to enter programmatically. They do not necessarily represent anything on their own.
                A = -(xc[i] - x[j]) * cos[j] -(zc[i] - z[j]) * sin[j]
                B = (xc[i] - x[j])**2 + (zc[i] - z[j])**2
                C = np.sin(phi[i] - phi[j])
                D = np.cos(phi[i] - phi[j])
                E = (xc[i] - x[j])*sin[j] - (zc[i] - z[j])*cos[j]
                F = np.log(1. + S[j]* (S[j]+2*A)/B )
                G = np.arctan2(E*S[j], B+A*S[j])

                P = (xc[i] - x[j]) * np.sin(phi[i] - 2*phi[j]) + (zc[i] - z[j]) * np.cos(phi[i] - 2*phi[j])
                Q = (xc[i] - x[j]) * np.cos(phi[i] - 2*phi[j]) - (zc[i] - z[j]) * np.sin(phi[i] - 2*phi[j])

                CN2[i, j] = D + 0.5*Q*F/S[j] - (A*C+D*E)*G/S[j]
                CN1[i, j] = 0.5*D*F + C*G - CN2[i,j]
                CT2[i, j] = C + 0.5*P*F/S[j] + (A*D-C*E)*G/S[j]
                CT1[i, j] = 0.5*C*F - D*G - CT2[i,j]

    AN = np.zeros((numPan+1, numPan+1))
    AT = np.zeros_like(AN)

    for i in range(numPan):
        AN[i,0] = CN1[i,0]
        AN[i,numPan] = CN2[i,-1]
        AT[i,0] = CT1[i,0]
        AT[i,numPan] = CT2[i,-1]

        for j in range(1,numPan):
            AN[i,j] = CN1[i,j] + CN2[i,j-1]
            AT[i,j] = CT1[i,j] + CT2[i,j-1]

    AN[-1, 0] = 1.
    AN[-1, -1] = 1.

    return AN, AT

def rhs_vec(a, phi):
    '''
    Function to initialise a right-hand-side vector
    :param a: Angle of attack [deg]
    :param phi: Panel orientation angle [rad]
    :return vec: right hand side vector of size numPan+1
    '''

    vec = np.sin(phi - a * np.pi/180)
    print(vec.shape)
    vec = np.append(vec, 0)
    print(vec.shape)
    return vec

def solver(AN, AT, rhs, phi, alpha):
    gammas = np.linalg.solve(AN,rhs)
    Vt = np.cos(phi-alpha) + np.matmul(AT, gammas)[:-1]
    Cp = 1 - Vt**2
    return Cp, gammas

# if __name__ == "__main__":



# def induced_velocity(xcontrol, ycontrol, ):
