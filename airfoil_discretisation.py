import numpy as np
import matplotlib.pyplot as plt
import math as m

def load_airfoil(filename): #filename is airfoil number as STRING
    loc = 'airfoils/'+filename+'.dat'
    data = np.genfromtxt(fname=loc,skip_header=1)
    x = data[:,0]
    y = data[:,1]
    return x, y

def gen_airfoil(airfoilname, npoints):
    b = np.linspace(0,np.pi,npoints)
    x = 0.5*(1-np.cos(b))
    t = float(airfoilname[-2:])/100
    m = float(airfoilname[0])/100
    p = float(airfoilname[1])/10
    print(m,p,t)




    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.126 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)

    if m==0.:
        yu = yt
        yl = -1*np.flip(yt)
        xu = x
        xl = np.flip(x)
        xf = np.append(xl,xu)
        yf = np.append(yl,yu)
    else:
        x1 = x[x <= p]
        x2 = x[x > p]
        yc1 = m / p ** 2 * (2 * p * x1 - x ** 2)
        yc2 = m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x ** 2)
        yc = np.append(yc1, yc2)
        dydx1 = 2*m/p**2 * (p-x1)
        dydx2 = 2*n/(1-p)**2 * (p-x2)
        dydx = np.append(dydx1,dydx2)
        theta = np.atan(dydx)

        xu = x - yt*np.sin(theta)
        yu = yc + yt*np.cos(theta)
        xl = x + yt*np.sin(theta)
        yl = yc - yt*np.cos(theta)

        xl = np.flip(xl)
        yl = np.flip(yl)

        xf = np.append(xl,xu)
        yf = np.append(yl,yu)

    return xf, yf

def discretization(x,y,AoA):

    # Number of panels
    numPan = len(x) - 1   # Number of panels

    # Check for direction of points
    edge = np.zeros(numPan)  # Initialize edge check value
    for i in range(numPan):  # Loop over all panels
        edge[i] = (x[i + 1] - x[i]) * (y[i + 1] - y[i])  # Compute edge value for each panel

    sumEdge = np.sum(edge)  # Sum all panel edge values

    # If panels are CCW, flip them (don't if CW)
    if (sumEdge < 0):  # If sum is negative
        print('Points are counter-clockwise.  Flipping.\n')  # Display message in console
        x = np.flipud(x)  # Flip the X boundary points array
        y = np.flipud(y)  # Flip the Y boundary points array
    elif (sumEdge > 0):  # If sum is positive
        print('Points are clockwise.  Not flipping.\n')  # Do nothing, display message in consolve

    # %% COMPUTE GEOMETRIC VARIABLES

    # Initialize variables
    x_control = np.zeros(numPan)  # Initialize X control points
    y_control = np.zeros(numPan)  # Initialize Y control points
    S = np.zeros(numPan)  # Initialize panel lengths
    phi = np.zeros(numPan)  # Initialize panel orientation angles

    # Find geometric quantities of the airfoil
    for i in range(numPan):  # Loop over all panels
        x_control[i] = 0.5 * (x[i] + x[i + 1])  # X control point coordinate
        y_control[i] = 0.5 * (y[i] + y[i + 1])  # Y control point coordinate
        dx = x[i + 1] - x[i]  # Panel X length
        dy = y[i + 1] - y[i]  # Panel Y length
        S[i] = (dx ** 2 + dy ** 2) ** 0.5  # Panel length

        phi[i] = m.atan2(dy, dx)  # Panel orientation angle [rad]
        if (phi[i] < 0):  # If panel orientation is negative
            phi[i] = phi[i] + 2 * np.pi  # Add 2pi to the panel angle

    # Compute angle of panel normal w.r.t. horizontal and include AoA
    delta = phi + (np.pi / 2)  # Compute panel normal angle [rad]
    beta = delta - (AoA * (np.pi / 180))  # Angle between freestream and panel normal [rad]

    # %% PLOTTING
    # Plot the paneled geometry
    fig = plt.figure(1)  # Create figure
    plt.cla()  # Get ready for plotting
    plt.fill(x, y, 'k')  # Plot polygon (circle or airfoil)
    X = np.zeros(2)  # Initialize panel X variable
    Y = np.zeros(2)  # Initialize panel Y variable
    for i in range(numPan):  # Loop over all panels
        X[0] = x_control[i]  # Panel starting X point
        X[1] = x_control[i] + S[i] * np.cos(delta[i])  # Panel ending X point
        Y[0] = y_control[i]  # Panel starting Y point
        Y[1] = y_control[i] + S[i] * np.sin(delta[i])  # Panel ending Y point
        if (i == 0):  # For first panel
            plt.plot(X, Y, 'b-', label='First Panel')  # Plot the first panel normal vector
        elif (i == 1):  # For second panel
            plt.plot(X, Y, 'g-', label='Second Panel')  # Plot the second panel normal vector
        else:  # For every other panel
            plt.plot(X, Y, 'r-')  # Plot the panel normal vector
    plt.xlabel('X-Axis')  # Set X-label
    plt.ylabel('Y-Axis')  # Set Y-label
    plt.title('Panel Geometry')  # Set title
    plt.axis('equal')  # Set axes equal
    plt.legend()  # Plot legend
    plt.show()  # Display plot

    return numPan, x, y, phi, S, x_control, y_control

x,y = gen_airfoil('0008',100)
# plt.plot(gen_airfoil('0008',100))
# plt.show()