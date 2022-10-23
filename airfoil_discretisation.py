import numpy as np
import matplotlib.pyplot as plt

def load_airfoil(filename): #filename is airfoil number as STRING
    loc = 'airfoils/'+filename+'.dat'
    data = np.genfromtxt(fname=loc,skip_header=1)
    x = data[:,0]
    z = data[:,1]
    return x, z #checked

def gen_airfoil(airfoilname, numPan):
    '''
    Function to generate a NACA 4-digit airfoil with a given number of panels.

    :param airfoilname: Name of airfoil as string.  Must be 4 digits only.
    :param numPan: Number of panels. First and last panel will have point j and j+1 as the same.

    :return xf: x coordinate of each panel boundary point
    :return zf: z coordinate of each panel boundary point
    '''

    #airfoil name must be 4-digits as STRING
    npoints = numPan/2+1

    #generate cosine distribution of x locations on airfoil chordline (LE and TE clustering)
    x = 0.5*(1-np.cos(np.linspace(0, np.pi, int(npoints))))

    t = float(airfoilname[-2:])/100 #max thickness
    m = float(airfoilname[0])/100 #maximum camber
    p = float(airfoilname[1])/10 #Max camber location
    print(m,p,t)

    zt = 5 * t * (0.2969 * np.sqrt(x) - 0.126 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4) #thickness function for NACA 4digit airfoils

    if m==0.: #if symmetric
        zu = zt
        zl = -1*np.flip(zt)
        xu = x
        xl = np.flip(x)
        xf = np.append(xl,xu[1:])
        zf = np.append(zl,zu[1:])
    else: #if asymmetric
        x1 = x[x <= p] #all x values on LE side of max camber location
        x2 = x[x > p] #all x values on TE side of max camber location

        #from NACA 4-series definition
        zc1 = m / p ** 2 * (2 * p * x1 - x1 ** 2)
        zc2 = m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x2 - x2 ** 2)
        zc = np.append(zc1, zc2)
        dzdx1 = 2*m/p**2 * (p-x1)
        dzdx2 = 2*m/(1-p)**2 * (p-x2)
        dzdx = np.append(dzdx1,dzdx2)
        theta = np.arctan(dzdx)

        xu = x - zt*np.sin(theta)
        zu = zc + zt*np.cos(theta)
        xl = x + zt*np.sin(theta)
        zl = zc - zt*np.cos(theta)

        xl = np.flip(xl)
        zl = np.flip(zl)

        xf = np.append(xl,xu[1:])
        zf = np.append(zl,zu[1:])

    return xf, zf #checked

def discretization(x,z, numPan):
    '''
    Function to use airfoil coordinates to calculate positions of control points and panel characteristic values.

    :param x: Panel boundary x coordinate array
    :param z: Panel boundary z coordinate array
    :param numPan: Number of panels

    :return numPan: Number of panels. Unchanged from input value.
    :return x: Panel boundary x coordinate array. Unchanged from input value.
    :return z: Panel boundary z coordinate array. Unchanged from input value.
    :return phi: Panel orientation angle. Defined as the angle measured counter-clockwise from +ve x-axis to the inside
    of panel.
    :return S: Panel length.
    :return x_control: Array of panel control point x coordinate.
    :return z_control: Array of panel control point z coordinate.
    '''

    # Check for direction of points
    edge = np.zeros(numPan)  # Initialize edge check value
    for i in range(numPan):  # Loop over all panels
        edge[i] = (x[i + 1] - x[i]) * (z[i + 1] - z[i])  # Compute edge value for each panel

    # Initialize variables
    x_control_up = np.zeros(int(numPan/2))  # Initialize X control points
    z_control_up = np.zeros(int(numPan/2))  # Initialize Z control points
    x_control_lo = np.zeros(int(numPan / 2))  # Initialize X control points
    z_control_lo = np.zeros(int(numPan / 2))
    S = np.zeros(numPan)  # Initialize panel lengths
    phi = np.zeros(numPan)  # Initialize panel orientation angles

    # Find geometric quantities of the airfoil
    for i in range(int(numPan/2)):  # Loop over all panels

        # Panel midpoint coordinates for upper and lower surface
        x_control_lo[i] = 0.5 * (x[i] + x[i + 1])
        z_control_lo[i] = 0.5 * (z[i] + z[i + 1])
        x_control_up[i] = 0.5 * (x[-i-1] + x[-i-2])
        z_control_up[i] = 0.5 * (z[-i-1] + z[-i-2])

        # Panel 75% chord coordinates for upper and lower surface
        # x_control_lo[i] = 0.5 * (x_half_lo + x[i])
        # z_control_lo[i] = 0.5 * (z_half_lo + z[i])
        # x_control_up[i] = 0.5 * (x_half_up + x[-i-1])
        # z_control_up[i] = 0.5 * (z_half_up + z[-i-1])

    for i in range(numPan):
        dx = x[i + 1] - x[i]  # Panel X length
        dz = z[i + 1] - z[i]  # Panel Z length
        S[i] = (dx * dx + dz * dz) ** 0.5  # Panel length
        phi[i] = np.arctan2(dz, dx)

        if phi[i] < 0:
            phi[i] = phi[i] + 2*np.pi

    x_control = np.append(x_control_lo,np.flip(x_control_up))
    z_control = np.append(z_control_lo,np.flip(z_control_up))
    # Compute angle of panel normal w.r.t. horizontal and include AoA
    # delta = phi + (np.pi / 2)  # Compute panel normal angle [rad]
    # beta = delta - (AoA * (np.pi / 180))  # Angle between freestream and panel normal [rad]
    z = np.round(z, 10)
    return numPan, x, z, phi ,S, x_control, z_control

# x,y = gen_airfoil('4412',100)
# plt.plot(x,y)
# plt.axis('equal')
# plt.show()