import numpy
import math
from circle import circle

# xrange just "range" in python3.
# This code means fastest implementation used in 2 and 3
try:
    xrange
except NameError:
    xrange = range

def phaseFromZernikes(zCoeffs, size, norm="noll"):
    """
    Creates an array of the sum of zernike polynomials with specified coefficeints

    Parameters:
        zCoeffs (list): zernike Coefficients
        size (int): Diameter of returned array
        norm (string, optional): The normalisation of Zernike modes. Can be ``"noll"``, ``"p2v"`` (peak to valley), or ``"rms"``. default is ``"noll"``.

    Returns:
        ndarray: a `size` x `size` array of summed Zernike polynomials
    """
    Zs = zernikeArray(len(zCoeffs), size, norm=norm)
    phase = numpy.zeros((size, size))
    for z in xrange(len(zCoeffs)):
        phase += Zs[z] * zCoeffs[z]

    return phase

def zernike(j, N):
    """
     Creates the Zernike polynomial with mode index j,
     where j = 1 corresponds to piston.

     Args:
        j (int): The noll j number of the zernike mode
        N (int): The diameter of the zernike more in pixels
     Returns:
        ndarray: The Zernike mode
     """

    n, m = zernIndex(j)
    return zernike_nm(n, m, N)

def zernike_nm(n, m, N):
    """
     Creates the Zernike polynomial with radial index, n, and azimuthal index, m.

     Args:
        j (int): The noll j number of the zernike mode
        N (int): The diameter of the zernike more in pixels
     Returns:
        ndarray: The Zernike mode
     """
    coords = numpy.linspace(-1, 1, N)
    X,Y = numpy.meshgrid(coords, coords)
    R = numpy.sqrt(X**2 + Y**2)
    theta = numpy.arctan2(Y, X)

    if m==0:
        Z = numpy.sqrt(n+1)*zernikeRadialFunc(n, 0, R)
    else:
        if m > 0: # j is even
            Z = numpy.sqrt(2*(n+1)) * zernikeRadialFunc(n, m, R) * numpy.cos(m*theta)
        else:   #i is odd
            m = abs(m)
            Z = numpy.sqrt(2*(n+1)) * zernikeRadialFunc(n, m, R) * numpy.sin(m * theta)


    return Z*circle(N/2., N)



def zernikeRadialFunc(n, m, r):
    """
    Fucntion to calculate the Zernike radial function

    Parameters:
        n (int): Zernike radial order
        m (int): Zernike azimuthal order
        r (ndarray): 2-d array of radii from the centre the array

    Returns:
        ndarray: The Zernike radial function
    """

    R = numpy.zeros(r.shape)
    for i in xrange(0,int((n-m)/2)+1):

        R += r**(n-2*i) * (((-1)**(i))*numpy.math.factorial(n-i)) / ( numpy.math.factorial(i) * numpy.math.factorial(0.5*(n+m)-i) * numpy.math.factorial(0.5*(n-m)-i) )

    return R



def zernIndex(j):
    """
    Find the [n,m] list giving the radial order n and azimuthal order
    of the Zernike polynomial of Noll index j.

    Parameters:
        j (int): The Noll index for Zernike polynomials

    Returns:
        list: n, m values
    """
    n = int((-1 + math.sqrt(8*(j-1)+1)) /2.)
    p = (j-(n*(n+1))/2.)
    k = n%2
    m = int((p+k)/2.)*2 - k

    if m!=0:
        if j%2==0:
            s=1
        else:
            s=-1
        m *= s

    return [n, m]


def zernikeArray(J, N, norm="noll"):
    """
    Creates an array of Zernike Polynomials

    Parameters:
        maxJ (int or list): Max Zernike polynomial to create, or list of zernikes J indices to create
        N (int): size of created arrays
        norm (string, optional): The normalisation of Zernike modes. Can be ``"noll"``, ``"p2v"`` (peak to valley), or ``"rms"``. default is ``"noll"``.

    Returns:
        ndarray: array of Zernike Polynomials
    """
    # If list, make those Zernikes
    try:
        nJ = len(J)
        Zs = numpy.empty((nJ, N, N))
        for i in xrange(nJ):
            Zs[i] = zernike(J[i], N)

    # Else, cast to int and create up to that number
    except TypeError:

        maxJ = int(numpy.round(J))
        N = int(numpy.round(N))

        Zs = numpy.empty((maxJ, N, N))

        for j in xrange(1, maxJ+1):
            Zs[j-1] = zernike(j, N)


    if norm=="p2v":
        for z in xrange(len(Zs)):
            Zs[z] /= (Zs[z].max()-Zs[z].min())

    elif norm=="rms":
        for z in xrange(len(Zs)):
            # Norm by RMS. Remember only to include circle elements in mean
            Zs[z] /= numpy.sqrt(
                    numpy.sum(Zs[z]**2)/numpy.sum(circle(N/2., N)))

    return  Zs

