import numpy as np
#test
def cartesian_to_spherical(x, y, z, vx, vy, vz):
    """
    Function to convert cartesian coordinates to spherical coordinates

    :param x: x coordinate
    :param y: y coordinate
    :param z: z coordinate
    :param vx: x velocity
    :param vy: y velocity
    :param vz: z velocity

    :return:
    r: spherical coordinates (array of shape Nx3 with r, phi, theta)
    v: spherical velocity components(array of shape Nx3 with vr, vphi, vtheta)

    """
    phi = np.arctan2(y, x)
    th = np.arccos(z / np.sqrt(x ** 2 + y ** 2 + z ** 2))
    r = np.sqrt(x ** 2. + y ** 2. + z ** 2.)

    # Code part from Mpetersen
    cost = z / r
    sint = np.sqrt(1. - cost * cost)
    cosp = np.cos(phi)
    sinp = np.sin(phi)

    vr = sint * (cosp * vx + sinp * vy) + cost * vz
    vphi = (-sinp * vx + cosp * vy)
    vtheta = (cost * (cosp * vx + sinp * vy) - sint * vz)

    return np.array([r, phi, th]), np.array([vr, vphi, vtheta])


def euler_xyz(phi, theta, psi=0., deg=False):
    """
    Function to compute the rotation matrix for a given set of euler angles

    :param phi: azimuthal rotation angle
    :param theta: inclination rotation angle
    :param psi: roll rotation angle
    :param deg: if True, angles are in degrees

    :return: Rotation matrix

    :note: The rotation matrix is defined as R = Rz(psi) * Ry(theta) * Rx(phi)
    :note: this code is adapted from the code by M. Petersen at https://github.com/michael-petersen/ReflexMotion/blob/master/reflexmotion/reflex.py
    """
    if deg == True:
        phi, theta = np.deg2rad(phi), np.rad2deg(theta)

    Rmatrix = np.array([[np.cos(theta) * np.cos(phi), \
                         np.cos(theta) * np.sin(phi), \
                         -np.sin(theta)], \
                        [np.sin(psi) * np.sin(theta) * np.cos(phi) - np.cos(psi) * np.sin(phi), \
                         np.sin(psi) * np.sin(theta) * np.sin(phi) + np.cos(psi) * np.cos(phi), \
                         np.cos(theta) * np.sin(psi)], \
                        [np.cos(psi) * np.sin(theta) * np.cos(phi) + np.sin(psi) * np.sin(phi), \
                         np.cos(psi) * np.sin(theta) * np.sin(phi) - np.sin(psi) * np.cos(phi), \
                         np.cos(theta) * np.cos(psi)]])
    return Rmatrix


def add_vtravel(vtravel, phi2):
    """
    Function to add the travel velocity to the spherical velocity components

    :param vtravel: travel velocity
    :param phi2: azimuthal angle in the rotated frame

    :return: array of shape Nx3 with vr, vphi, vtheta, the spherical velocity components in the rotated frame
    """

    vr = -vtravel * np.cos(phi2)
    vphi2 = +vtravel * np.sin(phi2)
    if len(phi2) >0:
        vphi1 = np.zeros_like(vr)
    else:
        vphi1 = 0.
    return np.array([vr, vphi1, vphi2])


def spherical_unit_vectors(phi, th):
    """
    Function to compute the spherical unit vectors

    :param phi: azimuthal angle
    :param th: polar angle

    :return: er, ephi, eth, the spherical unit vectors with shape (3, N) for each component
    """
    er = np.array([np.sin(th) * np.cos(phi), np.sin(th) * np.sin(phi), (np.cos(th))])
    ephi = np.array([-np.sin(phi), np.cos(phi), np.zeros_like(phi)])
    eth = np.array([np.cos(th) * np.cos(phi), np.cos(th) * np.sin(phi), -np.sin(th)])
    return er, ephi, eth


def spherical_to_cartesian(r, phi, th, vr, vphi, vth):
    """
    Function to convert spherical coordinates to cartesian coordinates

    :param r: radial distance
    :param phi: azimuthal angle
    :param th: polar angle
    :param vr: radial velocity
    :param vphi: azimuthal velocity
    :param vth: polar velocity

    :return:
    rcart: cartesian coordinates (array of shape Nx3 with x, y, z)
    vcart: cartesian velocity components(array of shape Nx3 with vx, vy, vz)
    """
    ur = np.array([np.sin(th) * np.cos(phi), np.sin(th) * np.sin(phi), (np.cos(th))])
    uphi = np.array([-np.sin(phi), np.cos(phi), np.zeros_like(r)])
    uth = np.array([np.cos(th) * np.cos(phi), np.cos(th) * np.sin(phi), -np.sin(th)])

    rcart = np.array(r * ur)

    vx = vr * ur[0] + vth * uth[0] + vphi * uphi[0]
    vy = vr * ur[1] + vth * uth[1] + vphi * uphi[1]
    vz = vr * ur[2] + vth * uth[2] + vphi * uphi[2]

    vcart = np.array([vx, vy, vz])

    return rcart, vcart