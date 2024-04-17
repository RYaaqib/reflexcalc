import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab as plt
import matplotlib.cm as cm

from coord import *
"""
This script is used to generate the reflex motion for halo stars in the MW. The best fit values are the same as those foundin
Yaaqib, Petersen and Penarrubia 2024. The reflex motion is generated for 4 different distances, 20-30, 30-40, 40-50 and 50+ kpc.
All values used here are generated from table 1 in the paper. The reflex motion is then plotted on the sky in Mollweide projection.

In order to use the file, simply provide a file with
"""
vsun_mw = np.array([11.1, 244.24, 7.25])  # km/s
rsun_mw = np.array([-8.3, 0., 0.02])  # kpc


def read_posterior(infile, raw=False, cosb=True):
    if raw == True:
        I = np.genfromtxt((infile + '.txt'))
        I = I[:, 2:]
    else:
        I = np.genfromtxt((infile + 'post_equal_weights.dat'))

    P = dict()
    P['l'] = np.rad2deg(I[:, 0])
    if cosb == True:
        P['b'] = 90 - np.rad2deg(np.arccos(I[:, 1]))
    else:
        P['b'] = 90 - np.rad2deg(I[:, 1])
    P['vtravel'] = I[:, 2]

    P['vr'] = I[:, 3]
    P['vphi'] = I[:, 4]
    P['vth'] = I[:, 5]

    P['sigvlos'] = np.sqrt(1. / I[:, 6])
    P['sigmul'] = np.sqrt(1. / I[:, 7])
    P['sigmub'] = np.sqrt(1. / I[:, 8])

    return P

# this is a slightly modified version of the get_v function in the reflex code
def get_v(cube, rgal, vgal, solar=False):
    # 2. compute spherical galactic coordinates
    rgalsph = np.zeros_like(rgal)
    # vgalsph = np.zeros((len(x), 3))
    tmp = cartesian_to_spherical(rgal[:, 0], rgal[:, 1], rgal[:, 2],
                                 vgal[:, 0], vgal[:, 1], vgal[:, 2])
    rgalsph[:, :] = tmp[0].T

    # 3. rotate cartesian coordinates such that the z axis points at lapex, bapex
    # through the euler angle rotation x-y-z
    bapex = cube[1]
    Rrot = euler_xyz(cube[0], bapex, deg=False)
    rp = np.zeros_like(rgal)
    vp = np.zeros_like(vgal)

    rp = np.dot(Rrot, rgal.T).T

    # 4. compute the spherical coordinates of r,v in rotated frame

    rpsph = np.zeros_like(rp)
    vpsph = np.zeros_like(vp)

    tmp = cartesian_to_spherical(rp[:, 0], rp[:, 1], rp[:, 2],
                                 vp[:, 0], vp[:, 1], vp[:, 2])
    rpsph[:, :] = tmp[0].T
    vpsph[:, :] = add_vtravel(cube[2], rpsph[:, 2]).T  # + tmp[1]
    # 5. translate the spherical coordinates in the rotated frame back to cartesian coordinates
    # note: we don't care about the position vector anymore, we only want the velocities

    rp1 = np.zeros_like(rp)
    vp1 = np.zeros_like(rp)

    temp = spherical_to_cartesian(rpsph[:, 0], rpsph[:, 1], rpsph[:, 2],
                                  vpsph[:, 0], vpsph[:, 1], vpsph[:, 2])
    rp1[:, :] = temp[0].T
    vp1[:, :] = temp[1].T
    # 6. now rotate the coordinates back to the galactic frame

    rgal1 = np.zeros_like(rgal)
    vgal1 = np.zeros_like(rgal)


    # for i in range(len(x)):
    rgal1[:, :] = np.dot(np.linalg.inv(Rrot), rp1[:, :].T).T
    vgal1[:, :] = np.dot(np.linalg.inv(Rrot), vp1[:, :].T).T
    # 7. compute the vector v and translate to be at the sun
    r = np.zeros_like(rp)
    v = np.zeros_like(vp)

    # 7a. convert bulk velocities to cartesian velocities
    tmp = spherical_to_cartesian(rgalsph[:, 0], rgalsph[:, 1],
                                 rgalsph[:, 2], cube[3], cube[4], cube[5])

    if solar == True:
        r[:, :] = rgal1[:, :] - rsun_mw
        v[:, :] = vgal1[:, :] + tmp[1].T - vsun_mw
    else:
        r[:, :] = rgal1[:, :] - rsun_mw
        v[:, :] = vgal1[:, :] + tmp[1].T
    # 8. find the galactic proper motions using the unit vectors of ul, ub, ulos
    # spherical unit vectors in the galactic frame..
    rsunsph = np.zeros_like(r)
    vsunsph = np.zeros_like(v)


    #vectorized
    temp = cartesian_to_spherical(r[:, 0], r[:, 1], r[:, 2],
                                  v[:, 0], v[:, 1], v[:, 2])
    rsunsph[:, :] = temp[0].T
    vsunsph[:, :] = temp[1].T
    # now multiply the vector v by the unit vectors in spherical coordinates centred at
    # the positions of the sun

    elos, ephi, eth = spherical_unit_vectors(rsunsph[:, 1], rsunsph[:, 2])
    elos, ephi, eth = elos.T, ephi.T, eth.T

    # helicocentric distance
    dist = np.linalg.norm(rsunsph, axis=1)
    fac = 4.74057 * dist

    dsun = rsunsph[:,0]
    l    = rsunsph[:,1]
    b    = rsunsph[:,2]
    vlos = np.einsum('ij,ij->i', v, elos)
    mul = np.einsum('ij,ij->i', v, ephi) * fac
    mub = -np.einsum('ij,ij->i', v, eth) * fac
    return l,b, dsun, vlos, mul, mub


def plot_reflex_model(cube, rgal, vgal, ax=None, quant="vlos", vlosnorm=None,mulnorm=None,mubnorm=None):

    #lapex, bapex, vtravel, vr, vphi, vtheta,siglos, sigl, sigb

    l,b,dist, vlos, mul, mub = get_v(cube,rgal*300., vgal*240./1.4, solar=False)

    if ax == None:
        fig, ax = plt.subplots(3, facecolor="white", figsize=(4, 8), subplot_kw={'projection': 'mollweide'})
    else:
        print("using axes defined outside function")
    # fig.patch.set_facecolor('black') #setting plot background to dark colour

    l[l>np.pi]-=2.*np.pi
    b = np.pi/2. - b

    if quant == "vlos":
        if vlosnorm is not None:
            print("using vlosnorm")
            ax.scatter(-l, b, c=cm.coolwarm((vlos -np.min(vlosnorm))/(np.max(vlosnorm) - np.min(vlosnorm))), cmap="seismic")
        else:
            ax.scatter(-l, b, c=cm.coolwarm((vlos -np.min(vlos))/(np.max(vlos) - np.min(vlos))), cmap="seismic")
    elif quant == "mul":
        if mulnorm is not None:
            ax.scatter(-l, b, c=cm.coolwarm((mul - np.min(mulnorm))/(np.max(mulnorm)-np.min(mulnorm))), cmap="seismic")
        else:
            ax.scatter(-l, b, c=cm.coolwarm((mul - np.min(mul))/(np.max(mul)-np.min(mul))), cmap="seismic")
    elif quant == "mub":
        if mubnorm is not None:
            ax.scatter(-l, b, c=cm.coolwarm((mul - np.min(mubnorm))/(np.max(mubnorm)-np.min(mubnorm))), cmap="seismic")
        else:
            ax.scatter(-l, b, c=cm.coolwarm((mub - np.min(mub))/(np.max(mub)-np.min(mub)), cmap="seismic"))

    return ax
