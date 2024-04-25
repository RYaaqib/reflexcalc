import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab as plt
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator
from coord import *
"""
This script is used to generate the reflex motion for halo stars in the MW. The best fit values are the same as those foundin
Yaaqib, Petersen and Penarrubia 2024. The reflex motion is generated for 4 different distances, 20-30, 30-40, 40-50 and 50+ kpc.
All values used here are generated from table 1 in the paper. The reflex motion is then plotted on the sky in Mollweide projection.

In order to use the file, simply provide a file with
"""
vsun_mw = np.array([11.1, 244.24, 7.25])  # km/s
rsun_mw = np.array([-8.3, 0., 0.02])  # kpc


def table2_results():
    """
    main results of paper in Table 2.
    """

    M = np.array([[120., 139., 64., 38.],
                  [19., -59., -47., -37.],
                  [16., 17.,23., 40.],
                  [8., -10., -27., -9.],
                  [-11., -13., -20., -24.],
                  [-9., 9., 20., 17.],
                  [103., 98., 100., 87.],
                  [76., 78., 85., 70.],
                  [61., 83., 96., 74.]]).T
    
    Eu = np.array([[9., 32., 23., 11.],
                   [11., 15., 20., 11.],
                   [2., 5., 7., 7.],
                   [3., 4., 7., 7.],
                   [2., 4., 6., 6.],
                   [3., 5., 8., 7.],
                   [2., 2., 4., 4.],
                   [1., 2., 4., 5.],
                   [1., 3., 4., 5.],]).T
    
    Ed = np.array([[-8., -31., -24., -11.],
                   [-12., -12., -16., -10.],
                   [-2., -5., -7., -7.],
                   [-3.,-4., -7., -6.],
                   [-2., -4., -6., -6.],
                   [-3., -5., -7., -7.],
                   [-2., -2., -4., -4.],
                   [-1., -2., -4., -4.],
                   [-1., -2., -4., -5.]]).T
    
    midpoints = np.array([23.85, 34.31, 44.14, 59.80])
    return M, Eu, np.abs(Ed), midpoints
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


def plot_reflex_model(cube, rgal, vgal, ax=None, quant="vlos", vlosnorm=None,mulnorm=None,mubnorm=None, nobflip=False):

    #lapex, bapex, vtravel, vr, vphi, vtheta,siglos, sigl, sigb

    l,b,dist, vlos, mul, mub = get_v(cube,rgal*300., vgal*240./1.4, solar=False)

    print("the first 5 values of dist, vlos, mul, mub are:")
    print("dist: ", dist[:5])
    print("vlos: ", vlos[:5])
    print("vl: ", mul[:5])
    print("vb: ", mub[:5])

    if ax == None:
        fig, ax = plt.subplots(3, facecolor="white", figsize=(4, 8), subplot_kw={'projection': 'mollweide'})
    else:
        print("using axes defined outside function")
    # fig.patch.set_facecolor('black') #setting plot background to dark colour

    #change phi range to be compatible with the mollwiede projection
    
    if nobflip == True:
        b = np.pi/2. - b
        l = l
    else:
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
            ax.scatter(-l, b, c=cm.coolwarm((mub - np.min(mubnorm))/(np.max(mubnorm)-np.min(mubnorm))), cmap="seismic")
        else:
            ax.scatter(-l, b, c=cm.coolwarm((mub - np.min(mub))/(np.max(mub)-np.min(mub)), cmap="seismic"))

    return ax


def make_apex_data(ax, color="k", marker="s", ebarc="k",label="YPP+24"):
    """
    Figure to reproduce figure 1. in the paper -- without the simulation lines!
    """
    M, Eu, Ed, midpoints = table2_results()

    #convert phi range to 0,2pi for the fitted values
    

    ax[1].set_ylabel(r"$\ell_{\rm apex}$ [deg]")
    
    #apex L
    ax[1].errorbar(midpoints, M[:,0], yerr = [Ed[:,0],Eu[:,0]], color=ebarc, fmt="none", capsize=3)
    ax[1].scatter(midpoints,  M[:,0], c=color, marker=marker,s=20,zorder=100, label=label)
    #V Travel
    ax[0].errorbar(midpoints, M[:,2], yerr = [Ed[:,2],Eu[:,2]], color=ebarc, fmt="none", capsize=3)
    ax[0].scatter(midpoints, M[:,2], c=color, marker=marker,s=20,zorder=100, label=label)
    #b_apex
    ax[2].scatter(midpoints,  M[:,1], c=color, marker=marker,s=20,zorder=100, label=label)
    ax[2].errorbar(midpoints, M[:,1], yerr = [Ed[:,1],Eu[:,1]], color=ebarc,fmt="none", capsize=3)


    ax[0].yaxis.set_minor_locator(MultipleLocator(5))
    ax[0].yaxis.set_major_locator(MultipleLocator(10))
    ax[0].tick_params(axis="x", which="both", labelbottom=False, direction="in")
    ax[0].set_ylabel(r"$v_{travel}$ kms$^{-1}$ ")
    ax[0].set_ylim(1,65)
    
    ax[1].yaxis.set_minor_locator(MultipleLocator(20))
    ax[1].yaxis.set_major_locator(MultipleLocator(40))
    ax[1].tick_params(axis="x", which="both", labelbottom=False, direction="in",top=True)
    
    ax[2].scatter(midpoints, M[:,1], c=color, marker=marker,s=20,zorder=100, label=label)
    ax[2].yaxis.set_minor_locator(MultipleLocator(10))
    ax[2].yaxis.set_major_locator(MultipleLocator(20))
    ax[2].set_ylabel(r"$b_{\rm apex}$ [deg] ")
    ax[2].set_xlabel(r"$r_{\rm galactocentric}$ kpc ")
    ax[2].tick_params(axis="x", which="both", direction="in",top=True)
    
    #setting x axis for all plots
    for i in ax:
        i.tick_params(axis="y", labelsize=8)
        i.set_xlim(15,100)
        i.xaxis.set_minor_locator(MultipleLocator(5))
        i.xaxis.set_major_locator(MultipleLocator(10))
    print(midpoints)
    return ax

def make_bulk_motion_data(ax, color="k", marker="s",ebarc="g",label="YPP+24"):
    M, Eu, Ed, midpoints = table2_results()
    #vr
    ax[0].errorbar(midpoints, M[:,3], yerr = [Ed[:,3],Eu[:,3]], color=ebarc,fmt="none", capsize=3)
    ax[0].scatter(midpoints, M[:,3], c=color, marker=marker,s=20,zorder=100, label=label)
    #vphi
    ax[1].errorbar(midpoints, M[:,4], yerr = [Ed[:,4],Eu[:,4]], color=ebarc,fmt="none", capsize=3)
    ax[1].scatter(midpoints, M[:,4], c=color, marker=marker,s=20,zorder=100, label=label)
    #vtheta
    ax[2].errorbar(midpoints, M[:,5], yerr = [Ed[:,5],Eu[:,5]], color=ebarc,fmt="none", capsize=3)
    ax[2].scatter(midpoints, M[:,5], c=color, marker=marker,s=20,zorder=100, label=label)



    ax[0].yaxis.set_minor_locator(MultipleLocator(5))
    ax[0].yaxis.set_major_locator(MultipleLocator(10))
    ax[0].set_ylabel(r"$v_{r}$ $\rm km ~s^{-1}$")
    ax[0].tick_params(axis="x", which="both", labelbottom=False, direction="in")
    
    ax[1].yaxis.set_minor_locator(MultipleLocator(5))
    ax[1].yaxis.set_major_locator(MultipleLocator(10))
    ax[1].set_ylabel(r"$v_{\phi}$ $\rm km ~s^{-1}$ ")
    ax[1].tick_params(axis="x", which="both", labelbottom=False, direction="in",top=True)
    
    ax[2].yaxis.set_minor_locator(MultipleLocator(5))
    ax[2].yaxis.set_major_locator(MultipleLocator(10))
    ax[2].set_ylabel(r"$v_{\theta}$ $\rm km ~s^{-1}$")
    ax[2].set_xlabel(r"$r_{\rm galactocentric}$ kpc ")
    ax[2].tick_params(axis="x", which="both", direction="in",top=True)

    for i in ax:
        i.tick_params(axis="y", labelsize=8)
        i.set_xlim(15,100)
        i.set_ylim(-40,40)
        i.xaxis.set_minor_locator(MultipleLocator(5))
        i.xaxis.set_major_locator(MultipleLocator(10))
    return ax