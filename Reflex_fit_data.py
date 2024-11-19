import numpy as np
from pymultinest.solve import solve
import json
import sys
from coord import *
"""
Important:

Spherical unit conventions:
- Phi: azimuthal angle
- Theta: polar angle

Input file should have the following columns, in the following order:
Columns 1-6: x, y, z, vx, vy, vz (galactic cartesian positions and velocities)
Columns 7-12: l, b, distance (heliocentric), line-of-sight velocity, mu_l, mu_b
Columns 13-17: distance_error, vlos_error, mul_error, mub_error, proper motion correlation

The units for the input file are as follows:
- Positions: kpc
- Velocities: km/s
- Distances: kpc
- Line-of-sight velocities: km/s
- Proper motions: mas/yr
- Errors: mas/yr
"""

# Define the solar motion and position in the galactocentric frame
vsun_mw = np.array([11.1, 244.24, 7.25])  # km/s
rsun_mw = np.array([-8.3, 0., 0.02])  # kpc

dfname = sys.argv[1]
prefix = sys.argv[2]


# 1. Define the galactic cartesian coordinates from input file of shape (cols,rows)
d = np.loadtxt(dfname).T
#uncomment this if file has shape rows,cols
# d = np.loadtxt(dfname)
x, y, z = d[:, 0], d[:, 1], d[:, 2]
vx, vy, vz = d[:, 3], d[:, 4], d[:, 5]

rgal = np.zeros((len(x), 3))
rgal[:, 0] = x
rgal[:, 1] = y
rgal[:, 2] = z

vgal = np.zeros((len(x), 3))
vgal[:, 0] = vx
vgal[:, 1] = vy
vgal[:, 2] = vz

# 2. compute spherical galactic coordinates
rgalsph = np.zeros((len(x), 3))
tmp = cartesian_to_spherical(rgal[:, 0], rgal[:, 1], rgal[:, 2],
                             vgal[:, 0], vgal[:, 1], vgal[:, 2])
rgalsph[:, :] = tmp[0].T

def get_v(cube, rgal, vgal):
    """
    Function to compute and add the reflex motion to the velocities of the stars from the hypercube parameters

    :param cube: Multinest hypercube (1xNparam)
    :param rgal: galactic cartesian coordinates
    :param vgal: galactic cartesian velocities

    :return: vlos, mul, mub, the line-of-sight velocity, proper motion in l and b with reflex motion and bulk motion added
    """
    # 3. rotate cartesian coordinates such that the z axis points at lapex, bapex
    # through the euler angle rotation x-y-z
    bapex = np.arccos(cube[1])
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
    # 4i. add vtravel in the rotated frame
    vpsph[:, :] = add_vtravel(cube[2], rpsph[:, 2]).T 

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


    rgal1[:, :] = np.dot(np.linalg.inv(Rrot), rp1[:, :].T).T
    vgal1[:, :] = np.dot(np.linalg.inv(Rrot), vp1[:, :].T).T

    # 7. compute the vector v and translate to be at the sun
    r = np.zeros_like(rp)
    v = np.zeros_like(vp)

    # 7i. convert bulk velocities to cartesian velocities
    tmp = spherical_to_cartesian(rgalsph[:, 0], rgalsph[:, 1],
                                 rgalsph[:, 2], cube[3], cube[4], cube[5])
    r[:, :] = rgal1[:, :] - rsun_mw
    # 7ii. sum the velocities from the bulk motion(tmp[1]) and reflex motion(vgal1) and correct for solar motion
    v[:, :] = vgal1[:, :] + tmp[1].T - vsun_mw

    # 8. find the galactic proper motions using the unit vectors of ul, ub, ulos
    # spherical unit vectors in the galactic frame..
    rsunsph = np.zeros_like(r)
    vsunsph = np.zeros_like(v)
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
    # 9. compute the observables
    vlos = np.einsum('ij,ij->i', v, elos)
    mul = np.einsum('ij,ij->i', v, ephi) / fac
    mub = -np.einsum('ij,ij->i', v, eth) / fac
    return vlos, mul, mub


# 9. Define the likelihood function for vlos
def like_vlos(cube, vlos_data, vlos_param, evlos):
    """
    Function to compute the line-of-sight velocity likelihood

    :param cube: Multinest hypercube (1xNparam)
    :param vlos_data: observational line-of-sight velocity
    :param vlos_param: model line-of-sight velocity
    :param evlos: observational line-of-sight velocity error

    :return: plos, the likelihood of the line-of-sight velocity
    """
    sigvlos = 1. / cube[6]
    evlos2 = evlos ** 2. + sigvlos
    exponent =  -0.5 * (1. /evlos2)  * ((vlos_data - vlos_param) ** 2.)
    ln = - 0.5 * np.log(2 * np.pi * evlos2)
    plos = ln + exponent

    return plos

def like_pms(cube, pml_data, pmb_data, dist, corr, epml, epmb, edist, pml_param, pmb_param):
    """
    Function to compute the proper motion likelihood

    :param cube: Multinest hypercube (1xNparam)
    :param pml_data: observational proper motion in l
    :param pmb_data: observational proper motion in b
    :param dist: heliocentric distance
    :param corr: correlation between proper motions
    :param epml: proper motion in l error
    :param epmb: proper motion in b error
    :param edist: distance error
    :param pml_param: model proper motion in l
    :param pmb_param: model proper motion in b

    :return: ppm, the likelihood of the proper motion
    """

    sigpml = 1. / cube[7]
    sigpmb = 1. / cube[8]
    fac = 4.74057 * dist
    elp2 = epml ** 2. + edist ** 2. * (np.abs(pml_data) ** 2. / dist ** 2.) + (sigpml / fac ** 2.)
    ebp2 = epmb ** 2. + edist ** 2. * (np.abs(pmb_data) ** 2. / dist ** 2.) + (sigpmb / fac ** 2.)

    cov = np.array([[elp2, epml * epmb * corr],
                    [epml * epmb * corr, ebp2]])

    kl = pml_data - pml_param
    kb = pmb_data - pmb_param

    #matrix version of above
    X = np.array([kl,kb])
    det = np.linalg.det(cov)
    exp = -0.5 * np.dot(X.T,np.dot(np.linalg.inv(cov), X))
    ln = -0.5 * np.log(((2 * np.pi) ** 2.) * np.linalg.det(cov))
    ppm = ln + exp

    return ppm

# 10. compute the likelihood function for the data and model..
def LogLikelihood(cube):
    """
    Function to compute the log likelihood of the model given the data

    :param cube: Multinest hypercube (1xNparam)

    :return: lnptot, the total log likelihood
    """

    lnptot = 0.

    vlos, mul, mub = get_v(cube, rgal, vgal)
    for i in range(len(x)):
        a = like_vlos(cube, d[i, 9], vlos[i], d[i, 13]) + \
            like_pms(cube, d[i, 10], d[i, 11], d[i, 8], d[i,16], d[i, 14], d[i, 15], d[i, 12], mul[i], mub[i])
        lnptot += a
        if np.isinf(lnptot):
            lnptot = 1.e-160
    return lnptot


# define prior
def Prior(cube):
    """
    Function to define the prior for the hypercube parameters

    :param cube: Multinest hypercube (1xNparam)

    :note: the prior parameters are defined as follows:
    :note: [0] l_apex (rad)
    :note: [1] cos(b_apex)
    :note: [2] v_travel (kms^-1)
    :note: [3] vr (kms^-1)
    :note: [4] vphi (kms^-1)
    :note: [5] vth (kms^-1)
    :note: [6] sigvlos (kms^-1)
    :note: [7] sigmul (mas/yr)
    :note: [8] sigmub (mas/yr)

    :note: the prior ranges are defined as follows:
    :note: [0] -pi < l_apex < pi                           flat prior in l_apex
    :note: [1] -1 < cos(b_apex) < 1                        flat prior in cos(b_apex)
    :note: [2] 0 < v_travel < 150                          flat prior in v_travel
    :note: [3,4,5] -250 < vr,vphi,vth < 250                flat priors in bulk motion
    :note: [6,7,8] 1.e-6 < sigvlos,sigmul,sigmub < 1.e-2   Jeffreys priors
    """
    # a + (b-a) * cube rane a < x < b
    # Order:  l,b,vtravel, vr,vphi,vth,svlos,sl,sb
    cube[0] = -np.pi + cube[0] *  2. * np.pi  # l_apex (rad)
    cube[1] = -1. + cube[1]*2.  # b_apex (rad)
    cube[2] = cube[2] * (150.)  # v_travel (kms^-1)
    cube[3] = - 250. + cube[3] * (500.)  # vr
    cube[4] = - 250. + cube[4] * (500.)  # vphi
    cube[5] = - 250. + cube[5] * (500.)  # vth
    maxsig2 = 1.e-1
    minsig2 = 1.e-4
    cube[6] = (minsig2) ** 2. + cube[6] * ((maxsig2) ** 2. - ((minsig2) ** 2.))  # sigvlos
    cube[7] = (minsig2) ** 2. + cube[7] * ((maxsig2) ** 2. - ((minsig2) ** 2.))  # sigpml
    cube[8] = (minsig2) ** 2. + cube[8] * ((maxsig2) ** 2. - ((minsig2) ** 2.))  # sipmb

    return cube


# Running the sampler
parameters = ["l", "b", "vtravel", "vr", "vphi", "vth", "sigvlos",
              "sigmul", "sigmub"]

n_params = len(parameters)
result = solve(LogLikelihood=LogLikelihood, Prior=Prior,
               n_dims=n_params, outputfiles_basename=prefix, verbose=True,
               resume=False, n_live_points=1000, wrapped_params=None,
               n_iter_before_update=100)

print()
print('evidence: %(logZ).1f +- %(logZerr).1f' % result)
print()
print('parameter values:')
for name, col in zip(parameters, result['samples'].transpose()):
    print('%15s : %.3f +- %.3f' % (name, col.mean(), col.std()))

# make marginal plots by running:
# $ python multinest_marginals.py chains/3-
# For that, we need to store the parameter names:
with open('%sparams.json' % prefix, 'w') as f:
    json.dump(parameters, f, indent=2)
