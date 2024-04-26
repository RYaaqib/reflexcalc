import numpy as np


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

def print_posterior(P):
    params = ["l", "b", "vtravel", "vr", "vphi", "vth", "sigvlos",
              "sigmul", "sigmub"]
    for p in params:
        print('{0:10s} '.format(p), end='')
        print('{0:10.3f}| '.format(np.nanmedian(P[p])), end='')
        # add uncertainties
        print('+{0:<10.3f}/'.format((np.abs(np.nanpercentile(P[p], 86.) - np.nanpercentile(P[p], 50.)))), end='')
        print('-{0:<10.3f} '.format((np.abs(np.nanpercentile(P[p], 14.) - np.nanpercentile(P[p], 50.)))), end='')
        print('')
    return


    
def get_binned_fit_medians(infile, raw=False, cosb=True):

    P = read_posterior(infile, raw, cosb)

    M = np.zeros(9)
    Eu = np.zeros(9)
    Ed = np.zeros(9)

    params = ["l", "b", "vtravel", "vr", "vphi", "vth", "sigvlos",
                "sigmul", "sigmub"]
    
    for j in range(9):
        if j==0:
            M[j] = np.nanmedian(P[params[j]])

            l_samples = np.deg2rad(P[params[j]]) #in degrees
            v_samples = P[params[j+2]]

            func = shift_samples(v_samples, l_samples)
            Eu[j] = func[0]
            Ed[j] = func[1]

        elif j==1:
            b_samples = np.deg2rad(P[params[j]]) #in degrees
            v_samples = P[params[j+1]]

            func = shift_samples(v_samples, b_samples)
            M[j] = np.nanmedian(P[params[j]])
            Eu[j] = func[0]
            Ed[j] = func[1]

        else:

            M[j] = np.nanmedian(P[params[j]])
            Eu[j] = np.abs(np.nanpercentile(P[params[j]],86.)-np.nanpercentile(P[params[j]],50.))
            Ed[j] = np.abs(np.nanpercentile(P[params[j]],14.)-np.nanpercentile(P[params[j]],50.))


    return M, Eu, Ed

def shift_samples(r, theta, all=None):

    th_shift  = np.zeros_like(r)

    #shift angle by the bin with the largest number of samples.
    N, edges = np.histogram(theta[~np.isnan(theta)], bins=np.linspace(-np.pi,np.pi, 37))
    centres = (edges[:-1] + edges[1:])/2.
    idx = np.where(N == np.max(N))[0]
    shift_angle = -centres[idx[0]]
    # shift_angle2 = -np.median(theta[~np.isnan(theta)])

    rot = np.array([[np.cos(shift_angle), -np.sin(shift_angle)],
                     [np.sin(shift_angle), np.cos(shift_angle)]])

    #for each sample compute the x,y,z coordiante.
    for k in range(len(r)):
        rs = np.array([r[k]*np.cos(theta[k]),
          r[k]*np.sin(theta[k])])

        rsp = np.dot(rot, rs)
        #compute the new angles in the rotated frame
        sph = cartesian_to_spherical(rsp[0], rsp[1], 0.)

        th_shift[k] = np.rad2deg(sph[1])

    upper_percentile = np.abs(np.nanpercentile(th_shift,86.)-np.nanpercentile(th_shift,50.))
    lower_percentile = np.abs(np.nanpercentile(th_shift,14.)-np.nanpercentile(th_shift,50.))
    if all == True:
        return th_shift
    else:
        return np.array([upper_percentile, lower_percentile])

def cartesian_to_spherical(x, y, z):
    r = np.sqrt(x ** 2. + y ** 2. + z ** 2.)
    phi = np.arctan2(y, x)
    th = np.arccos(z / r)
    
    return np.array([r, phi, th])

