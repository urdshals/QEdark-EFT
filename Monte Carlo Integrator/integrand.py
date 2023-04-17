import numpy as np
import parameters as p
from calculate_M import M_prime_squared
from scipy import ndimage
import units as u
import sys

def find_whole(array, values):
    ixs = np.ceil(find_fractional(array, values))
    return ixs


def find_fractional(array, values):
    step = (array[-1] - array[0])/(len(array) - 1)
    ixs = (values - array[0])/step
    return ixs


def heaviside(v,time):
    argument = p.v_esc - np.linalg.norm(v + p.v_target(time), axis=1, keepdims=True)
    return np.heaviside(argument, 0)


def exponential(v,time):
    argument = -np.linalg.norm(v + p.v_target(time), axis=1, keepdims=True)**2/p.v0**2
    return np.exp(argument)



def W_interpolate(momentum, energy, W_data, interpolation_type='linear'):
    if p.interpolation_type == 'none':
        index_x = find_whole(W_data['kgvec'], momentum[:,[0]])
        index_y = find_whole(W_data['kgvec'], momentum[:,[1]])
        index_z = find_fractional(np.flip(W_data['kgvec'][W_data['kgvec'] >= 0]), abs(momentum[:, [2]]))
        index_E = find_whole(W_data['Evec'], energy)

        indexes = np.c_[index_x, index_y, index_z, index_E]

        W_interpolated = ndimage.map_coordinates(W_data['W'], indexes.T, order=1,mode='nearest')
    elif p.interpolation_type == 'linear':
        index_x = find_fractional(W_data['kgvec'], momentum[:,[0]])
        index_y = find_fractional(W_data['kgvec'], momentum[:,[1]])
        index_z = find_fractional(np.flip(W_data['kgvec'][W_data['kgvec'] >= 0]), abs(momentum[:, [2]]))
        index_E = find_fractional(W_data['Evec'], energy)

        indexes = np.c_[index_x, index_y, index_z, index_E]
        W_interpolated = ndimage.map_coordinates(W_data['W'], indexes.T, order=1,mode='nearest')
    else:
        raise Exception("You specified an incorrect interpolation type!")

    return W_interpolated.reshape(len(W_interpolated), 1)

def W_dummy(momentum):
    argument=-(momentum[:,[0]]**2 + momentum[:,[1]]**2 + 0.5*momentum[:,[2]]**2)/(p.m_e*p.alpha_em)**2
    return np.exp(argument)/p.m_e**3

def integrand(vnorm, v, qnorm,q, k_prime,k_prime_abs, deltaE, W_data_low,W_data_high, m_x, time, c_s, c_l, angles , tubes,size):
    list_of_q_dots = qnorm*qnorm
    list_of_v_dots = vnorm*vnorm
    kqnorm=np.linalg.norm(k_prime - q, axis=1, keepdims=True)
    kqnorm=np.reshape(kqnorm,len(kqnorm))
    lowkmax=(W_data_low['kgvec'])[len(W_data_low['kgvec'])-2]
    linalg_precomputes = {
        'qsq': list_of_q_dots,
        'vsq': list_of_v_dots
    }
    m_term = M_prime_squared(v, q, k_prime, deltaE, linalg_precomputes, m_x, c_s, c_l)
    if tubes:
        starting_momentum = k_prime - q
        rotated_momentum = np.zeros_like(starting_momentum)
        rotated_momentum[:, 0] = starting_momentum[:, 0] * np.cos(angles[:]) + starting_momentum[:, 2] * np.sin(angles[:])
        rotated_momentum[:, 1] = starting_momentum[:, 0] * np.sin(angles[:]) + starting_momentum[:, 2] * np.cos(angles[:])
        rotated_momentum[:, 2] = (starting_momentum)[:, 1]
        if p.dummy_W:
            W_term = W_dummy(rotated_momentum)
        else:
            W_term = np.zeros_like(deltaE)
            if len(W_term[kqnorm<lowkmax])>0:
                W_term[kqnorm<lowkmax] = W_interpolate(rotated_momentum[kqnorm<lowkmax,:], (k_prime_abs**2/(2.*p.m_e) - deltaE)[kqnorm<lowkmax], W_data_low, interpolation_type=p.interpolation_type)
            if len(W_term[kqnorm>=lowkmax])>0:
                W_term[kqnorm>=lowkmax] = W_interpolate(rotated_momentum[kqnorm>=lowkmax,:], (k_prime_abs**2/(2.*p.m_e) - deltaE)[kqnorm>=lowkmax], W_data_high, interpolation_type=p.interpolation_type)
    else:
        if p.dummy_W:
            W_term = W_dummy(k_prime - q)
        else:
            W_term = np.zeros_like(deltaE)
            if len(W_term[kqnorm<lowkmax])>0:
                W_term[kqnorm<lowkmax] = W_interpolate((k_prime - q)[kqnorm<lowkmax,:], (k_prime_abs**2/(2.*p.m_e) - deltaE)[kqnorm<lowkmax], W_data_low, interpolation_type=p.interpolation_type)
            if len(W_term[kqnorm>=lowkmax])>0:
                W_term[kqnorm>=lowkmax] = W_interpolate((k_prime - q)[kqnorm>=lowkmax,:], (k_prime_abs**2/(2.*p.m_e) - deltaE)[kqnorm>=lowkmax], W_data_high, interpolation_type=p.interpolation_type)
    return u.yr * p.integration_factor(m_x)*exponential(v,time) * heaviside(v,time) * m_term * W_term


if __name__ == '__main__':
    from monte_carlo import load_W
    kgvec_out, Evec_out, W_out = load_W()
    W_data = {
        'kgvec': kgvec_out,
        'Evec': Evec_out,
        'W': W_out
    }
    print(integrand([0.], [0.], np.array([[0., 0., 10000.]]), np.array([0., 0., 1152.]), 10., W_data))
    # print(integrand([0.], [0.], np.array([[0., 0., 100000.]]), np.array([0., 0., 13934.]), 200., W_data))




