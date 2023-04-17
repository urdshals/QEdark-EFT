from integrand import integrand
from integration_bounds import give_integration_bounds_spherical
import numpy as np
import parameters as p
import time
import os
from scipy.ndimage import rotate
import sys

def load_W(high):
    # data_folder, factor_W = 'data', 2
    if high:
        data_folder = 'data/high_ecut'
    else:
        data_folder = 'data/low_ecut'

    kgvec = np.loadtxt(data_folder + '/kG.dat')
    Evec = np.loadtxt(data_folder + '/E.dat') - p.graphene_work_function 
    if high:
        Evec -= p.fermi_energy_high
    else:
        Evec -= p.fermi_energy_low
    if p.dummy_W:
        return kgvec, Evec, 1
    W = np.loadtxt(data_folder + '/W.dat')
    W = np.reshape(W, (100, 100, 50, 100))
    Wsum1 = np.sum(W)
    W[100 - 1, :, :, :] = 0
    W[:, 100 - 1, :, :] = 0
    Wsum2 = np.sum(W)
    # print(Wsum1 / Wsum2)
    W *= Wsum1 / Wsum2

    return kgvec, Evec, W


def mc_integrator(distribution, k_prime,k_prime_abs,nthetas,nphis, W_data_low,W_data_high, m_x, time, c_s, c_l, tubes, size=1000, seed=True):
    # Set the random seed.
    if seed:
        np.random.seed(0.)
    costheta_v_min, costheta_v_max, costheta_q_min, costheta_q_max, phi_q_min, phi_q_max, phi_v_min, phi_v_max, q_min, q_max, delta_E_min, delta_E_max = give_integration_bounds_spherical(k_prime_abs, m_x, time,np.max(W_data_high['kgvec']),np.min([np.min(W_data_low['Evec']),np.min(W_data_high['Evec'])]))

        # Calculate the volume of the hypercube.
    if p.dummy_W:
        integral_volume = (-costheta_v_min + costheta_v_max) * (-costheta_q_min + costheta_q_max) * (-phi_q_min + phi_q_max) * (-phi_v_min + phi_v_max) * (-q_min + q_max)
    else:
        integral_volume = (-costheta_v_min + costheta_v_max) * (-costheta_q_min + costheta_q_max) * (-phi_q_min + phi_q_max) * (-phi_v_min + phi_v_max) * (-q_min + q_max) * (-delta_E_min + delta_E_max)
        # Generate random samples for these variables
    costheta_v = np.random.uniform(low=costheta_v_min, high=costheta_v_max, size=(size, 1))
    phi_v = np.random.uniform(low=phi_v_min, high=phi_v_max, size=(size, 1))
    phi_q = np.random.uniform(low=phi_q_min, high=phi_q_max, size=(size, 1))
    costheta_q = np.random.uniform(low=costheta_q_min, high=costheta_q_max, size=(size, 1))

    
    costheta_vq = np.sqrt(1.-costheta_q*costheta_q)  * np.sqrt(1.-costheta_v*costheta_v) * ( np.cos(phi_q) * np.cos(phi_v) + np.sin(phi_q) * np.sin(phi_v) ) + costheta_q * costheta_v
    costheta_v = costheta_v[costheta_vq>0]
    phi_v = phi_v[costheta_vq>0]
    phi_q = phi_q[costheta_vq>0]
    costheta_q = costheta_q[costheta_vq>0]
    costheta_vq = costheta_vq[costheta_vq>0]

    reduced_size=len(costheta_vq)
    relsize=float(reduced_size)/float(size)
    kplength=np.shape(k_prime)[0]
    nkp = divmod(reduced_size,kplength)[0]
    reduced_size=nkp*kplength
    k_prime_long=np.zeros((reduced_size,3))
    for ikp in range(kplength):
        k_prime_long[ikp*nkp:nkp*(ikp+1),:]=k_prime[ikp,:]


    costheta_v = costheta_v[0:reduced_size].reshape((reduced_size,1))
    phi_v = phi_v[0:reduced_size].reshape((reduced_size,1))
    phi_q = phi_q[0:reduced_size].reshape((reduced_size,1))
    costheta_q = costheta_q[0:reduced_size].reshape((reduced_size,1))
    costheta_vq = costheta_vq[0:reduced_size].reshape((reduced_size,1))

    qnorm = np.random.uniform(low=q_min, high=q_max, size=(reduced_size, 1))
    if p.dummy_W:
        delta_Es = k_prime_abs**2/(2*p.m_e) - p.initial_energy
    else:
        delta_Es = np.random.uniform(low=delta_E_min, high=delta_E_max, size=(reduced_size, 1))


    sintheta_q = np.sqrt(1.-costheta_q*costheta_q)
    sintheta_v = np.sqrt(1.-costheta_v*costheta_v)


    vnorm = (delta_Es + qnorm*qnorm/(2*m_x))/(qnorm*costheta_vq)

    v_x = vnorm * np.cos(phi_v) * sintheta_v
    v_y = vnorm * np.sin(phi_v) * sintheta_v
    v_z = vnorm * costheta_v

    vs = np.c_[v_x, v_y, v_z]

    q_x = qnorm * np.cos(phi_q) * sintheta_q
    q_y = qnorm * np.sin(phi_q) * sintheta_q
    q_z = qnorm * costheta_q
    if tubes:
        angles = np.random.uniform(0, 2*np.pi, size=reduced_size)
    else:
        angles=1
    qs = np.c_[q_x, q_y, q_z]
    evaluated_points_long = distribution(vnorm, vs, qnorm, qs, k_prime_long,k_prime_abs, delta_Es, W_data_low,W_data_high, m_x, time, c_s, c_l, angles, tubes, reduced_size) * qnorm  * vnorm * vnorm / costheta_vq
    rates=np.zeros((nthetas,nphis,1+p.nsplit))
    ikp=0
    for itheta in range(nthetas):
        for iphi in range(nphis):
            evaluated_points=evaluated_points_long[ikp*nkp:nkp*(ikp+1)]
            ikp+=1

            integral_sum = np.sum(evaluated_points)
            parts=np.zeros(p.nsplit+1)

            for isplit in range(p.nsplit):
                parts[isplit+1]=integral_volume*relsize*np.sum(evaluated_points[int(isplit*nkp/float(p.nsplit)) : int((isplit+1)*nkp/float(p.nsplit))])/len(evaluated_points[int(isplit*nkp/float(p.nsplit)) : int((isplit+1)*nkp/float(p.nsplit))])
            parts[0]=integral_volume*integral_sum*relsize/nkp
            rates[itheta,iphi,:]=parts
    return rates


def evaluate_rate(size):
    print("Mc integrator intiated")

    kgvec_out, Evec_out, W_out = load_W()
    W_data = {
        'kgvec': kgvec_out,
        'Evec': Evec_out,
        'W': W_out
    }
    print("W data from QE dark loaded")

    start = time.time()

    phi = 0.
    number_of_thetas = 15

    rate = np.zeros(number_of_thetas)
    thetas = np.linspace(0., np.pi/2., number_of_thetas)
    k_prime_abs = 5000.

    for theta_ind in range(number_of_thetas):
        theta = thetas[theta_ind]
        k_prime = np.array([k_prime_abs*np.cos(phi)*np.sin(theta), k_prime_abs*np.sin(phi)*np.sin(theta), k_prime_abs*np.cos(theta)])
        rate[theta_ind] = mc_integrator(integrand, k_prime, W_data, size=size, seed=False)

    end = time.time()
    print("Elapsed time for ", size, " points ", (end - start)/60.)
    print(thetas, rate)
    return thetas, rate, (end - start)/60.


def plot_and_store_rate(size):
    print('the sample size per theta is ', size)
    thetas, rate, time = evaluate_rate(size)

    np.savetxt('time_' + str(time) + '_size_' + str(size) + '.txt', thetas)
    np.savetxt('plots/rate_' + str(size) + '.txt', rate)
    np.savetxt('plots/thetas_' + str(size) + '.txt', thetas)


if __name__ == '__main__':
    for size in [10**7]:
        plot_and_store_rate(size)

