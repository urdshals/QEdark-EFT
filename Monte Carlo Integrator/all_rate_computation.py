import numpy as np
from monte_carlo import mc_integrator, load_W
from integrand import integrand
import parameters as p
import units as u
import convtests as c
import sys
import os
import time as t
from joblib import Parallel, delayed
import mass_list_merger as ma


t0=t.time()

def get_grid(number_of_thetas, number_of_phis, number_of_Es,mass):
    thetas = np.linspace(0., np.pi, number_of_thetas)
    phis = np.linspace(-np.pi, np.pi, number_of_phis)
    Es = np.linspace(0, 0.5*mass*u.MeV*(p.v_esc+p.v_target_abs)**2 - p.graphene_work_function, number_of_Es+2)

    np.savetxt('data/Esfolder/Es'+str(mass)+'.txt', Es[1:number_of_Es+1])
    np.savetxt('data/thetas.txt', thetas)
    np.savetxt('data/phis.txt', phis)
    return thetas, phis, Es[1:number_of_Es+1]


def all_rate_computation(Es, thetas, phis, W_data_low,W_data_high, m_x, time, c_s, c_l, tubes,sizefactor):
    rate = np.zeros((len(Es), len(thetas), len(phis), 1 + p.nsplit))
    for E_ind in range(len(Es)):
        k_prime_abs = np.sqrt(2. * p.m_e * Es[E_ind])
        k_prime=np.zeros((len(thetas)*len(phis),3))
        k_prime_index=0
        for theta_ind in range(len(thetas)):
            theta = thetas[theta_ind]
            for phi_ind in range(len(phis)):
                phi = phis[phi_ind]
                k_prime[k_prime_index,:] = np.array([k_prime_abs*np.cos(phi)*np.sin(theta), k_prime_abs*np.sin(phi)*np.sin(theta), k_prime_abs*np.cos(theta)])
                k_prime_index+=1
        
        nrounds=divmod(sizefactor,p.sizemax)[0]+1
        size = int(sizefactor*10**6/nrounds)
        for i in range(int(nrounds)):
            rate[E_ind, :, :,:] += mc_integrator(integrand, k_prime,k_prime_abs,len(thetas),len(phis), W_data_low,W_data_high, m_x, time, c_s, c_l, tubes, size=size, seed=False)
    return rate/nrounds


def run_calculation():
    kgvec_out, Evec_out, W_out = load_W(False)
    W_data_low = {
        'kgvec': kgvec_out,
        'Evec': Evec_out,
        'W': W_out
    }
    kgvec_out, Evec_out, W_out = load_W(True)
    W_data_high = {
        'kgvec': kgvec_out,
        'Evec': Evec_out,
        'W': W_out
    }
    # for each of the calculations
    sizefactor=1
    cont=False
    tubes
    runname
    if p.dummy_W: 
        runname+='_dummy_W'
    print('Tubes: ', tubes)
    if os.path.isfile('temp_rates/'+runname+'_temp.npy') and p.mode!='from_scratch':
        all_rates = np.load('temp_rates/'+runname+'_temp.npy',allow_pickle=True)
        all_rates = all_rates.flat[0]
    else:
        all_rates = {}
    for coupling_nr in [1, 16, 17]:
        #mass_list_cont = np.logspace(np.log10(2*p.graphene_work_function/(p.v_esc+p.v_target_abs)**2/u.MeV), 2, 30)
        #mass_list_disc = [2., 5., 10., 20., 50., 100.]
        #mass_list = ma.listmerger(mass_list_cont,mass_list_disc)
        mass_list = [2.]
        for mass in mass_list:
            thetas, phis, Es = get_grid(50, 50, 64, mass)
            for interaction in ['l', 's']:
                for time in range(24):
                    # set the parameters
                    m_x = mass * u.MeV
                    c_s = np.zeros(16)
                    c_l = np.zeros(16)
                    if interaction == 's' and coupling_nr<16:
                        c_s[coupling_nr] = 1.
                    elif interaction == 'l' and coupling_nr<16:
                        c_l[coupling_nr] = 1.
                    coupling_run_name = 'c' + str(coupling_nr) + '_' + interaction

                    if coupling_nr == 16:
                        c_s = np.zeros(len(c_s))
                        c_l = np.zeros(len(c_l))
                        if interaction == 's':
                            return True
                            coupling_run_name = 'hochberg'
                            c_s[1] = np.sqrt(10**-37 * u.cm**2 * 16. * np.pi * m_x**2 * p.m_e**2 / p.mu(m_x)**2)
                        elif interaction == 'l':
                            coupling_run_name = 'anapole'
                            glam = 1 * u.GeV ** -2
                            c_s[8] = 8. * u.ElementaryCharge * p.m_e * m_x * glam
                            c_s[9] = -8. * u.ElementaryCharge * p.m_e * m_x * glam
                    elif coupling_nr == 17:
                        c_s = np.zeros(len(c_s))
                        c_l = np.zeros(len(c_l))
                        if interaction == 's':
                            coupling_run_name = 'magnetic_dipole'
                            glam = 1 *  u.GeV ** -1
                            c_s[1] = 4. * u.ElementaryCharge * p.m_e * glam
                            c_s[4] = 16. * u.ElementaryCharge * m_x * glam
                            c_l[5] = 16. * u.ElementaryCharge * p.m_e**2 * m_x * glam/(p.alpha_em * p.m_e)**2
                            c_l[6] = -16. * u.ElementaryCharge * p.m_e**2 * m_x * glam/(p.alpha_em * p.m_e)**2
                        elif interaction == 'l':
                            coupling_run_name = 'electric_dipole'
                            glam = 1 * u.GeV ** -1
                            c_l[11] = 16. * u.ElementaryCharge * p.m_e**2 * m_x * glam/(p.alpha_em * p.m_e)**2
                    # run the calculation
                    paralell = True
                    if not coupling_run_name in all_rates:
                        all_rates[coupling_run_name] = {}
                    if not str(mass) in all_rates[coupling_run_name]:
                        all_rates[coupling_run_name][str(mass)] = {}
                    if not str(time) in  all_rates[coupling_run_name][str(mass)] or p.mode=='improve':

                        if mass==mass_list[0] and cont==True:
                            rate = np.zeros((len(Es), len(thetas), len(phis)))
                            rate = np.array(rate).squeeze()
                        else:
                            if paralell:
                                if p.mode=='improve':
                                    rate,sizefactor=all_rates[coupling_run_name][str(mass)][str(time)]
                                    converged=False
                                else:
                                    if sizefactor<p.sizemin:
                                        sizefactor=p.sizemin
                                    new_rate=Parallel(n_jobs=len(Es))(delayed(all_rate_computation)([E], thetas, phis, W_data_low,W_data_high, m_x, float(time), c_s, c_l, tubes,sizefactor) for E in Es)
                                    new_rate = np.array(new_rate).squeeze()
                                    rate=new_rate[:,:,:,0]
                                    if np.sum(rate)>0:
                                        converged = c.convergedtest(rate, new_rate[:,:,:,1:p.nsplit+1], Es, thetas, phis,time,tubes)
                                    else:
                                        converged = False
                                while not converged:
                                    sizefactor *= p.sizestep
                                    new_rate = Parallel(n_jobs=len(Es))(delayed(all_rate_computation)([E], thetas, phis, W_data_low,W_data_high, m_x, float(time), c_s, c_l, tubes,sizefactor) for E in Es)
                                    new_rate = np.array(new_rate).squeeze()
                                    rate += new_rate[:, :, :, 0] * p.sizestep
                                    rate /= float(p.sizestep + 1)
                                    if np.sum(rate)>0:
                                        converged = c.convergedtest(rate, new_rate[:,:,:,1:p.nsplit+1], Es, thetas, phis,time,tubes)
                                    if sizefactor>3000:
                                        converged=True
                                        print('Failed to reach convergence')
                                sizefactor /= float(p.sizestep)
                            else:
                                rate = all_rate_computation(Es, thetas, phis, W_data, m_x, float(time), c_s, c_l)

                        all_rates[coupling_run_name][str(mass)][str(time)] = rate, sizefactor*(p.sizestep + 1)
                    # tell me where I am
                        print(mass, coupling_nr, interaction, time, sizefactor, np.sum(rate), t.time()-t0)
                        sys.stdout.flush()
                        np.save('temp_rates/'+runname+'_temp.npy',all_rates)
        np.save('final_rates/'+runname+'_final.npy',all_rates)
    return True


if __name__ == '__main__':
    run_calculation()

