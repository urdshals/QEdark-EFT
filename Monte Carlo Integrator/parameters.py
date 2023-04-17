import numpy as np
import units as u
from math import erf


# Test ############################################################################################################

dummy_W=False
initial_energy = -9 * u.eV #Initial energy of the electrons

# MC parameters ####################################################################################################

tol=0.005
nsplit=10
sizestep=3
sizemax=58
sizemin=0.1
mode='improve'
# Velocities #######################################################################################################
fermi_energy_high = -2.7
fermi_energy_low = -4.14
v0 = 238. * u.km/u.sec  # in km/s
v_target_abs = 250.5 * u.km/u.sec  # in km/s
alpha_target = 42. * np.pi / 180.  # in radian
v_esc = 544. * u.km/u.sec  # in km/s
m_Sheet = 1.*u.kg
m_e = 510.99895 * u.keV  # mass of electron in eV
m_cell = 2. * 12. * u.AMU
j_x = 1./2.
interpolation_type = 'linear'
alpha_em = 1./137.
graphene_work_function = 4.3

# Calculated variables ############################################################################################


def n_x(m_x):
    return 0.4*(u.GeV/u.cm**3)/m_x


def mu(m_x):
    return m_x*m_e/(m_x + m_e)


def N_esc():
    return erf(v_esc/v0) - 2.*(v_esc/v0)*np.exp(-v_esc**2/v0**2)/np.sqrt(np.pi)


def N_v():
    return 1./(N_esc()*np.pi**(3./2.)*v0**3)  # calculate from the velocity prefactor.


def beta(time):
    return 2. * np.pi * time/24.


def v_target_x(time):
    return np.sin(alpha_target) * np.sin(beta(time))


def v_target_y(time):
    return np.sin(alpha_target) * np.cos(alpha_target) * (np.cos(beta(time)) - 1.)


def v_target_z(time):
    return np.cos(alpha_target)**2 + np.sin(alpha_target)**2 * np.cos(beta(time))


def v_target(time):
    return v_target_abs * np.array([v_target_x(time), v_target_y(time), v_target_z(time)])


def integration_factor(m_x):
    return N_v() * n_x(m_x) * m_Sheet / (32. * np.pi**2 * m_x**2 * m_e**2 * m_cell)





