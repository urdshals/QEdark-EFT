import units as u
import parameters as p
import numpy as np


def give_integration_bounds_spherical(k, m_x, time,kgmax,Emin):
    costheta_v_min, costheta_v_max = -1.,1.
    costheta_q_min, costheta_q_max = -1.,1.
    phi_q_min, phi_q_max = 0., 2.*np.pi
    phi_v_min, phi_v_max = 0., 2.*np.pi
    vmax=p.v_target_abs + p.v_esc
    if p.dummy_W:
        q_min, q_max = m_x*vmax - np.sqrt(m_x**2 * vmax**2 + 2.0 * m_x * p.initial_energy), m_x*vmax + np.sqrt(m_x**2 * vmax**2 + 2.0 * m_x * p.initial_energy)
    else:
        q_min, q_max = m_x*vmax - np.sqrt(m_x**2 * vmax**2 - 2.0 * m_x * p.graphene_work_function), min(k + np.sqrt(3)*kgmax,  m_x*vmax + np.sqrt(m_x**2 * vmax**2 - 2.0 * m_x * p.graphene_work_function))
    delta_E_min, delta_E_max = p.graphene_work_function + k**2/(2. * p.m_e), min(k**2/(2. * p.m_e) - Emin, m_x * vmax**2 / 2)
    return costheta_v_min, costheta_v_max, costheta_q_min, costheta_q_max, phi_q_min, phi_q_max, phi_v_min, phi_v_max, q_min, q_max, delta_E_min, delta_E_max
