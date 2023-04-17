import numpy as np
import parameters as p
from CsParams import cs


def M_squared(v, q, k_prime, deltaE, VHelpers, Cs, lg):
    q_cross_vperp_sq = VHelpers['q_cross_vel_sq']
    q_me_sq = VHelpers['q_me_sq']
    v_perp_sq = VHelpers['v_perp_sq']
    q_dot_vperp_sq = VHelpers['v_perp_q']**2/p.m_e**2

    first_four = Cs['c1']**2 + Cs['c3']**2/4. * q_cross_vperp_sq + Cs['c7']**2/4. * v_perp_sq + Cs['c10']**2/4.*q_me_sq
    factor = p.j_x*(p.j_x + 1.)/12.
    first_four_inside = 3.*Cs['c4']**2 + (4.*Cs['c5']**2 - 2.*Cs['c12']*Cs['c15']) * q_cross_vperp_sq + Cs['c6']**2 * q_me_sq**2 + (4.*Cs['c8']**2 + 2.*Cs['c12']**2)*v_perp_sq
    second_four_inside = (2.*Cs['c9']**2 + 4.*Cs['c11']**2 + 2.*Cs['c4']*Cs['c6']) * q_me_sq + (Cs['c13']**2 + Cs['c14']**2) * q_me_sq * v_perp_sq + Cs['c15']**2 * q_me_sq * q_cross_vperp_sq + 2.*Cs['c13']*Cs['c14']*q_dot_vperp_sq

    return first_four + factor * (first_four_inside + second_four_inside)


def me_real_part(v, q, k_prime, deltaE, VHelpers, Cs, lg):
    v_perp = VHelpers['v_perp']
    q_me_sq = VHelpers['q_me_sq']
    q_dot_vperp = VHelpers['v_perp_q']/p.m_e
    qk_prime = VHelpers['qk_prime']

    first_bracket = Cs['c3']**2/2.*(q_dot_vperp * q/p.m_e - q_me_sq*v_perp) - Cs['c7']**2/2.*v_perp
    first_bracket_inside = (4.*Cs['c5']**2 + Cs['c15']**2*q_me_sq) * (q_dot_vperp*q/p.m_e - q_me_sq*v_perp) - (4.*Cs['c8']**2 + 2.*Cs['c12']**2 + (Cs['c13']**2 + Cs['c14']**2)*q_me_sq)*v_perp
    second_bracket_inside = 2.*Cs['c12']*Cs['c15']*(q_dot_vperp*q/p.m_e - q_me_sq*v_perp)
    third_bracket_inside = 2.*Cs['c13']*Cs['c14']*q_dot_vperp*q/p.m_e

    inside_total_vec = first_bracket_inside - second_bracket_inside - third_bracket_inside
    return np.sum((first_bracket + p.j_x*(p.j_x + 1.)/6. * inside_total_vec) * qk_prime, axis=1, keepdims=True)


def me_squared_part(v, q, k_prime, deltaE, VHelpers, Cs, lg):
    q_me_sq = VHelpers['q_me_sq']
    qk_prime = VHelpers['qk_prime']
    qk_prime_sq = np.linalg.norm(qk_prime, axis=1, keepdims=True)**2
    q_qk_prime = np.sum(q / p.m_e * qk_prime, axis=1, keepdims=True)

    first_two_terms = (Cs['c3']**2/4. * q_me_sq + Cs['c7']**2/4.)*qk_prime_sq - Cs['c3']**2/4. * q_qk_prime**2
    first_inside = qk_prime_sq*((4.*Cs['c5']**2 + Cs['c13']**2 + Cs['c14']**2 - 2.*Cs['c12']*Cs['c15'])*q_me_sq + 4.*Cs['c8']**2 + 2.*Cs['c12']**2 + Cs['c15']**2*q_me_sq**2)
    second_inside = q_qk_prime**2 * (-4.*Cs['c5']**2 - Cs['c15']**2*q_me_sq + 2.*Cs['c12']*Cs['c15'] + 2.*Cs['c13']*Cs['c14'])
    return first_two_terms + p.j_x*(p.j_x + 1.)/12. * (first_inside + second_inside)


def M_prime_squared(v, q, k_prime, deltaE, lg, m_x, c_s, c_l):
    Cs = cs(lg, c_s, c_l)

    VHelpers = {
        'v_perp_sq': lg['vsq'] + lg['qsq'] / (4. * p.mu(m_x) ** 2) * (m_x - p.m_e) / (m_x + p.m_e) - deltaE / p.mu(m_x),
        'v_perp': v - q / (2. * p.mu(m_x)),
        'v_perp_q': deltaE - lg['qsq'] / (2. * p.m_e),
        'q_me_sq': lg['qsq'] / (p.m_e ** 2),
        'qk_prime': (q - k_prime) / p.m_e
    }

    VHelpers['q_cross_vel_sq'] = lg['qsq'] / (p.m_e ** 2) * VHelpers['v_perp_sq'] - (np.sum(q * VHelpers['v_perp'], axis=1, keepdims=True) / p.m_e) ** 2

    first = M_squared(v, q, k_prime, deltaE, VHelpers, Cs, lg)
    second = me_real_part(v, q, k_prime, deltaE, VHelpers, Cs, lg)
    third = me_squared_part(v, q, k_prime, deltaE, VHelpers, Cs, lg)
    return first + second + third

