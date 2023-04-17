import parameters as p


def cs(lg, c_s, c_l):
     cs = {
         'c1': c_s[1] + c_l[1] * (p.alpha_em * p.m_e) ** 2 / lg['qsq'],
         'c3': c_s[3] + c_l[3] * (p.alpha_em * p.m_e) ** 2 / lg['qsq'],
         'c4': c_s[4] + c_l[4] * (p.alpha_em * p.m_e) ** 2 / lg['qsq'],
         'c5': c_s[5] + c_l[5] * (p.alpha_em * p.m_e) ** 2 / lg['qsq'],
         'c6': c_s[6] + c_l[6] * (p.alpha_em * p.m_e) ** 2 / lg['qsq'],
         'c7': c_s[7] + c_l[7] * (p.alpha_em * p.m_e) ** 2 / lg['qsq'],
         'c8': c_s[8] + c_l[8] * (p.alpha_em * p.m_e) ** 2 / lg['qsq'],
         'c9': c_s[9] + c_l[9] * (p.alpha_em * p.m_e) ** 2 / lg['qsq'],
         'c10': c_s[10] + c_l[10] * (p.alpha_em * p.m_e) ** 2 / lg['qsq'],
         'c11': c_s[11] + c_l[11] * (p.alpha_em * p.m_e) ** 2 / lg['qsq'],
         'c12': c_s[12] + c_l[12] * (p.alpha_em * p.m_e) ** 2 / lg['qsq'],
         'c13': c_s[13] + c_l[13] * (p.alpha_em * p.m_e) ** 2 / lg['qsq'],
         'c14': c_s[14] + c_l[14] * (p.alpha_em * p.m_e) ** 2 / lg['qsq'],
         'c15': c_s[15] + c_l[15] * (p.alpha_em * p.m_e) ** 2 / lg['qsq']
     }
     return cs

