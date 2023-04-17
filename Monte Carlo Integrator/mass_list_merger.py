import numpy as np
def min_index(vec):
    value=vec[0]
    index=0
    for i in range(len(vec)):
        if vec[i]<value:
            value=vec[i]
            index=i
    return index

def listmerger(mass_list_cont,mass_list_disc):
    for i in range(len(mass_list_disc)):
        index=min_index(abs(np.log10(mass_list_cont)-np.log10(mass_list_disc[i])))
        mass_list_cont[index]=mass_list_disc[i]
    return mass_list_cont



