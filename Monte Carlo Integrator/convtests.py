import numpy as np
import parameters as p


def R(rate,Es,thetas,phis,thetastar,time):
    if time==0:
        result=np.asarray([0,0])
    else:
        result=np.asarray([0,1])
    for iE in range(len(Es)):
        prefac=p.m_e*np.sqrt(2*p.m_e*Es[iE])
        for it in range(len(thetas)):
            theta=thetas[it]
            if theta>thetastar:
                result[0]+=prefac*np.sin(theta)*np.average(rate[iE,it,0:len(phis)-1])*(np.max(Es)-np.min(Es))*2*np.pi**2/((len(thetas))*len(Es))
            elif theta<np.pi-thetastar and time==0:
                result[1]+=prefac*np.sin(theta)*np.average(rate[iE,it,0:len(phis)-1])*(np.max(Es)-np.min(Es))*2*np.pi**2/((len(thetas))*len(Es))
    return result


def convergedtest(rate, part, Es, thetas, phis, time, tubes):
    if tubes:
        thetastar=np.pi-80.0*np.pi/180.0
    else:
        thetastar=np.pi/2.0
    total=R(rate,Es,thetas,phis,thetastar,time)
    errors=np.zeros(p.nsplit)
    for i in range(p.nsplit):
        partrate=R(part[:,:,:,i], Es, thetas, phis, thetastar, time)
        errors[i]=np.sum(abs(total-partrate)/(total+partrate))
    if np.max(errors)>p.tol:
        return False
    else:
        return True
