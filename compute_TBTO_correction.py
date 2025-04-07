import numpy as np
import matplotlib.pyplot as plt
import scipy.io

def compute_exponential_fit(v1,v2,tau,t2,t):
    
    #let the exponential fit function be as follows:
    # v=alpha(1-e(-(t+beta)/tau))
    # with v the parameter value, t the time axis, tau the t value when the tangent at the origin crosses v(t=inf)
    # alpha and beta are the function coefficient
    # let v1(t1=0), v2(t2) and tau:
    vr=v1/v2
    beta=-tau*np.log((vr-1)/(vr*np.exp(-t2/tau)-1))
    alpha=v1/(1-np.exp(-beta/tau))
    
    v_fit=alpha*(1-np.exp(-(t+beta)/tau))
    
    return v_fit


wmo='3901999'
psal_on_theta_dir="C:/Users/ddobler/Documents/08_DD_scripts/06_OWC_matlab/05_OWC_OUTPUTS/3901999_minus_0.067/"
psal_on_theta_ds = scipy.io.loadmat(psal_on_theta_dir + 'cal_psal_on_theta_'+wmo+'.mat')
vint=psal_on_theta_ds["Sint"][:]
ntheta,ncycle=vint.shape
print("ntheta,ncycle=",ntheta,ncycle)

tau=3
t2=7
    
correction=np.nan*np.ones((ntheta,t2))

for itheta in range(ntheta):
    print("itheta=",itheta)
    vint_th=vint[itheta,:]

    if ~np.isnan(np.min(vint_th[:t2])):
        # theta level 2 (from savings of Sint in OWC plot function)
        # # vint=np.array([34.6579984741211, 34.6606178917931, 34.6620000915527, 34.6627138550698, 34.6629940890433, 34.6639989929199, 34.6639989929199, 34.6639989929199, 34.6639989929199,34.6639989929199,34.6639989929199,34.6649984436035,34.6639989929199,34.6647237894753,34.6649984436035,34.6642507980733,34.6649984436035,34.6649984436035,34.6649984436035,34.6651971028820,34.6639989929199,34.6639989929199,34.6639989929199,34.6649984436035])
        t=np.arange(len(vint_th))


        v_fit=compute_exponential_fit(vint_th[0],vint_th[t2],tau,t2,t)
        rmse=np.sqrt((1/t2)*np.sum(v_fit[:t2]-vint_th[:t2])**2)
        print("rmse=",rmse)

        correction[itheta,:] = np.nanmean(vint_th[t2:])-v_fit[:t2]
        print(correction[itheta,:])

        # plt.figure()
        # plt.plot(t,v_fit,c='b')
        # plt.scatter(t,vint_th,c='k')
        # plt.grid()
        # plt.title("itheta=" + str(itheta) + "\n" + "fit_rmse={:.3f}".format(rmse))
        # plt.show()
        
print("correction=",correction)

print("np.nanmean(correction,axis=0)")
print(np.nanmean(correction,axis=0))
print("np.nanstd(correction,axis=0)")
print(np.nanstd(correction,axis=0))

