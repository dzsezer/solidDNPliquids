import numpy as np

def delta2(closest_distance_nm, radical_conc_molar):
    #dipolar constant
    D_dipolar = 2*np.pi * (79.066) # nm^3/ms
    Ddip_per_us = 1e-3 * D_dipolar # nm^3/us

    #converstion from moles/lt to 1/nm^3
    number_per_nm3 = 6.02214e23*radical_conc_molar/1e24
    N_3b3 = number_per_nm3/(3 * closest_distance_nm**3) # 1/nm^6

    delta2 = 6*np.pi/5 * Ddip_per_us**2 * N_3b3 # 1/us^2
    return delta2


def omega1_from_B1(B1):
    #takes B1 in Gauss
    #returns omega1 in M rad/s
    return 2*np.pi* 2.8 * B1


def bloch_ss(R1,R2,omega1,Omega):
    #Steady state of the Bloch equations
    # Equations (11) and (15)
    fy = 1/(R2 + Omega**2/R2)
    fx = Omega/R2 * fy
    fz  = 1/(R1 + omega1**2 * fy)
    return fx, fy, fz


def new_bloch_ss(R1,R2,omega1,omegaI,Omega):
    #Steady state of new Bloch equations
    # Equations (38), (40) and (53)
    Fy = 1/(R2+1j*omegaI + Omega**2/(R2+1j*omegaI))
    Fx = Omega/(R2+1j*omegaI) * Fy
    Fz = 1/( R1+1j*omegaI + omega1**2 *Fy )
    Fy_prime = ( 1 + R2/(R2+1j*omegaI) ) * Fy    
    return Fx, Fy, Fz, Fy_prime


def approxFz(R1,R2,omega1,omega_I,Omega):
    #Approximate Fz from Eq. (51)
    w_eff = np.sqrt(omega1**2 + Omega**2)
    sin = omega1/w_eff
    cos = Omega/w_eff
    sin2 = sin**2
    cos2 = cos**2
    
    R1til = (R1*cos2 + R2*sin2) 
    R2til = (R2*(1+cos2) + R1*sin2)/2
    
    lambda_0 = R1til + 1j*omega_I
    lambda_m = R2til + 1j*(omega_I - w_eff)
    lambda_p = R2til + 1j*(omega_I + w_eff)
    
    Fz = cos2/lambda_0 + 0.5*sin2*(1/lambda_m + 1/lambda_p)
    return Fz


def TiTz_exact(R1,R2,omega1,omega_I,Omega):
    #New Bloch
    Fx, Fy, Fz, Fyp = new_bloch_ss(R1,R2,omega1,omega_I,Omega)
    #Bloch
    fx, fy, fz = bloch_ss(R1,R2,omega1,Omega)

    #
    Ti = np.real(Fz)

    Fz0 = 1/(R1+1j*omega_I)
    Ti0 = np.real(Fz0)
    vp_d2 = Ti - Ti0
    
    Txp = omega1 * np.real(1j*Fy*Fz)
    Typ = omega1 * np.real(-1j*Fx*Fz)

    Tx = omega1 * np.real(1j*Fyp*Fz)
    Tz = (omega1*fx)*Tx
        
    return Ti, Ti0, vp_d2, Tz, Txp, Typ, Tx#, vp_d2


def vZvD_classical(R1,R2,omega1,omega_I,Omega):

    ReLp = np.real(1/(R2 +1j*(omega_I + Omega)))
    ReLm = np.real(1/(R2 +1j*(omega_I - Omega)))
    
    vZ = 1/2 * (omega1/omega_I)**2 * ReLp
    vD = 1/2 * (omega1/omega_I)**2 * ReLm
    return vZ, vD

########################
#######################

def lambdas2model(lambdas,model,tau):
    j11s = []

    for lam in lambdas:
        if model == 'solid':
            j11s.append(1/lam)
        elif model == 'exp':
            j = tau*1/(1+lam*tau)
            j11s.append(j)
        elif model == 'ffhs':
            x = np.sqrt(lam*tau)
            j = tau*(x+4)/(x**3 + 4*x**2 + 9*x +9)
            j11s.append(j)

    return j11s

############################
def eig_numerical(R1,R2,omega1,omega_I,Omega,model,tau):

    h001 = np.array([[0,0,1]])
    v001 = h001.T
    v100 = 1j*np.array([[-1],[0],[0]])
    v010 = 1j*np.array([[0],[1],[0]])

    iFz0 = R1 + 1j*omega_I

    Ti = []
    vp  = []
    Txp = []
    Typ = []

    for i,W in enumerate(Omega):
    
        mxB = np.array([
            [R2 + 1j*omega_I,      W,               0],
            [-W,       R2 + 1j*omega_I,          omega1],
            [0,                         -omega1, R1 + 1j*omega_I]
        ])

        val,vec = np.linalg.eig(mxB)
        
        lambdas=[iFz0, val[0], val[1], val[2]]
        j_iFz0, j_0, j_1, j_2 = lambdas2model(lambdas,model,tau)
        
        Ti0 = np.real(j_iFz0)
        
        val_d = np.diag(np.array([j_0, j_1, j_2]))
        
        i_vec = np.linalg.inv(vec)

        iBnew = vec @ val_d @ i_vec
        
        
        ti = np.real(h001 @ iBnew @ v001 ).flatten()[0]
        tx = np.real(h001 @ iBnew @ v010 ).flatten()[0]
        ty = np.real(h001 @ iBnew @ v100 ).flatten()[0]
        
        #Re_Cs = np.real( horiz @ val_d @ ((vert_Sx * Sx[i]) + vert_iSy * iSy[i]) )
        
        Ti.append(ti)
        vp.append(ti - Ti0)
        Txp.append(tx)
        Typ.append(ty)

    return np.array(Ti), Ti0, np.array(vp), np.array(Txp), np.array(Typ)


############################
def eig_perturbative(R1,R2,omega1,omega_I,Omega,model,tau):

    w_eff = np.sqrt(omega1**2 + Omega**2)
    sin = omega1/w_eff
    cos = Omega/w_eff
    sin2 = sin**2
    cos2 = cos**2
    
    R1til = (R1*cos2 + R2*sin2) 
    R2til = (R2*(1+cos2) + R1*sin2)/2
    
    iFz0     = R1    + 1j*omega_I
    lambda_0 = R1til + 1j*omega_I
    lambda_m = R2til + 1j*(omega_I - w_eff)
    lambda_p = R2til + 1j*(omega_I + w_eff)

    
    lambdas=[iFz0,lambda_0, lambda_m, lambda_p]
    j_iFz0, j_lambda_0, j_lambda_m, j_lambda_p = lambdas2model(lambdas,model,tau)

    j_0  = j_lambda_0
    j_sum_hf = (j_lambda_m + j_lambda_p)/2
    j_dif_hf = (j_lambda_m - j_lambda_p)/2
    
    Ti0 = np.real(j_iFz0)
    Ti  = np.real(j_0*cos2 + j_sum_hf*sin2)
    vp = Ti - Ti0

    X = R2*(j_0 - j_sum_hf) -1j*w_eff*j_dif_hf
    Txp = omega1/(w_eff**2)*np.real(1j*X)

    return Ti, Ti0, vp, Txp


########################
#######################

def process_2columns(lines,scalings=[1,1]):
    data = []
    
    for line in lines:
        x, y = line[:-1].split('\t')
        data.append([ scalings[0] * float(x.replace(',', '.')), 
                     scalings[1] * float(y.replace(',', '.')) ])
        
    return np.array(data)

################################
def field_to_frequency(fields_G,key):
    '''
    Inp: magntic field in Gauss
    Out: frequency in MHz
    '''
    if key == 'epr':
        center_G = 94185.01 # G for EPR spectrum 
    elif key == 'dnp':
        center_G = 94182.15 # G for DNP spectrum of BDPA

    D_fields_T = (fields_G - center_G)* 1e-4
    
    conversion_GHz_per_T = 13.99625 # GHz/T
    g = 2.003
    
    rescale_A_to_T = 1.055
    
    return  rescale_A_to_T * (g * conversion_GHz_per_T * D_fields_T) * 1e3


def Lorentzian_derivative(x,center,width):

    R = width
    R2 = R**2
    Omega = x - center
    O2 = Omega**2

    denum = np.pi*(R2 + O2)**2

    return -2*R*Omega/denum, -(R2-O2)/denum


def calculate_enhancement(model,frequencies,omega_I,B1,T2e,R1n_del2,r1,r2,r3):

    Omega = 2*np.pi*frequencies
    omega1 = omega1_from_B1(B1)
    R2 = 1/T2e
    R1 = R2/r1

    tau = T2e/r2

    fx, fy, fz = bloch_ss(R1,R2,omega1,Omega)
    p =  R1 * fz

    Ti,Ti0,vp_d2, Txp, Typ = eig_numerical(R1,R2,omega1,omega_I,Omega,model,tau)
    Tz = omega1*(fx*Txp - fy*Typ)
    
    return 658*p*Tz/(R1n_del2/r3+vp_d2)


def RMSDerror(y_true,y_pred):
    return np.sum(np.abs(y_true)**1.9*(y_pred-y_true)**2)
