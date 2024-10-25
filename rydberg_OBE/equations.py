# Contains different differential equations we use in this project



import numpy as np

def rabi_oscillation(t, y, Omega, omega0, omega):
    #Rabi oscillation equations
    #Omega is the Rabi frequency, omega0 is the transition frequency, omega is the laser frequency
    delta = omega - omega0
    c1=y[0]
    c2=y[1]
    dc1dt=-1j/2*Omega*np.exp(1j*delta*t)*c2
    dc2dt=-1j/2*Omega*np.exp(-1j*delta*t)*c1
    dydt=[dc1dt,dc2dt]
    return dydt

def optical_bloch_eqs(t,y,Omega, Gamma, Delta, gamma_perp):
    #Two level optical Bloch equations
    #Omega rabi frequency, omega laser frequency, Gamma decay rate, 
    # Delta detuning, gamma_perp dephasing rate
    
    
    
    rho_gg = y[0]
    rrho_ge = y[1] # rho_ge in the rotating frame
    rrho_eg = y[2] # rho_eg in the rotating frame
    rho_ee = y[3]

    drho_eedt = 1j*Omega/2*(rrho_eg - rrho_ge) - Gamma*rho_ee
    drho_ggdt = -1j*Omega/2*(rrho_eg- rrho_ge) + Gamma*rho_ee
    drrho_gedt = -(gamma_perp + 1j*Delta)*rrho_ge  - 1j*Omega/2*(rho_ee-rho_gg)
    drrho_egdt = -(gamma_perp - 1j*Delta)*rrho_eg  + 1j*Omega/2*(rho_ee-rho_gg)
    dydt = [drho_ggdt, drrho_gedt, drrho_egdt, drho_eedt]
    return dydt

def three_level_OBEs(t,y,gamma2,Omega2,gamma1,Omega1,b1,b2,delta1,delta2):
    # Write down the three level optical bloch equations
    # gamma1/2 decay rate
    # Omega1/2 rabi frequency
    # b1/2 branching ratio
    # delta1/2 detuning
    # under the rotating wave approximation, so results independent of omega

    [rho_gg, rho_ii, rho_ee, rho_gi, rho_ig, 
     rho_ge, rho_eg, rho_ie, rho_ei] = y
    # Population equations
    drho_eedt = -gamma2*rho_ee +1j/2*Omega2*(rho_ei-rho_ie)
    drho_iidt = -gamma1*rho_ii +b2*gamma2*rho_ee-1j/2*Omega2*(rho_ei-rho_ie)+1j/2*Omega1*(rho_ig-rho_gi)
    drho_ggdt = b1*gamma1*rho_ii-1j/2*Omega1*(rho_ig-rho_gi)

    # Coherence equations
    drho_eidt = (-(gamma1+gamma2)/2+1j*delta2)*rho_ei+1j/2*Omega2*(rho_ee-rho_ii)+1j/2*Omega1*rho_eg
    drho_egdt = (-gamma2/2 + 1j*(delta1+delta2))*rho_eg+1j/2*Omega1*rho_ei-1j/2*Omega2*rho_ig
    drho_igdt = (-gamma1/2 + 1j*delta1)*rho_ig+1j/2*Omega1*(rho_ii-rho_gg)-1j/2*Omega2*rho_eg

    # Complex conjugates
    drho_gedt = np.conj(drho_egdt)
    drho_gidt = np.conj(drho_igdt)
    drho_iedt = np.conj(drho_eidt)
 

    dydt = [drho_ggdt, drho_iidt, drho_eedt, drho_gidt, drho_igdt, drho_gedt, drho_egdt, drho_iedt, drho_eidt]

    return dydt