# -*- coding: utf-8 -*-

# This subfunction defines mass of particles based on PIP measurement
# utilizing hydrodynamic method of MH Mitchell and Heymsfield, 2005 parametrizized by Szyrmer & Zawadzki 2010

def MassEstimate(n, D_PIP, time_vector_par, envr, PIPD_par, PIPV_par, PIPEmajmax, PIPAR, kD, diacorr, foldername):



    # Define the area ratio as with the minimum diameter and area as
    # circumscribing sphere
    Area_ratio_sph_max = PIPAR/(pi*(PIPEmajmax)^2)
    Area_ratio_sph_max(Area_ratio_sph_max>1) = 1
    
    # Check that velocities are larger than 0
    PIPV_par(PIPV_par<0) = np.NaN
    
    # Exclude all the NaN values
    PIPV_par(np.isnan(PIPEmajmax)) = []
    PIPD_par(np.isnan(PIPEmajmax)) = []
    Area_ratio_sph_max(np.isnan(PIPEmajmax)) = []
    time_vector_par(np.isnan(PIPEmajmax)) = []
    PIPEmajmax(np.isnan(PIPEmajmax)) = []
    
    PIPD_par(np.isnan(PIPV_par)) = []
    Area_ratio_sph_max(np.isnan(PIPV_par)) = []
    time_vector_par(np.isnan(PIPV_par)) = []
    PIPEmajmax(np.isnan(PIPV_par)) = []
    PIPV_par(np.isnan(PIPV_par)) = []
    
    PIPD_par(np.isnan(Area_ratio_sph_max)) = []
    time_vector_par(np.isnan(Area_ratio_sph_max)) = []
    PIPEmajmax(np.isnan(Area_ratio_sph_max)) = []
    PIPV_par(np.isnan(Area_ratio_sph_max)) = []
    Area_ratio_sph_max(np.isnan(Area_ratio_sph_max)) = []
    
    # Compute the mass
    # Define the environmental parameters 
    epsilon = 0.622# ratio of molecular weight of water and air
    Rcon = 8.3144*10^3# universal gas constant in cm^3kPa/(K*mole)
    Ma = 28.9644 # molecular weight of air in g/mole
    g = 981 #
            
    # Environment parameters
    nu = 1.8325e-4*((296.16+120)/(envr.temp+273.16+120))*((envr.temp+273.16)/(296.16))^(3/2)# dynamic viscosity [g/cms]
    # CALCULATE THE VAPOR DENSITY AT WATER SATURATION AT TEMP OF ENVIRONMENT
    # Goff-Gratch formulation
    a = 373.16/(envr.temp + 273.16)
    b1 = 7.90298*(a-1)
    b2 = 5.02808*log10(a)
    b3 = (1.3816e-7)*(10^(11.344*(1-1/a))-1)
    b4 = 8.1328e-3*(10^(-(3.49149*(a-1)))-1)
    ew = 10^(b2-b3+b4+3.0057166-b1) #saturation vapor pressure over water ar temperature T, kPa
    rhoa = Ma*(0.1*envr.press-(1-epsilon)*envr.rh*ew)/(Rcon*(envr.temp+273.16)) #density of moist air in g/cm^3 
    # Interpolate the environmental values for the calculation        
    nu_int = interp1(envr.time,nu,time_vector_par)
    rhoa_int = interp1(envr.time,rhoa,time_vector_par) 
    nu_int = double(nu_int)
    rhoa_int = double(rhoa_int)
    C0 = 0.6 # inviscid drag coeffient 
    delta0 = 5.83 # a constant related to thickness of boundary layer
    Rey_maxmaxD = (2*rhoa_int*PIPEmajmax*0.1*PIPV_par*100)/nu_int
    # constant correction factor for aspect ratio 0.8
    corr_fac_const08 = 1.11
    Rey_maxmaxD_corrconst08 = (2*rhoa_int*corr_fac_const08*PIPEmajmax*0.1*PIPV_par*100)/nu_int
    # constant correction factor for aspect ratio 0.6
    corr_fac_const06 = 1.22
    Rey_maxmaxD_corrconst06 = (2*rhoa_int*corr_fac_const06*PIPEmajmax*0.1*PIPV_par*100)/nu_int
    # constant correction factor for aspect ratio 0.4
    corr_fac_const04 = 1.42
    Rey_maxmaxD_corrconst04 = (2*rhoa_int*corr_fac_const04*PIPEmajmax*0.1*PIPV_par*100)/nu_int
    ############################################################################
    # MH05 as parametrized with Szyrmer and Zawadzki 2010
    X_maxmaxD_MH05 = 10^(3.8233*log10(Rey_maxmaxD)-1.5211*(log10(Rey_maxmaxD)^2)+0.30065*(log10(Rey_maxmaxD)^3)
        -0.06104*(log10(Rey_maxmaxD)^4)+0.13074*(log10(Rey_maxmaxD)^5)-0.073429*(log10(Rey_maxmaxD)^6)
        +0.016006*(log10(Rey_maxmaxD)^7)-0.0012483*(log10(Rey_maxmaxD)^8))
    m_maxmaxD_MH05 = (pi*nu_int^2*X_maxmaxD_MH05)/(8*g*rhoa_int)*Area_ratio_sph_max^(1/4)
    # MH05 as parametrized with Szyrmer and Zawadzki 2010 with constant
    # correction for aspect ratio of 0.8
    X_maxmaxD_MH05_corrconst08 = 10^(3.8233*log10(Rey_maxmaxD_corrconst08)-1.5211*(log10(Rey_maxmaxD_corrconst08)^2)+0.30065*(log10(Rey_maxmaxD_corrconst08)^3)
        -0.06104*(log10(Rey_maxmaxD_corrconst08)^4)+0.13074*(log10(Rey_maxmaxD_corrconst08)^5)-0.073429*(log10(Rey_maxmaxD_corrconst08)^6) 
        +0.016006*(log10(Rey_maxmaxD_corrconst08)^7)-0.0012483*(log10(Rey_maxmaxD_corrconst08)^8))
    m_maxmaxD_MH05_corrconst08 = (pi*nu_int^2*X_maxmaxD_MH05_corrconst08)/(8*g*rhoa_int)*Area_ratio_sph_max^(1/4)
    # MH05 as parametrized with Szyrmer and Zawadzki 2010 with constant
    # correction for aspect ratio of 0.6
    X_maxmaxD_MH05_corrconst06 = 10^(3.8233*log10(Rey_maxmaxD_corrconst06)-1.5211*(log10(Rey_maxmaxD_corrconst06)^2)+0.30065*(log10(Rey_maxmaxD_corrconst06)^3) 
        -0.06104*(log10(Rey_maxmaxD_corrconst06)^4)+0.13074*(log10(Rey_maxmaxD_corrconst06)^5)-0.073429*(log10(Rey_maxmaxD_corrconst06)^6) 
        +0.016006*(log10(Rey_maxmaxD_corrconst06)^7)-0.0012483*(log10(Rey_maxmaxD_corrconst06)^8))
    m_maxmaxD_MH05_corrconst06 = (pi*nu_int^2*X_maxmaxD_MH05_corrconst06)/(8*g*rhoa_int)*Area_ratio_sph_max^(1/4)
    # MH05 as parametrized with Szyrmer and Zawadzki 2010 with constant
    # correction for aspect ratio of 0.4
    X_maxmaxD_MH05_corrconst04 = 10^(3.8233*log10(Rey_maxmaxD_corrconst04)-1.5211*(log10(Rey_maxmaxD_corrconst04)^2)+0.30065*(log10(Rey_maxmaxD_corrconst04)^3) 
        -0.06104*(log10(Rey_maxmaxD_corrconst04)^4)+0.13074*(log10(Rey_maxmaxD_corrconst04)^5)-0.073429*(log10(Rey_maxmaxD_corrconst04)^6) 
        +0.016006*(log10(Rey_maxmaxD_corrconst04)^7)-0.0012483*(log10(Rey_maxmaxD_corrconst04)^8))
    m_maxmaxD_MH05_corrconst04 = (pi*nu_int^2*X_maxmaxD_MH05_corrconst04)/(8*g*rhoa_int)*Area_ratio_sph_max^(1/4)
    
    time_massrel = time_vector_par
    
    PIPD = PIPD_par # the correction is devided before MassRel
    PIPV = PIPV_par
    diacorr_rel = 1
    flag.massrel_facname = 'MH05_maxmaxD'
    [time_mass.MH05_maxmaxD, amass_PIP.MH05_maxmaxD, bmass_PIP.MH05_maxmaxD]  =  MassRel_GPM(n, m_maxmaxD_MH05, PIPD, PIPV,time_massrel, D_PIP, flag, kD, diacorr_rel, foldername)
    
    diacorr_rel =  1.11
    flag.massrel_facname = 'MH05_maxmaxD_corrconst08' 
    [time_mass.MH05_maxmaxD_corrconst08, amass_PIP.MH05_maxmaxD_corrconst08, bmass_PIP.MH05_maxmaxD_corrconst08]  =  MassRel_GPM(n, m_maxmaxD_MH05_corrconst08, PIPD, PIPV,time_massrel, D_PIP, flag, kD, diacorr_rel, foldername)
    
    diacorr_rel =  1.22
    flag.massrel_facname = 'MH05_maxmaxD_corrconst06'
    [time_mass.MH05_maxmaxD_corrconst06, amass_PIP.MH05_maxmaxD_corrconst06, bmass_PIP.MH05_maxmaxD_corrconst06]  =  MassRel_GPM(n, m_maxmaxD_MH05_corrconst06, PIPD, PIPV,time_massrel, D_PIP, flag, kD, diacorr_rel, foldername)
    
    diacorr_rel =  1.42
    flag.massrel_facname = 'MH05_maxmaxD_corrconst04'
    [time_mass.MH05_maxmaxD_corrconst04, amass_PIP.MH05_maxmaxD_corrconst04, bmass_PIP.MH05_maxmaxD_corrconst04]  =  MassRel_GPM(n, m_maxmaxD_MH05_corrconst04, PIPD, PIPV,time_massrel, D_PIP, flag, kD, diacorr_rel, foldername)
    
    return time_mass, amass_PIP, bmass_PIP
