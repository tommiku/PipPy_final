# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 08:54:01 2019

@author: kumlin


This is the main function for analyzing PIP data in Hyyti�l�
for BAECC 2014 and consecutive winters

"""
import os
import ConfigParser
import numpy as np
import shelve
import traceback
from datetime import datetime, timedelta


def shelving(my_shelf):
    for key in dir():
        try:
            my_shelf[key] = globals()[key]
        except TypeError:
            #
            # __builtins__, my_shelf, and imported modules can not be shelved.
            #
        #XARRAY -> netcdf
            print('ERROR shelving: {0}'.format(key))
        
def main():
    #Definitions of the event fetched from config.cfg
    config = ConfigParser.RawConfigParser()
    config.read("config.cfg")
    
    startdate = config.get("input",   "startdate")
    enddate = config.get("input",   "enddate")
    n = config.getint("input",   "n")
    diacorr = config.getfloat("input",   "diacorr")
    
    #Convert dates to datetime
    try:
        startdt = datetime.strptime(startdate, "%Y%m%d%H%M")
        enddt = datetime.strptime(enddate, "%Y%m%d%H%M")
    except Exception:
        traceback.print_exc()
    
    foldername = os.path.join(config.get("input",   "dataloc"), 
                              datetime.strftime(startdt, '%Y%m%d'))
    if not os.path.exists(foldername):                          
        os.mkdir(foldername)
    
    paraname =  os.path.join(foldername, 
                             datetime.strftime(startdt, '%H%M') + '_' + datetime.strftime(enddt, '%H%M'))
    print(paraname)                         
    parafile = shelve.open(paraname)
                             
    #Flags for optional data
    pluvio200 = config.getboolean("flags",   "pluvio200")
    pluvio400 = config.getboolean("flags",   "pluvio400")
    wind_GILL = config.getboolean("flags",   "wind_GILL")
    FMI_met = config.getboolean("flags",   "FMI_met")
    PIP_psd = config.getboolean("flags",   "PIP_psd")
    PIP_vel = config.getboolean("flags",   "PIP_vel")
    PIP_par = config.getboolean("flags",   "PIP_par")
    PIP_parrel = config.getboolean("flags",   "PIP_parrel")
    PIP_mass = config.getboolean("flags",   "PIP_mass")
    PIP_minute = config.getboolean("flags",   "PIP_minute")
    imag = config.getboolean("flags",   "imag")
    version = config.getboolean("flags",   "version")
    
    #Read precipitation data
    if (pluvio200 == 1 or pluvio400 == 1):
        print("Reading precipitation data (PLUVIO)")
        tv_PL200, PL200_acc, PL200_rr, tv_PL400, PL400_acc, PL400_rr = PluvioIntensity_GPM(n, startdt, enddt)
    else:
        tv_PL200 = []
        PL200_acc = []
        PL200_rr = []
        tv_PL400 = []
        PL400_acc = []
        PL400_rr = []
    #Convert dates to datetime
    try:
        startdt = datetime.strptime(startdate, "%Y%m%d%H%M")
        enddt = datetime.strptime(enddate, "%Y%m%d%H%M")
    except Exception:
        traceback.print_exc()
    
    #Save parameters to file for later examination
    ##parafile.write(tv_PL200, ' ' , PL200_acc, ' ' , PL200_rr , ' ' , tv_PL400 , ' ' , PL400_acc , ' ' , PL400_rr )
    #shelving(parafile)
    
    #Read data from FMI-met
    if (FMI_met == 1):
        print("Reading FMI MET data")
        tv_FMI, temp_FMI, rh_FMI, sn_FMI, rr_FMI, press_FMI = Read_FMIstation_GPM(startdt, enddt, n)
    else:
        tv_FMI = []
        temp_FMI = []
        rh_FMI = []
        sn_FMI = []
        rr_FMI = []
        press_FMI = []
        
    #parafile.write(tv_FMI , ' ' , temp_FMI, ' ' , rh_FMI, ' ' , sn_FMI, ' ' , rr_FMI, ' ' , press_FMI)
    #shelving(parafile)
    
    #Read and plot GILL wind data
    ##TODO: Sanity check
    if (wind_GILL == 1):
        print("Reading wind data (GILL)")
        time_vector_GILL, mean_vel_GILL, mode_dir_GILL =  WindComparison_GPM(n, startdt, enddt)
    else:
        time_vector_GILL = []
        mean_vel_GILL = []
        mode_dir_GILL = []
        
    #parafile.write(time_vector_GILL, ' ', mean_vel_GILL, ' ', mode_dir_GILL)
    #shelving(parafile)
    
    #Read the PSD tables (size distribution)
    if (PIP_psd == 1):
        print("Reading PSD tables")
        D_PIP, PIP_PSD, PIPtime_psd, N_mean_PIP, Dm-PIP, N0_PIP, lambda_PIP, 
        mu_PIP, Nw_PIP, D02_exp_PIP, lambda2_exp_PIP, PIPtime_N = ReadPIP_PSD_GPM(
        startdt, enddt, n, foldername)
    else:
        D_PIP = []
        PIP_PSD = []
        PIPtime_psd = []
        N_mean_PIP = [] 
        Dm_PIP = []
        N0_PIP = []
        lambda_PIP = []
        mu_PIP = []
        Nw_PIP = []
        D02_exp_PIP = []
        lambda2_exp_PIP = [] 
        PIPtime_N = []
        
    #parafile.write(D_PIP, ' ', PIP_PSD, ' ', PIPtime_psd, ' ',N_mean_PIP , ' ',Dm_PIP , ' ',N0_PIP , ' ', 
    #                ' ',lambda_PIP , ' ',mu_PIP , ' ',Nw_PIP , ' ',D02_exp_PIP , ' ',lambda2_exp_PIP  , ' ',PIPtime_N)
    #shelving(parafile)
    
    #Read the velocity tables and perform a fit as function of D_PIP
    if (PIP_vel == 1):
        print("Reading the velocity tables (PIP_vel)")
        D_PIP, PIPtime_vel, PIPD_vel, PIPV_vel = ReadPIP_vel_GPM(startdt, enddt)
        avel_DPIP_vel, bvel_DPIP_vel, PIPtime_vel_n = PIPVelRel_GPM(startdt, enddt, D_PIP, PIPtime_vel, PIPD_vel, PIPV_vel, n, foldername)
    else:
        PIPtime_vel = []
        PIPD_vel = []
        PIPV_vel = []
        avel_DPIP = []
        bvel_DPIP = []
        DPIP_V_vel_n = []
        PIPtime_vel_n = []
        
    #parafile.write(D_PIP, ' ', PIPtime_vel, ' ', PIPD_vel, ' ', PIPV_vel, ' ', avel_DPIP_vel, ' ', bvel_DPIP_vel, ' ', PIPtime_vel_n)
    #shelving(parafile)
    
    #Read the particle tables
    if (PIP_par == 1):
        D_PIP, PIPD_par, PIPV_par, PIPtime_par, PIPEmaj, PIPEmajmax, PIPEmin, PIPAR, PIPOR, PIPLen, PIPHig = ReadPIP_par_GPM(startdt, enddt)
        temp_PIPD = PIPD_par
        temp_Dmax = 2*PIPEmajmax
        temp_PIPD[np.isnan(tempDmax)] = []
        temp_Dmax[np.isnan(tempDmax)] = []
        temp_PIPD[np.isinf(tempDmax)] = []
        temp_Dmax[np.isinf(tempDmax)] = []
        
        temp_Dmax[np.isnan(temp_PIPD)] = []
        temp_PIPD[np.isnan(temp_PIPD)] = []
        temp_Dmax[np.isinf(temp_PIPD)] = []
        temp_PIPD[np.isinf(temp_PIPD)] = []
        
        #Area equivalent of whole snow event
        ##TODO
        #C1 =
        #kD = C1[1]
        
    else:
        D_PIP = []
        PIPD_par = []
        PIPV_par = []
        PIPtime_par = []
        PIPEmaj = []
        PIPEmajmax = []
        PIPEmin = []
        PIPAR = []
        PIPlonX = []
        PIPDia = []
        PIPOR = []
        PIPLen = []
        PIPHig = []
        kD = []
        
    #parafile.write(D_PIP, ' ', PIPtime_vel, ' ', PIPD_vel, ' ', PIPV_vel, ' ', avel_DPIP_vel, ' ', bvel_DPIP_vel, ' ', PIPtime_vel_n)
    #shelving(parafile)
    
    if (PIP_mass == 1):
        temp_FMI[np.isnan(rh_FMI)] = []
        press_FMI[np.isnan(rh_FMI)] = []
        tv_FMI[np.isnan(rh_FMI)] = []
        rh_FMI[np.isnan(rh_FMI)] = []
        press_FMI[np.isnan(temp_FMI)] = []
        tv_FMI[np.isnan(temp_FMI)] = []
        rh_FMI[np.isnan(temp_FMI)] = []
        temp_FMI[np.isnan(temp_FMI)] = [] 
        tv_FMI[np.isnan(press_FMI)] = []
        rh_FMI[np.isnan(press_FMI)] = []
        temp_FMI[np.isnan(press_FMI)] = []
        press_FMI[np.isnan(press_FMI)] = []
        envr.temp = temp_FMI
        envr.press = press_FMI
        envr.time = tv_FMI
        envr.rh = rh_FMI/100
        
    
        time_mass, amass_PIP, bmass_PIP = MassEstimate(n,D_PIP,PIPtime_par, envr, PIPD_par, PIPV_par, PIPEmajmax, PIPAR, kD, diacorr, foldername)
    else:
        time_mass = []
        amass_PIP = []
        bmass_PIP = []
    
    
    
    
    """
    Part 2
    """
    
    #Calculate the accumulation and reflectivity from the mass estimate
    master_time_vector = np.arange(startdt, enddt, timedelta(minutes=n))
    
    #PSD properties
    N0_mtv = []
    lambda_mtv = []
    mu_mtv = []
    
    N02_exp_mtv = []
    lambda2_exp_mtv = []
    D02_exp_mtv = []
    
    N_mean_PIP_mtv = []
    Dmax_mtv = []
    
    #Selected a, b (mass factors) and av, bv (velocity factors)
    amass_mtv = []
    bmass_mtv = []
    avel_mtv_max = []
    bvel_mtv_max = []
    """
    #Choose used diameter correction
    if diacorr == 0.9:
        time = time_mass.MH05_maxmaxD_corrconst08
        amass = amass_PIP.MH05_maxmaxD_corrconst08
        bmass = bmass_PIP.MH05_maxmaxD_corrconst08
    elif diacorr == 0.82:
        time = time_mass.MH05_maxmaxD_corrconst06
        amass = amass_PIP.MH05_maxmaxD_corrconst06
        bmass = bmass_PIP.MH05_maxmaxD_corrconst06
    elif diacorr == 0.7:
        time = time_mass.MH05_maxmaxD_corrconst04
        amass = amass_PIP.MH05_maxmaxD_corrconst04
        bmass = bmass_PIP.MH05_maxmaxD_corrconst04
    else:
        time_mtv = time_mass.MH05_maxmaxD
        amass = amass_PIP.MH05_maxmaxD
        bmass = bmass_PIP.MH05_maxmaxD
    """
    accum_PIP_fac_MH05_maxmaxD = []
    accum_PIP_fac_MH05_maxmaxD_corrconst08 = []
    accum_PIP_fac_MH05_maxmaxD_corrconst06 = []
    accum_PIP_fac_MH05_maxmaxD_corrconst04 = []
            
    Ze_PIP_fac_MH05_maxmaxD = []
    Ze_PIP_fac_MH05_maxmaxD_corrconst08 = []
    Ze_PIP_fac_MH05_maxmaxD_corrconst06 = []
    Ze_PIP_fac_MH05_maxmaxD_corrconst04 = []
    
    diff_PIP = np.zeros(np.shape(D_PIP))
    diff_PIP[0] = D_PIP[0]
    diff_PIP[1:] = np.diff(D_PIP)
    
    diacorr_no = kD/1
    diacorr08 = kD/0.9
    diacorr06 = kD/0.82
    diacorr04 = kD/0.70
    
    for dd in (np.arange(1, np.shape(master_time_vector)[1])):
        d_ind = np.where(time_mass.MH05_maxmaxD > (master_time_vector[dd]-30*(1/(24*3600))) and time_mass.MH05_maxmaxD < (master_time_vector[dd]+30*(1/(24*3600))))        
        #d_ind = np.where(time_mass.MH05_maxmaxD > (master_time_vector[dd]-30*(1/(24*3600))) & time_mass.MH05_maxmaxD < (master_time_vector[dd]+30*(1/(24*3600))));
        dn_ind = np.where(PIPtime_n > (master_time_vector[dd]-30*(1/24*3600)) and PIPtime_n < (master_time_vector[dd]+30*(1/(24*3600))))
        #dn_ind = np.where(PIPtime_n  > (master_time_vector(dd)-30*(1/(24*3600))) & PIPtime_n < (master_time_vector(dd)+30*(1/(24*3600))));
        dv_ind = np.where(time_vector_parrel_max > (master_time_vector[dd]-30*(1/(24*3600)) and time_vector_parrel_max < (master_time_vector[dd]+30*(1/(24*3600)))))        
        #dv_ind = np.where(time_vector_parrel_max > (master_time_vector(dd)-30*(1/(24*3600))) & time_vector_parrel_max < (master_time_vector(dd)+30*(1/(24*3600))));        
        if (np.isempty(d_ind) == 0  and np.isempty(dn_ind)==0 and np.isempty(dv_ind)==0):
            #PSD parameters
            N0_mtv = np.column_stack(N0_mtv, N0_PIP(dn_ind))
            lambda_mtv = np.column_stack(lambda_mtv, lambda_PIP(dn_ind))
            mu_mtv = np.column_stack(mu_mtv, mu_PIP(dn_ind))
                
            N02_exp_mtv = np.column_stack(N02_exp_mtv, N02_exp_PIP[dn_ind])
            lambda2_exp_mtv = np.column_stack(lambda2_exp_mtv, lambda2_exp_PIP[dn_ind])
            D02_exp_mtv = np.column_stack(D02_exp_mtv, D02_exp_PIP[dn_ind])
            N_mean_PIP_mtv = np.row_stack(N_mean_PIP_mtv, N_mean_PIP[dn_ind,:])
            Dmax_mtv = np.column_stack(Dmax_mtv, Dmax_PIP(dn_ind))
            
            #particle relations
            avel_mtv_max = np.column_stack(avel_mtv_max, avel_PIP_par_max(dv_ind))
            bvel_mtv_max = np.column_stack(bvel_mtv_max, bvel_PIP_par_max(dv_ind))
            
            accum_PIP_fac_MH05_maxmaxD = np.column_stack(accum_PIP_fac_MH05_maxmaxD, 
                                                         n*60*10^(-3)*np.nansum(amass_PIP.MH05_maxmaxD[d_ind]*(diacorr_no*0.1*D_PIP)^(bmass_PIP.MH05_maxmaxD[d_ind])*avel_PIP_par_max[dv_ind]*(diacorr_no*D_PIP)^(bvel_PIP_par[dv_ind])*N_mean_PIP[dn_ind,:]*diff_PIP))
            if (bmass_PIP.MH05_maxmaxD(d_ind)>=1 and bmass_PIP.MH05_maxmaxD(d_ind)< 3.5) and (bvel_PIP_par_max(dv_ind)>=0):
                Ze_PIP_fac_MH05_maxmaxD = np.column_stack(Ze_PIP_fac_MH05_maxmaxD, 10^6*1.2076*0.2/0.93*(6/pi)^2*nansum((amass_PIP.MH05_maxmaxD(d_ind)*(diacorr_no*0.1*D_PIP)^(bmass_PIP.MH05_maxmaxD(d_ind)))^2*N_mean_PIP[dn_ind,:]*diff_PIP))
                amass_mtv = np.column_stack(amass_mtv, amass(d_ind))
                bmass_mtv = np.column_stack(bmass_mtv, bmass(d_ind))       
            else:
                Ze_PIP_fac_MH05_maxmaxD = np.column_stack(Ze_PIP_fac_MH05_maxmaxD, 0)
                amass_mtv = np.column_stack(amass_mtv, 0)
                bmass_mtv = np.column_stack(bmass_mtv, 0)
        else:
            accum_PIP_fac_MH05_maxmaxD =  np.column_stack(accum_PIP_fac_MH05_maxmaxD, 0)
            Ze_PIP_fac_MH05_maxmaxD = np.column_stack(Ze_PIP_fac_MH05_maxmaxD, 0)
            amass_mtv = np.column_stack(mass_mtv, 0)
            bmass_mtv = np.column_stack(bmass_mtv, 0)
            N0_mtv = np.column_stack(N0_mtv, 0)
            lambda_mtv = np.column_stack(lambda_mtv, 0)
            mu_mtv = np.column_stack(mu_mtv, 0)
            N02_exp_mtv = np.column_stack(N02_exp_mtv, 0)
            D02_exp_mtv = np.column_stack(D02_exp_mtv, 0)
            lambda2_exp_mtv = np.column_stack(lambda2_exp_mtv, 0)
            N_mean_PIP_mtv = np.row_stack(N_mean_PIP_mtv, np.zeros(size(D_PIP)))
            Dmax_mtv = np.column_stack(Dmax_mtv, 0)
            avel_mtv_max = np.column_stack(avel_mtv_max, 0)
            bvel_mtv_max = np.column_stack(bvel_mtv_max, 0)
       
        d_ind = np.where(time_mass.MH05_maxmaxD_corrconst08 > (master_time_vector[dd]-30*(1/(24*3600))) and time_mass.MH05_maxmaxD_corrconst08 < (master_time_vector[dd]+30*(1/(24*3600))))
        dn_ind = np.where(PIPtime_n  > (master_time_vector(dd)-30*(1/(24*3600))) & PIPtime_n < (master_time_vector(dd)+30*(1/(24*3600))));
        dv_ind = np.where(time_vector_parrel_max > (master_time_vector(dd)-30*(1/(24*3600))) & time_vector_parrel_max < (master_time_vector(dd)+30*(1/(24*3600))));
        if np.isempty(d_ind) == 0  and np.isempty(dn_ind)==0 and np.isempty(dv_ind)==0 :
           accum_PIP_fac_MH05_maxmaxD_corrconst08 = np.column_stack(accum_PIP_fac_MH05_maxmaxD_corrconst08, n*60*10^(-3)*np.nansum(amass_PIP.MH05_maxmaxD_corrconst08[d_ind]*(diacorr08*0.1*D_PIP)^(bmass_PIP.MH05_maxmaxD_corrconst08[d_ind])*avel_PIP_par_max[dv_ind]*(D_PIP*diacorr08)^(bvel_PIP_par_max[dv_ind])*N_mean_PIP[dn_ind,:]*diff_PIP))
           if (bmass_PIP.MH05_maxmaxD_corrconst08(d_ind)>=1 and bmass_PIP.MH05_maxmaxD_corrconst08(d_ind)< 3.5) and (bvel_PIP_par_max(dv_ind)>=0):
               Ze_PIP_fac_MH05_maxmaxD_corrconst08 = np.column_stack(Ze_PIP_fac_MH05_maxmaxD_corrconst08, 10^6*1.2076*0.2/0.93*(6/pi)^2*np.nansum((amass_PIP.MH05_maxmaxD_corrconst08[d_ind]*(diacorr08*0.1*D_PIP)^(bmass_PIP.MH05_maxmaxD_corrconst08[d_ind]))^2*N_mean_PIP[dn_ind,:]*diff_PIP))
           else:
               Ze_PIP_fac_MH05_maxmaxD_corrconst08 = np.column_stack(Ze_PIP_fac_MH05_maxmaxD_corrconst08, 0)
        else:
           accum_PIP_fac_MH05_maxmaxD_corrconst08 =  np.column_stack(accum_PIP_fac_MH05_maxmaxD_corrconst08, 0)
           Ze_PIP_fac_MH05_maxmaxD_corrconst08 = np.column_stack(Ze_PIP_fac_MH05_maxmaxD_corrconst08, 0)
           
        d_ind = np.where(time_mass.MH05_maxmaxD_corrconst06 > (master_time_vector[dd]-30*(1/(24*3600))) and time_mass.MH05_maxmaxD_corrconst06 < (master_time_vector[dd]+30*(1/(24*3600))))
        dn_ind = np.where(PIPtime_n  > (master_time_vector[dd]-30*(1/(24*3600))) and PIPtime_n < (master_time_vector[dd]+30*(1/(24*3600))))
        dv_ind = np.where(time_vector_parrel_max > (master_time_vector[dd]-30*(1/(24*3600))) and time_vector_parrel_max < (master_time_vector[dd]+30*(1/(24*3600))))
        if np.isempty(d_ind) == 0  and np.isempty(dn_ind)==0 and np.isempty(dv_ind)==0 :
           accum_PIP_fac_MH05_maxmaxD_corrconst06 = np.column_stack(accum_PIP_fac_MH05_maxmaxD_corrconst06, n*60*10^(-3)*np.nansum(amass_PIP.MH05_maxmaxD_corrconst06[d_ind]*(diacorr06*0.1*D_PIP)^(bmass_PIP.MH05_maxmaxD_corrconst06[d_ind])*avel_PIP_par_max[dv_ind]*(diacorr06*D_PIP)^(bvel_PIP_par_max[dv_ind])*N_mean_PIP[dn_ind,:]*diff_PIP))
           if (bmass_PIP.MH05_maxmaxD_corrconst06[d_ind]>=1 and bmass_PIP.MH05_maxmaxD_corrconst06[d_ind]< 3.5) and (bvel_PIP_par_max[dv_ind]>=0):
               Ze_PIP_fac_MH05_maxmaxD_corrconst06 = np.column_stack(Ze_PIP_fac_MH05_maxmaxD_corrconst06, 10^6*1.2076*0.2/0.93*(6/pi)^2*np.nansum((amass_PIP.MH05_maxmaxD_corrconst06[d_ind]*(diacorr06*0.1*D_PIP)^(bmass_PIP.MH05_maxmaxD_corrconst06[d_ind]))^2*N_mean_PIP[dn_ind,:]*diff_PIP))   
           else:
               Ze_PIP_fac_MH05_maxmaxD_corrconst06 = np.column_stack(Ze_PIP_fac_MH05_maxmaxD_corrconst06, 0)
        else:
           accum_PIP_fac_MH05_maxmaxD_corrconst06 =  np.column_stack(accum_PIP_fac_MH05_maxmaxD_corrconst06, 0)
           Ze_PIP_fac_MH05_maxmaxD_corrconst06 = np.column_stack(Ze_PIP_fac_MH05_maxmaxD_corrconst06, 0)
        d_ind = np.where(time_mass.MH05_maxmaxD_corrconst04 > (master_time_vector[dd]-30*(1/(24*3600))) and time_mass.MH05_maxmaxD_corrconst04 < (master_time_vector[dd]+30*(1/(24*3600))))
        dn_ind = np.where(PIPtime_n  > (master_time_vector[dd]-30*(1/(24*3600))) & PIPtime_n < (master_time_vector[dd]+30*(1/(24*3600))))
        dv_ind = np.where(time_vector_parrel_max > (master_time_vector[dd]-30*(1/(24*3600))) and time_vector_parrel_max < (master_time_vector[dd]+30*(1/(24*3600))))
        if (np.isempty(d_ind) == 0  and np.isempty(dn_ind)==0 and np.isempty(dv_ind)==0):
           accum_PIP_fac_MH05_maxmaxD_corrconst04 = np.column_stack(accum_PIP_fac_MH05_maxmaxD_corrconst04, n*60*10^(-3)*np.nansum(amass_PIP.MH05_maxmaxD_corrconst04[d_ind]*(diacorr04*0.1*D_PIP)^(bmass_PIP.MH05_maxmaxD_corrconst04[d_ind])*avel_PIP_par_max[dv_ind]*(diacorr04*D_PIP)^(bvel_PIP_par[dv_ind])*N_mean_PIP[dn_ind,:]*diff_PIP))
           if (bmass_PIP.MH05_maxmaxD_corrconst04[d_ind]>=1 and bmass_PIP.MH05_maxmaxD_corrconst04[d_ind]< 3.5) and (bvel_PIP_par_max[dv_ind]>=0):
               Ze_PIP_fac_MH05_maxmaxD_corrconst04 = np.column_stack(Ze_PIP_fac_MH05_maxmaxD_corrconst04, 10^6*1.2076*0.2/0.93*(6/pi)^2*np.nansum((amass_PIP.MH05_maxmaxD_corrconst04[d_ind]*(diacorr04*0.1*D_PIP)^(bmass_PIP.MH05_maxmaxD_corrconst04[d_ind]))^2*N_mean_PIP[dn_ind,:]*diff_PIP))
           else:
               Ze_PIP_fac_MH05_maxmaxD_corrconst04 = np.column_stack(Ze_PIP_fac_MH05_maxmaxD_corrconst04, 0)
        else:
           accum_PIP_fac_MH05_maxmaxD_corrconst04 =  np.column_stack(accum_PIP_fac_MH05_maxmaxD_corrconst04, 0)
           Ze_PIP_fac_MH05_maxmaxD_corrconst04 = np.column_stack(Ze_PIP_fac_MH05_maxmaxD_corrconst04, 0)
    
    #shelving(parafile)
    
    #Plot the summary of the event (IMPORTANT)
    Plot_EventSummary_GPM()
    
    #Define Z(S) with Rayleigh approximation
    #Choose the relation
    if (diacorr == 0.9):
        accum = accum_PIP_fac_MH05_maxmaxD_corrconst08
        Ze = Ze_PIP_fac_MH05_maxmaxD_corrconst08
        fname1 = os.path.join(foldername, '\ZeS_corrconst08_',startdt,'_',enddt)
        fname2 = os.path.join(foldername, '\ZeS_timeseries_corrconst08_',startdt,'_',enddt)
    elif (diacorr == 0.82):
        accum = accum_PIP_fac_MH05_maxmaxD_corrconst06
        Ze = Ze_PIP_fac_MH05_maxmaxD_corrconst06
        fname1 = os.path.join(foldername, '\ZeS_corrconst06_',datestr(event_start_time,30),'_',datestr(event_end_time,30))
        fname2 = os.path.join(foldername, '\ZeS_timeseries_corrconst06_',datestr(event_start_time,30),'_',datestr(event_end_time,30))
    elif (diacorr == 0.7):
        accum = accum_PIP_fac_MH05_maxmaxD_corrconst04
        Ze = Ze_PIP_fac_MH05_maxmaxD_corrconst04
        fname1 = (foldername, '\ZeS_corrconst04_',datestr(event_start_time,30),'_',datestr(event_end_time,30))
        fname2 = (foldername, '\ZeS_timeseries_corrconst04_',datestr(event_start_time,30),'_',datestr(event_end_time,30))
    else:
        accum = accum_PIP_fac_MH05_maxmaxD
        Ze = Ze_PIP_fac_MH05_maxmaxD
        fname1 = (foldername, '\ZeS_nocorr_',datestr(event_start_time,30),'_',datestr(event_end_time,30))
        fname2 = (foldername, '\ZeS_timeseries_nocorr_',datestr(event_start_time,30),'_',datestr(event_end_time,30))
        

    parafile.close()


if __name__ == "__main__":
    main()