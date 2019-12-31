# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 11:11:16 2019

@author: kumlin
"""
import numpy as np
import os
from datetime import datetime, timedelta
from glob import glob
# read FMI tables for temperature humidity and rain rate

def Read_FMIstation_GPM(event_start_time, event_end_time, n):

    event_start_day = event_start_time.replace(hour=0, minute=0)
    event_end_day = event_end_time.replace(hour=0, minute=0)
    
    date_t_FMI = []
    temp_FMI = []
    rh_FMI = []
    date_sn_FMI = []
    sn_temp_FMI = []
    sn_rh_FMI = []
    sn_FMI = []
    date_rr_FMI = []
    rr_FMI = []
    press_FMI_J = []
    date_press_FMI_J = []
    temp_FMI_tv = []
    rh_FMI_tv = [] 
    sn_FMI_tv = [] 
    rr_FMI_tv = []
    press_FMI_J_tv = []
    file_year = []
    file_month = []
    file_name= []
    
    """
    Read data into arrays for events later than 11.1.2017
    """
        
    # open FMI observations of the Hyyti�l� site
    #TODO: READ FROM CONFIG
    folder_FMI = '/home/kumlin/Koodit/SnowRetrievals/Data_malli_10062019/FMI_MET/'
    if (event_start_day < datetime(2016,04,30)):
        # search the temperature and relative humidity files
        #np.where FILENAME HERE
        flist_FMI = glob(os.path.join(folder_FMI,'FMI2303_20*_T_TD_TW_RH.csv'))
        # Select files that fall within start and end times
        file_name, indx1 = select_files(flist_FMI, event_start_day, event_end_day)
        jj = 0
    
        date_t_FMI = np.empty([1,5])
        #TÄMÄ OSA TARKISTETTU
        for jj in range(0,indx1):
            fname_FMI = os.path.join(folder_FMI,str(file_name[jj]))
            Data = np.loadtxt(fname_FMI, dtype=str, delimiter=',')
            #print(Data)
            date_tmp = str(Data[1:len(Data[2])-1,2])
            #print(date_tmp[2:3])
            day_FMI = str(date_tmp[2:4])
            month_FMI = str(date_tmp[5:7])
            year_FMI = str(date_tmp[8:12])
            time_tmp = str(Data[1:len(Data[3])-1,3])
            #print(time_tmp)
            hour_FMI = str(time_tmp[2:4])
            min_FMI = str(time_tmp[5:7])
            date_t_FMI = np.vstack((date_t_FMI, np.array([year_FMI , month_FMI , day_FMI , hour_FMI , min_FMI])))
            if temp_FMI in dir():
                
                temp_FMI = np.vstack((temp_FMI, Data[:,5]))
            else:
                
                temp_FMI = np.empty_like(Data[:,5])
                temp_FMI = np.vstack((temp_FMI, Data[:,5]))
            if rh_FMI in dir():
                rh_FMI = np.vstack((rh_FMI, Data[:,8]))
            else:
                rh_FMI = np.empty_like(Data[:,8])
                rh_FMI = np.vstack((rh_FMI, Data[:,8]))
            #print(date_t_FMI)
            #print(temp_FMI)
            #print(rh_FMI)
            
    
        # search the temperature and relative humidity with snowdepth files
        flist_FMI = glob(os.path.join(folder_FMI,'FMI2303_20*_T_RH_SN.csv'))
        #flist_FMI =  glob(os.path.join(folder_FMI,filename_FMI))
    
        # Select files that fall within start and end times
        file_name, indx1 = select_files(flist_FMI, event_start_day, event_end_day)
        jj = 0
    
        date_sn_FMI = np.empty([1,5])
        
        for jj in range(0,indx1):
            fname_FMI = os.path.join(folder_FMI,str(file_name[jj]))
            Data = np.loadtxt(fname_FMI, dtype=str, delimiter=',')
            #print(Data)
            date_tmp = str(Data[1:len(Data[2])-1,2])
            #print(date_tmp[2:3])
            day_FMI = str(date_tmp[2:4])
            month_FMI = str(date_tmp[5:7])
            year_FMI = str(date_tmp[8:12])
            time_tmp = str(Data[1:len(Data[3])-1,3])
            #print(time_tmp)
            hour_FMI = str(time_tmp[2:4])
            min_FMI = str(time_tmp[5:7])
            date_sn_FMI = np.vstack((date_sn_FMI, np.array([year_FMI , month_FMI , day_FMI , hour_FMI , min_FMI])))
            if sn_temp_FMI in dir():
                sn_temp_FMI = np.vstack((sn_temp_FMI, Data[:,5]))
            else:
                sn_temp_FMI = np.empty_like(Data[:,5])
                temp_FMI = np.vstack((sn_temp_FMI, Data[:,5]))
            if sn_rh_FMI in dir():
                sn_rh_FMI = np.vstack((sn_rh_FMI, Data[:,6]))
            else:
                sn_rh_FMI = np.empty_like(Data[:,6])
                sn_rh_FMI = np.vstack((sn_rh_FMI, Data[:,6]))
            if sn_FMI in dir():
                sn_FMI = np.vstack((sn_FMI, Data[:,11]))
            else:
                sn_FMI = np.empty_like(Data[:,11])
                sn_FMI = np.vstack((sn_FMI, Data[:,11]))
            #print(date_sn_FMI)
            #print(sn_temp_FMI)
            #print(sn_rh_FMI)
            
    
    
        # search the rainrate
        #np.where FILENAME HERE
        flist_FMI =  glob(os.path.join(folder_FMI,'FMI2303_20*_RI.csv'))
        # Select files that fall within start and end times


        
        file_name, indx1 = select_files(flist_FMI, event_start_day, event_end_day)
        jj = 0
        date_rr_FMI = np.empty([1,5])
    
        for jj in range(0,indx1):
            fname_FMI = os.path.join(folder_FMI,str(file_name[jj]))
            Data = np.loadtxt(fname_FMI, dtype=str, delimiter=',')
            #print(Data)
            date_tmp = str(Data[1:len(Data[2])-1,2])
            #print(date_tmp[2:3])
            day_FMI = str(date_tmp[2:4])
            month_FMI = str(date_tmp[5:7])
            year_FMI = str(date_tmp[8:12])
            time_tmp = str(Data[1:len(Data[3])-1,3])
            #print(time_tmp)
            hour_FMI = str(time_tmp[2:4])
            min_FMI = str(time_tmp[5:7])
            date_rr_FMI = np.vstack((date_rr_FMI, np.array([year_FMI , month_FMI , day_FMI , hour_FMI , min_FMI])))
            if rr_FMI in dir():
                print("Dir")
                rr_FMI = np.vstack((rr_FMI, Data[:,4]))
            else:
                print("Tämä vain kerran")
                rr_FMI = np.empty_like(Data[:,4])
                print(np.shape(rr_FMI))
                rr_FMI = np.vstack((rr_FMI, Data[:,4]))    
    
        
        # search the pressure data from J�ms�
        flist_FMI = flist_FMI = glob(os.path.join(folder_FMI,'FMI2301_20*_T_TD_TW_RH_P.csv'))
        # Select files that fall within start and end times
        file_name, indx1 = select_files(flist_FMI, event_start_day, event_end_day)
        date_press_FMI_J = np.empty([1,5])
        
        for jj in range(0,indx1):
            fname_FMI = os.path.join(folder_FMI,str(file_name[jj]))
            Data = np.loadtxt(fname_FMI, dtype=str, delimiter=',')
            print(np.shape(Data))
            #print(Data)
            date_tmp = str(Data[1:len(Data[2])-1,2])
            #print(date_tmp[2:3])
            day_FMI_J = str(date_tmp[2:4])
            month_FMI_J = str(date_tmp[5:7])
            year_FMI_J = str(date_tmp[8:12])
            time_tmp = str(Data[1:len(Data[3])-1,3])
            #print(time_tmp)
            hour_FMI_J = str(time_tmp[2:4])
            min_FMI_J = str(time_tmp[5:7])
            date_press_FMI_J = np.vstack((date_press_FMI_J, np.array([year_FMI_J , month_FMI_J , day_FMI_J , hour_FMI_J , min_FMI_J])))
            if press_FMI_J in dir():
                print("Dir")
                press_FMI_J = np.vstack((press_FMI_J, Data[:,9].astype(float)))
            else:
                print("Testitulostus_press")
                press_FMI_J = np.empty_like(Data[:,9])
                press_FMI_J = np.vstack((press_FMI_J, Data[:,9].astype(float)))
                print(press_FMI_J)
        
    '''
    elif (event_start_day > datetime.datetime(2016,04,30) and (event_start_day < datetime.datetime(2017,01,11)):
        
        # search the temperature and relative humidity, precipitation rate and snow depth files
        #np.where FILENAME HERE        
        flist_FMI = os.listdir(os.path.join(folder_FMI,'FMI2303_20*_T_TD_RH_RI_SD.csv'))
        findx = 0
        indx1 = 0
        jj = 0
        file_name = {}
        # Select files that fall within start and end times
        for (findx in range(0,len(flist_FMI))):
            filename = flist_FMI(findx,1).name
            file_year(findx) = int([filename(9:12)])
            file_month(findx) = int([filename(13:14)])
            if (file_month(findx) == int(datestr(event_start_day, 'mm')) and file_year(findx) == int(datestr(event_start_day, 'yyyy'))):
                indx1      = indx1+1
                file_name(indx1)= {flist_FMI(findx,1).name}
            elif (file_month(findx) == int(datestr(event_end_day, 'mm')) and file_year(findx) == int(datestr(event_end_day, 'yyyy'))):
                indx1      = indx1+1
                file_name(indx1)= {flist_FMI(findx,1).name}
    
        for (jj in range(0,indx1):
            fname_FMI = [folder_FMI,str(file_name(jj))]
            Data = importdata(fname_FMI)
            date_tmp = str( Data.textdata(2:end-1,3))
            day_FMI = int(date_tmp(:,1:2))
            month_FMI = int(date_tmp(:,4:5))
            year_FMI = int(date_tmp(:,7:10))
            time_tmp = str( Data.textdata(2:end-1,4))
            hour_FMI = int(time_tmp(:,1:2))
            min_FMI = int(time_tmp(:,4:5))
            date_t_FMI = [date_t_FMI; datenum([year_FMI, month_FMI,day_FMI, hour_FMI, min_FMI,zeros(len(min_FMI),1)])]
            temp_FMI = [temp_FMI; int(str(Data.textdata(2:end-1,6)))]
            rh_FMI = [rh_FMI; int(str(Data.textdata(2:end-1,9)))]
            date_sn_FMI = date_t_FMI
            sn_temp_FMI = temp_FMI
            sn_rh_FMI = rh_FMI
            sn_FMI = [sn_FMI; Data.data(:,2)]
            date_rr_FMI = date_t_FMI
            rr_FMI = [rr_FMI; Data.data(:,1)]        
        end
    
        findx = 0
        indx1 = 0
        jj = 0
        file_name = {}
        # search the pressure data parafile.writefrom J�ms�
        flist_FMI = os.listdir(os.path.join(folder_FMI,'FMI2301_20*_T_TD_RH_P.csv'));
        # Select files that fall within start and end times
        for (findx in range(0,len(flist_FMI))):
            filename = flist_FMI(findx,1).name
            file_year(findx) = int([filename(9:12)])
            file_month(findx) = int([filename(13:14)])
            if file_month(findx) == int(datestr(event_start_day, 'mm')) and file_year(findx) == int(datestr(event_start_day, 'yyyy'))
                indx1      = indx1+1
                file_name(indx1)= {flist_FMI(findx,1).name}
            elif file_month(findx) == int(datestr(event_end_day, 'mm')) and file_year(findx) == int(datestr(event_end_day, 'yyyy'))
                indx1      = indx1+1
                file_name(indx1)= {flist_FMI(findx,1).name}
        for (jj in range(0,indx1):
            fname_FMI = [folder_FMI,str(file_name(jj))]
            Data = np.loadtxt(fname_FMI)
            date_tmp = str( Data.textdata(2:end-1,3))
            day_FMI_J = int(date_tmp(:,1:2))
            month_FMI_J = int(date_tmp(:,4:5))
            year_FMI_J = int(date_tmp(:,7:10))
            time_tmp = str( Data.textdata(2:end-1,4))
            hour_FMI_J = int(time_tmp(:,1:2))
            min_FMI_J = int(time_tmp(:,4:5))
            date_press_FMI_J = [date_press_FMI_J; datenum([year_FMI_J, month_FMI_J,day_FMI_J, hour_FMI_J, min_FMI_J,zeros(len(min_FMI_J),1)])]
            press_FMI_J = [press_FMI_J; int(str(Data.textdata(2:end-1,9)))]    
            '''
    """
    Read data into arrays for events later than 11.1.2017
    """
    """
    #elif
    if event_start_day > datetime(2017,01,11):
        
        # search the temperature and relative humidity, precipitation rate and snow depth files
        #flist_FMI = glob(os.path.join(folder_FMI,'FMI2303_20*_TA_TG_TD_RH_PRON_RI_SND.csv'))
        flist_FMI = glob(os.path.join(folder_FMI,'FMI2303_20*_T_TD_TW_RH.csv'))
        print(flist_FMI)
        jj = 0
        # Select files that fall within start and end times
        file_name, indx1 = select_files(flist_FMI, event_start_day, event_end_day)
    
        for jj in range(0,indx1):
            fname_FMI = os.path.join(folder_FMI,str(file_name[jj]))
            Data = np.loadtxt(fname_FMI, dtype=str, delimiter=',')
            print(Data)
            date_tmp = str(Data[1:len(Data[2])-1,2])
            print(date_tmp[2:3])
            day_FMI = str(date_tmp[2:4])
            month_FMI = str(date_tmp[5:7])
            year_FMI = str(date_tmp[8:12])
            time_tmp = str(Data[1:len(Data[3])-1,3])
            print(time_tmp)
            hour_FMI = str(time_tmp[2:4])
            min_FMI = str(time_tmp[5:7])
            date_t_FMI = np.vstack((date_t_FMI, np.array([year_FMI , month_FMI , day_FMI , hour_FMI , min_FMI])))
            if temp_FMI in dir():
                temp_FMI = np.vstack((temp_FMI, Data[:,2]))
            else:
                temp_FMI = np.empty_like(Data[:,2])
                temp_FMI = np.vstack((temp_FMI, Data[:,2]))
            if rh_FMI in dir():
                rh_FMI = np.vstack((rh_FMI, Data[:,5]))
            else:
                rh_FMI = np.empty_like(Data[:,5])
                rh_FMI = np.vstack((rh_FMI, Data[:,5]))
            date_sn_FMI = date_t_FMI
            sn_temp_FMI = temp_FMI
            sn_rh_FMI = rh_FMI
            if sn_FMI in dir():
                sn_FMI = np.vstack((sn_FMI, Data[:,8]))
            else:
                sn_FMI = np.empty_like(Data[:,8])
                sn_FMI = np.vstack((sn_FMI, Data[:,8]))
            date_rr_FMI = date_t_FMI
            if rr_FMI in dir():
                rr_FMI = np.vstack((rr_FMI, Data[:,1]))
            else:
                rr_FMI = np.empty_like(Data[:,1])
                rr_FMI = np.vstack((rr_FMI, Data[:,1]))
            

        jj = 0
        # search the pressure data from J�ms�
        flist_FMI = glob(os.path.join(folder_FMI,'FMI2301_20*_TA_RH_TD_PA_PA0_WS_WD_WG_VIS_WAWA.csv'))
        # Select files that fall within start and end times
        file_name, indx1 = select_files(flist_FMI, event_start_day, event_end_day)
            
        for jj in range(0,indx1):
            fname_FMI = [folder_FMI,str(file_name(jj))]
            delimiter = ','
            startRow = 2
            formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]'
            fileID = fopen(fname_FMI,'r')
            Data = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n')
            fclose(fileID)
            #Data = importdata(fname_FMI, delimiterIn,headerlinesIn)
            date_tmp = str(Data(3))
            day_FMI_J = int(date_tmp[:,0:2])
            month_FMI_J = int(date_tmp[:,3:5])
            year_FMI_J = int(date_tmp[:,6:10])
            time_tmp = str(Data(4))
            hour_FMI_J = int(time_tmp[:,0:2])
            min_FMI_J = int(time_tmp[:,3:5])
            date_press_FMI_J = np.vstack(date_press_FMI_J, np.row_stack(year_FMI_J, month_FMI_J,day_FMI_J, hour_FMI_J, min_FMI_J,zeros(len(min_FMI_J),1)))
            press_FMI_J = np.vstack(press_FMI_J, int(str(Data(9))))
    """        
    
    """
    Vectorization and return
    """
    sn_FMI = np.array(sn_FMI)
    sn_FMI[sn_FMI <0]=0 
    print(press_FMI_J)
    # Divide observations to n-minute vector
    #print(n/(24.*60.))
    tv_FMI = np.arange(event_start_time, event_end_time, timedelta(minutes=n)) 
    #print(tv_FMI)
    
    press_FMI_J_tv = press_FMI_J
    temp_FMI_tv = np.empty_like(temp_FMI, dtype=float)
    rh_FMI_tv = np.empty_like(rh_FMI, dtype=float)
    sn_FMI_tv = np.empty_like(sn_FMI, dtype=float)
    rr_FMI_tv = np.empty_like(rr_FMI, dtype=float)
    
    for vv in range(0,len(tv_FMI)):
        ww = np.where(date_press_FMI_J > tv_FMI[vv-1] and date_press_FMI_J <= tv_FMI[vv])
        print("Kierros " + str(vv))
        print(np.shape(press_FMI_J_tv))
        if np.all(ww == 0) == False:
            press_FMI_J_tv[0, vv] = np.NaN
        else:
            press_FMI_J_tv[0, vv] = press_FMI_J[ww]
        
        xx = np.where(date_t_FMI > tv_FMI[vv-1] and date_t_FMI <= tv_FMI[vv])
        if np.all(xx == 0) == False:
            temp_FMI_tv[0, vv] = np.NaN
            rh_FMI_tv[0, vv] = np.NaN
        else:
            temp_FMI_tv[0, vv] = temp_FMI[0,xx]
            rh_FMI_tv[0, vv] = rh_FMI[0,xx]
        
        yy = np.where(date_sn_FMI > tv_FMI[vv-1] and date_sn_FMI <= tv_FMI[vv])
        if np.all(yy == 0) == False:
            sn_FMI_tv[0, vv] = np.NaN
        else:
            sn_FMI_tv[0, vv] = sn_FMI(yy[len(yy)])
        
        zz = np.where(date_rr_FMI > tv_FMI[vv-1] and date_rr_FMI <= tv_FMI[vv])
        if np.all(zz == 0) == False:
            rr_FMI_tv[0, vv] = np.NaN
        else:
            rr_FMI_tv[0, vv] = rr_FMI
    print(tv_FMI, temp_FMI_tv, rh_FMI_tv, sn_FMI_tv, rr_FMI_tv, press_FMI_J_tv)
    
    np.savetxt("readFMIstation_presstest.txt", press_FMI_J_tv, fmt="%s")
    np.savetxt("readFMIstation_temptest.txt", temp_FMI_tv, fmt="%s")
    np.savetxt("readFMIstation_rhtest.txt", rh_FMI_tv, fmt="%s")
    np.savetxt("readFMIstation_sntest.txt", sn_FMI_tv, fmt="%s")
    np.savetxt("readFMIstation_rrtest.txt", rr_FMI_tv, fmt="%s")
    return tv_FMI, temp_FMI_tv, rh_FMI_tv, sn_FMI_tv, rr_FMI_tv, press_FMI_J_tv
    
def select_files(flist_FMI, event_start_day, event_end_day):
    findx = 0
    indx1 = 0
    jj = 0
    file_year = []
    file_month = []
    file_name= []
    
    for findx in range(0,len(flist_FMI)):
        filename = flist_FMI[findx]
        #print(os.path.basename(filename)[8:12].strip())
        file_year.append(int(os.path.basename(filename)[8:12].strip()))
        file_month.append(int(os.path.basename(filename)[12:14].strip()))
        if file_month[findx] == event_start_day.month and file_year[findx] == event_start_day.year:
            indx1 = indx1+1
            file_name.append(os.path.basename(flist_FMI[findx]))
        elif file_month[findx] == event_end_day.month and file_year[findx] == event_end_day.year:
            indx1 = indx1+1
            file_name.append(os.path.basename(flist_FMI[findx]))
            
    return file_name, indx1
         
    
Read_FMIstation_GPM(datetime.strptime("201403100000", "%Y%m%d%H%M"),datetime.strptime("201403100500", "%Y%m%d%H%M"), 5) 