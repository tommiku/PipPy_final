# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 10:53:41 2019

@author: kumlin
"""

# Compare the wind data from Hotplate, 3D - Metek inside the fence and ARM 10 m mast
# Gill outside the fence with 60 s averaged wind speed and direction


def WindComparison_GPM(n, flag,event_start_time,event_end_time):

    if (flag.wind_GILL == 1):
        # GILL in the BAECC measurement field
        # FIND FOLDER HERE
        folder_GILL  = 'G:\GPM_overpass_data\GILL\'
     
        DELIMITER = ' '
        HEADERLINES = 2
    
        indx1 = 0
        flist_GILL = os.listdir(os.path.join(folder_GILL,'GILL_HY20*.wnd'))
        # Select files that fall within start and end times
        
        for (findx = xrange(0,length(flist_GILL)):
            filename = flist_GILL(findx,1).name
            file_time(findx) = datenum([str2num([filename(8:11)]), str2num(filename(12:13)), str2num(filename(14:15)),0,0,0])
            if file_time(findx) >= floor(event_start_time) and file_time(findx)<ceil(event_end_time):
                indx1      = indx1+1
                file_name(indx1)= {flist_GILL(findx,1).name}
        date_GILL = []
        dir_GILL = []
        vel_GILL = []
        temp_GILL = []
        velz_GILL = []
    
        # direction correction in deg
        dd_GILL = 13
    
        # Import the file
        for (ii in xrange(0:indx1)):
            Data = np.loadtxt([folder_GILL,str(file_name(ii))], DELIMITER, HEADERLINES)
            # time/h dd/degs ff/(m/s) ca/(m/s) tv/(C) w/(m/s) sample
            # time in h from midnight
            filename = str(file_name(ii))
            time_yyyy = str2num(filename(8:11))
            time_mm = str2num(filename(12:13))
            time_dd = str2num(filename(14:15))
            time_tmp = Data.data(:,1)
            time_h = floor(time_tmp)
            time_min = floor((time_tmp-time_h)*60)
            time_s = floor((time_tmp-time_h-time_min/60)*3600)
            date_GILL = [date_GILL datenum([time_yyyy*ones(size(time_h,1),1), time_mm*ones(size(time_h,1),1),time_dd*ones(size(time_h,1),1), time_h, time_min, time_s])]
            dir_GILL = [dir_GILL Data.data(:,2)]
            vel_GILL = [vel_GILL Data.data(:,3)]
            temp_GILL = [temp_GILL Data.data(:,5)]
            velz_GILL = [velz_GILL Data.data(:,6)]
        # average every n minutes
        # Define time vector
        if (np.isempty(date_GILL) == 0):
            min_start = floor(str2num(datestr(date_GILL(2), 'MM'))/n)
            min_end = floor(str2num(datestr(date_GILL(end), 'MM'))/n)
            time_start = datenum([str2num(datestr(date_GILL(2), 'yyyy')), str2num(datestr(date_GILL(2), 'mm')), str2num(datestr(date_GILL(2), 'dd')),str2num(datestr(date_GILL(2),'HH')),(min_start*n), 0])
            time_end = datenum([str2num(datestr(date_GILL(end), 'yyyy')),str2num(datestr(date_GILL(end), 'mm')),str2num(datestr(date_GILL(end), 'dd')),str2num(datestr(date_GILL(end),'HH')),(min_end*n), 0])
            time_vector_GILL = time_start:n/(24*60):time_end
            mean_vel_GILL = zeros(1,length(time_vector_GILL))
            mode_dir_GILL = zeros(1,length(time_vector_GILL))
            for (aa = 1:size(time_vector_GILL,2)-1):
                indx = find(date_GILL > time_vector_GILL(aa)+30/(60*60*24) & date_GILL <= time_vector_GILL(aa+1)+30/(60*60*24))
                mean_vel_GILL(aa) = nansum(vel_GILL(indx))/n
                mode_dir_GILL(aa) = mode(dir_GILL(indx))
        mean_vel_GILL = mean_vel_GILL(time_vector_GILL >= event_start_time & time_vector_GILL <= event_end_time)
        mode_dir_GILL = mode_dir_GILL(time_vector_GILL >= event_start_time & time_vector_GILL <= event_end_time)
        time_vector_GILL = time_vector_GILL(time_vector_GILL >= event_start_time & time_vector_GILL <= event_end_time)
    else:
        time_vector_GILL = []
        mean_vel_GILL = []
        mode_dir_GILL = []
    
    if (flag.wind_Metek == 1):
        # 3D-Metek
        # calculate the day nr
        # organize the data to vectors and take no data entries away
        folder_3DMetek  = 'G:\GPM_overpass_data\Metek\'
        
        yy_s = datestr(event_start_time,'yyyy')
        mm_s = datestr(event_start_time,'mm')
        day_s = datestr(event_start_time,'dd')
        
        yy_e = datestr(event_end_time,'yyyy')
        mm_e = datestr(event_end_time,'mm')
        day_e = datestr(event_end_time,'dd')
        
        if (str2num(yy_s) == 2014):
            day0 = datenum([2014,1,1,00,00,00])
            deltaday = datenum([2014,1,2,00,00,00])-day0
        elif (str2num(yy_s) == 2015):
            day0 = datenum([2015,1,1,00,00,00])
            deltaday = datenum([2015,1,2,00,00,00])-day0
        day_Ms = datenum([str2num(yy_s),str2num(mm_s),str2num(day_s),00,00,00])
        day_Me = datenum([str2num(yy_e),str2num(mm_e),str2num(day_e),00,00,00])
        day_nr_s = floor((day_Ms-day0)/deltaday)+1
        day_nr_e = floor((day_Me-day0)/deltaday)+1
        if (str2num(yy_s) == 2014):
            if day_nr_s < 10 :
                fname_Metek_s = ['G:\GPM_overpass_data\Metek\V1400',num2str(day_nr_s),'_wd.txt']
            elif day_nr_s < 100 :
                fname_Metek_s = ['G:\GPM_overpass_data\Metek\V140',num2str(day_nr_s),'_wd.txt']
            else:
                fname_Metek_s = ['G:\GPM_overpass_data\Metek\V14',num2str(day_nr_s),'_wd.txt']
        elif (str2num(yy_s) == 2015):
            if day_nr_s < 10:
                fname_Metek_s = ['G:\GPM_overpass_data\Metek\V1500',num2str(day_nr_s),'_wd.txt']
            elif day_nr_s < 100:
                fname_Metek_s = ['G:\GPM_overpass_data\Metek\V150',num2str(day_nr_s),'_wd.txt']
            else:
                fname_Metek_s = ['G:\GPM_overpass_data\Metek\V15',num2str(day_nr_s),'_wd.txt']
        elif (str2num(yy_s) == 2016):
            if day_nr_s < 10:
                fname_Metek_s = ['G:\GPM_overpass_data\Metek\V1600',num2str(day_nr_s),'_wd.txt']
            elif day_nr_s < 100:
                fname_Metek_s = ['G:\GPM_overpass_data\Metek\V160',num2str(day_nr_s),'_wd.txt']
            else:
                fname_Metek_s = ['G:\GPM_overpass_data\Metek\V16',num2str(day_nr_s),'_wd.txt']
        if (str2num(yy_e) == 2014):
            if day_nr_e < 10:
                fname_Metek_e = ['G:\GPM_overpass_data\Metek\V1400',num2str(day_nr_e),'_wd.txt']
            elif day_nr_e < 100:
                fname_Metek_e = ['G:\GPM_overpass_data\Metek\V140',num2str(day_nr_e),'_wd.txt']
            else:
                fname_Metek_e = ['G:\GPM_overpass_data\Metek\V14',num2str(day_nr_e),'_wd.txt']
        elif (str2num(yy_e) == 2015):
            if day_nr_e < 10:
                fname_Metek_e = ['G:\GPM_overpass_data\Metek\V1500',num2str(day_nr_e),'_wd.txt']
            elif day_nr_e < 100:
                fname_Metek_e = ['G:\GPM_overpass_data\Metek\V150',num2str(day_nr_e),'_wd.txt']
            else:
                fname_Metek_e = ['G:\GPM_overpass_data\Metek\V15',num2str(day_nr_e),'_wd.txt']
        elif (str2num(yy_e) == 2016):
            if day_nr_e < 10:
                fname_Metek_e = ['G:\GPM_overpass_data\Metek\V1600',num2str(day_nr_e),'_wd.txt']
            elif day_nr_e < 100:
                fname_Metek_e = ['G:\GPM_overpass_data\Metek\V160',num2str(day_nr_e),'_wd.txt']
            else:
                fname_Metek_e = ['G:\GPM_overpass_data\Metek\V16',num2str(day_nr_e),'_wd.txt']
        # averaged over n minutes. Observations done every 3 seconds    
        [date_Metek_s, temp_Metek_s, vel_Metek_s, dir_Metek_s, vel_z_Metek_s] = Read2DVDWind(fname_Metek_s,n)
        [date_Metek_e, temp_Metek_e, vel_Metek_e, dir_Metek_e, vel_z_Metek_e] = Read2DVDWind(fname_Metek_e,n)
        cc = setdiff(date_Metek_s, date_Metek_e) 
        if np.isempty(cc):
            time_Metek = [date_Metek_s(date_Metek_s>=event_start_time & date_Metek_s<=event_end_time)]
            vel_Metek = [vel_Metek_s(date_Metek_s>=event_start_time & date_Metek_s<=event_end_time)]
            dir_Metek = [dir_Metek_s(date_Metek_s>=event_start_time & date_Metek_s<=event_end_time)]
        else:
            # Select the time of the event
            time_Metek = [date_Metek_s(date_Metek_s>=event_start_time & date_Metek_s<=event_end_time) date_Metek_e(date_Metek_e>=event_start_time & date_Metek_e<=event_end_time)]
            vel_Metek = [vel_Metek_s(date_Metek_s>=event_start_time & date_Metek_s<=event_end_time) vel_Metek_e(date_Metek_e>=event_start_time & date_Metek_e<=event_end_time)]
            dir_Metek = [dir_Metek_s(date_Metek_s>=event_start_time & date_Metek_s<=event_end_time) dir_Metek_e(date_Metek_e>=event_start_time & date_Metek_e<=event_end_time)]
     else:
        time_Metek = []
        vel_Metek = []
        dir_Metek = []

