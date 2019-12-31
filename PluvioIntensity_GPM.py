# -*- coding: utf-8 -*-
import glob
#VEDEN TAI LUMEN PAINO
def PluvioIntensity_GPM(n, yy_s, mm_s, day_s, yy_e, mm_e, day_e, event_start_time, event_end_time):
    # Accumulation comparison
    # Plot the n accumulated liquid precipitation data of Pluvios,
    
    # Open Pluvios
    # TODO: Read from config
    folder_Pluvio_200 = '/home/kumlin/Koodit/SnowRetrievals/Data_malli_10062019/Pluvio200'
    folder_Pluvio_400 = '/home/kumlin/Koodit/SnowRetrievals/Data_malli_10062019/Pluvio400'
    
    flist_200_1 = glob.glob(os.path.join(folder_Pluvio_200,['pluvio200_02_',yy_s,mm_s,day_s,'*.txt']))
    flist_200_2 = glob.glob(os.path.join(folder_Pluvio_200,['pluvio200_02_',yy_e,mm_e,day_e,'*.txt']))
    flist_400_1 = glob.glob(os.path.join(folder_Pluvio_400,['pluvio400_01_',yy_s,mm_s,day_s,'*.txt']))
    flist_400_2 = glob.glob(os.path.join(folder_Pluvio_400,['pluvio400_01_',yy_e,mm_e,day_e,'*.txt']))
    
    indx_200_1 = 0
    indx_200_2 = 0
    indx_400_1 = 0
    indx_400_2 = 0
    
    event_start_day = floor(event_start_time)
    event_end_day = ceil(event_end_time)
    # Select files that fall within start and  times
    for findx = 0:len(flist_200_1):
        filename = flist_200_1(findx,1).name
        file_time.PL200_1(findx) = datenum([str2num([filename(14:17)]), str2num(filename(18:19)), str2num(filename(20:21)),str2num(filename(22:23)), 59, 59])
        if file_time.PL200_1(findx) >= event_start_day && file_time.PL200_1(findx)<=event_end_day:
            indx_200_1      = indx_200_1+1
            file_name.PL200_1(indx_200_1)= {flist_200_1(findx,1).name}
        
    
    
    
    # Select files that fall within start and  times
    for findx = 0:len(flist_200_2):
        filename = flist_200_2(findx,1).name
        file_time.PL200_2(findx) = datenum([str2num([filename(14:17)]), str2num(filename(18:19)), str2num(filename(20:21)),str2num(filename(22:23)), 0, 0])
        if file_time.PL200_2(findx) >= event_start_day && file_time.PL200_2(findx)<=event_end_day:
            indx_200_2      = indx_200_2+1
            file_name.PL200_2(indx_200_2)= {flist_200_2(findx,1).name}
        
    
    
    # Select files that fall within start and  times
    for findx = 0:len(flist_400_1):
        filename = flist_400_1(findx,1).name
        file_time.PL400_1(findx) = datenum([str2num([filename(14:17)]), str2num(filename(18:19)), str2num(filename(20:21)),str2num(filename(22:23)), 59, 59])
        if file_time.PL400_1(findx) >= event_start_day && file_time.PL400_1(findx)<=event_end_day:
            indx_400_1      = indx_400_1+1
            file_name.PL400_1(indx_400_1)= {flist_400_1(findx,1).name}
        
    
    
    
    # Select files that fall within start and  times
    for findx = 0:len(flist_400_2):
        filename = flist_400_2(findx,1).name
        file_time.PL400_2(findx) = datenum([str2num([filename(14:17)]), str2num(filename(18:19)), str2num(filename(20:21)), str2num(filename(22:23)), 0, 0])
        if file_time.PL400_2(findx) >= event_start_day && file_time.PL400_2(findx)<=event_end_day:
            indx_400_2      = indx_400_2+1
            file_name.PL400_2(indx_400_2)= {flist_400_2(findx,1).name}
        
    
    
    # Data is restored every minute and collected to 1 hour .txt file
    # RT=Real-Time (within a 1 minute of the precipitation event occurring) 
    # NRT=Non-Real-Time (5min after the precipitation event occurring)'
    # Each row, delimiter +
    # timestamp yyyymmddhhmmss
    # 2 Intensity RT  [mm/h or mm/min] Moving precipitation growth over the last minute before the sample interval
    # 3 Accumulated RT/NRT [mm] sum of the correct volumes of precipitation over the sample interval
    # 4 Accumulated NRT [mm] sum of the correct amounts of precipitation over the sample interval with a fixed output delay of 5 minutes
    # 5 Accumulated total NRT [mm] sum of the correct amounts of precipitation since the last device start with a fixed output delay of 5 minutes.
    # 6 Bucket RT [mm] current unfiltered bucket content
    # 7 Bucket NRT [mm] current measured filtered bucket content
    # 8 Temperature load cell [degC] internal temperature of the load cell used for compensating temperature changes
    # 9 Heating status
    # 10 Status
    # 11 Temperature electronics unit
    # 12 Supply Voltage
    # 13 Temperature orfice ring rim
    #20140110130001+0000.00+0000.00+0000.00+0059.69+0331.38+0331.36+00.6+128+000+00.6+24.4+00.0
    
    DELIMITER = ';'
    HEADERLINES = 0
    
    PL200_1_rr =[]
    PL200_2_rr =[]
    PL400_1_rr =[]
    PL400_2_rr =[]
    PL200_1_acc =[]
    PL200_2_acc =[]
    PL400_1_acc =[]
    PL400_2_acc =[]
    PL200_1_t =[]
    PL200_2_t =[]
    PL400_1_t =[]
    PL400_2_t =[]
    
    for ii = 0:indx_200_1:
        # Import the file
        rawData = importdata([str(folder_Pluvio_200),str(file_name.PL200_1(ii))], DELIMITER, HEADERLINES)
        time_vector = num2str(rawData(:,1))
        time = datenum( [str2num(time_vector(:,1:4)),
        str2num(time_vector(:,5:6)),str2num(time_vector(:,7:8)),
        str2num(time_vector(:,9:10)),str2num(time_vector(:,11:12)),
        str2num(time_vector(:,13:14))])
        dataPluvio = rawData(:,2)
        PL200_1_rr =[PL200_1_rr dataPluvio]
        dataPluvio = rawData(:,7)
        PL200_1_acc = [PL200_1_acc dataPluvio]
        PL200_1_t =[PL200_1_t time]
    
    
    for ii = 0:indx_200_2:
        # Import the file
        rawData = importdata([str(folder_Pluvio_200),str(file_name.PL200_2(ii))], DELIMITER, HEADERLINES)
        if np.isempty(rawData) ==0:
            time_vector = num2str(rawData(:,1))
            time = datenum( [str2num(time_vector(:,1:4)),
            str2num(time_vector(:,5:6)),str2num(time_vector(:,7:8)),
            str2num(time_vector(:,9:10)),str2num(time_vector(:,11:12)),
            str2num(time_vector(:,13:14))])
            dataPluvio = rawData(:,2)
            PL200_2_rr =[PL200_2_rr dataPluvio]
            dataPluvio = rawData(:,7)
            PL200_2_acc =[PL200_2_acc dataPluvio]
            PL200_2_t =[PL200_2_t time]
        
    
    
    for ii = 0:indx_400_1:
        # Import the file
        rawData = importdata([str(folder_Pluvio_400),str(file_name.PL400_1(ii))], DELIMITER, HEADERLINES)
        if np.isempty(rawData) == 0:
            time_vector = num2str(rawData(:,1))
            time = datenum( [str2num(time_vector(:,1:4)),
            str2num(time_vector(:,5:6)),str2num(time_vector(:,7:8)),
            str2num(time_vector(:,9:10)),str2num(time_vector(:,11:12)),
            str2num(time_vector(:,13:14))])
            dataPluvio = rawData(:,2)
            PL400_1_rr =[PL400_1_rr dataPluvio]
            dataPluvio = rawData(:,7)
            PL400_1_acc =[PL400_1_acc dataPluvio]
            PL400_1_t =[PL400_1_t time]
        
    
    
    for ii = 0:indx_400_2:
        # Import the file
        rawData = importdata([str(folder_Pluvio_400),str(file_name.PL400_2(ii))], DELIMITER, HEADERLINES)
        if np.isempty(rawData) == 0:
            time_vector = num2str(rawData(:,1))
            time = datenum( [str2num(time_vector(:,1:4)),
            str2num(time_vector(:,5:6)),str2num(time_vector(:,7:8)),
            str2num(time_vector(:,9:10)),str2num(time_vector(:,11:12)),
            str2num(time_vector(:,13:14))])
            dataPluvio = rawData(:,2)
            PL400_2_rr =[PL400_2_rr dataPluvio]
            dataPluvio = rawData(:,7)
            PL400_2_acc =[PL400_2_acc dataPluvio]
            PL400_2_t =[PL400_2_t time]
        
    
    
    cc = setdiff(PL200_1_t, PL200_2_t) 
    if np.isempty(cc):
        
        PL200_rr = [PL200_1_rr]
        PL400_rr = [PL400_1_rr]
        PL200_acc = [PL200_1_acc]
        PL400_acc = [PL400_1_acc] 
        PL200_t = [PL200_1_t]
        PL400_t = [PL400_1_t]
    else:
        PL200_rr = [PL200_1_rrPL200_2_rr]
        PL400_rr = [PL400_1_rrPL400_2_rr]
        PL200_acc = [PL200_1_accPL200_2_acc]
        PL400_acc = [PL400_1_accPL400_2_acc]
        PL200_t = [PL200_1_tPL200_2_t]
        PL400_t = [PL400_1_tPL400_2_t]
    
    
    
    time_start = datenum([str2num(datestr(PL200_t(1), 'yyyy')), str2num(datestr(PL200_t(1), 'mm')), str2num(datestr(PL200_t(1), 'dd')),0,0, 0])
    time_ = datenum([str2num(datestr(PL200_t(), 'yyyy')),str2num(datestr(PL200_t(), 'mm')),str2num(datestr(PL200_t(), 'dd'))+1,0,0, 0])
    N = floor((time_ - time_start)*(24*60/n))
    time_vectorPL200 = linspace(time_start, time_, N+1)
    
    for aa = 0:np.shape(time_vectorPL200,1)-1:
        indx = find(PL200_t > time_vectorPL200(aa)+30/(60*60*24) & PL200_t <= time_vectorPL200(aa+1)+30/(60*60*24))
        if np.isempty(indx):
            max_PL200_acc(aa) = 0
            mean_PL200_rr(aa) = 0
        else:
            max_PL200_acc(aa) = max(PL200_acc(indx))
            mean_PL200_rr(aa) = (60/len(indx)*(max(PL200_acc(indx))-min(PL200_acc(indx))))
        
    
    
    time_vectorPL200 = time_vectorPL200(2:)
    tv_PL200 = time_vectorPL200(time_vectorPL200<=event_end_time & time_vectorPL200>=event_start_time)
    PL200_acc = max_PL200_acc(time_vectorPL200<=event_end_time & time_vectorPL200>=event_start_time)
    PL200_rr = mean_PL200_rr(time_vectorPL200<=event_end_time & time_vectorPL200>=event_start_time)
    
    time_start = datenum([str2num(datestr(PL400_t(1), 'yyyy')), str2num(datestr(PL400_t(1), 'mm')), str2num(datestr(PL400_t(1), 'dd')),0,0, 0])
    time_ = datenum([str2num(datestr(PL400_t(), 'yyyy')),str2num(datestr(PL400_t(), 'mm')),str2num(datestr(PL400_t(), 'dd'))+1,0,0, 0])
    N = floor((time_ - time_start)*(24*60/n))
    time_vectorPL400 = linspace(time_start, time_, N+1)
    
    for aa = 0:np.shape(time_vectorPL400,1)-1:
        indx = find(PL400_t > time_vectorPL400(aa)+30/(60*60*24) & PL400_t <= time_vectorPL400(aa+1)+30/(60*60*24))
        if np.isempty(indx):
            max_PL400_acc(aa) = 0
            mean_PL400_rr(aa) = 0
        else:
            max_PL400_acc(aa) = max(PL400_acc(indx))
            mean_PL400_rr(aa) = (60/len(indx)*(max(PL400_acc(indx))-min(PL400_acc(indx))))
        
    
    time_vectorPL400 = time_vectorPL400(2:)
    
    # Select the data from the event
    tv_PL400 = time_vectorPL400(time_vectorPL400<=event_end_time & time_vectorPL400>=event_start_time)
    PL400_acc = max_PL400_acc(time_vectorPL400<=event_end_time & time_vectorPL400>=event_start_time)
    PL400_rr = mean_PL400_rr(time_vectorPL400<=event_end_time & time_vectorPL400>=event_start_time)
    
    return tv_PL200, PL200_acc, PL200_rr, tv_PL400, PL400_acc, PL400_rr
    
