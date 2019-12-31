# -*- coding: utf-8 -*-
 
# Read the velocity data of PIP or take it from the mat-file

def ReadPIP_vel_GPM(event_start_time,event_end_time):

    yy_s = datestr(event_start_time, 'yyyy')
    mm_s = datestr(event_start_time, 'mm')
    day_s = datestr(event_start_time, 'dd')
    yy_e = datestr(event_end_time, 'yyyy')
    mm_e = datestr(event_end_time, 'mm')
    day_e = datestr(event_end_time, 'dd')
    
    
    
    
    # File format changed 
    #upgrade_day = datenum(2014,11,23,0,0,0)
    upgrade_day = datenum(2014,1,31,0,0,0)
    # Open PIP files
    if event_start_time < upgrade_day:
        D_PIP = 0.125:0.25:26
        D_PIP = [D_PIP 26]
    else:
        D_PIP = 0.1:0.2:26
        D_PIP = [D_PIP 26]
    
    # for velocity
    if event_start_time < upgrade_day:
        folder_PIP = 'M:\GPM_overpass_data\PIP\f_2_2_Velocity_Tables\'
        subfolder1 = [folder_PIP,'004', yy_s,mm_s,day_s]
        flist_1 = os.listdir(fullfile(subfolder1,'*_a_v_2.dat'))
        subfolder2 = [folder_PIP,'004', yy_e,mm_e,day_e]
        flist_2 = os.listdir(fullfile(subfolder2,'*_a_v_2.dat'))
    else:
        folder_PIP = 'M:\GPM_overpass_data\PIP\f_2_2_Velocity_Tables\'
        #subfolder1 = [folder_PIP,'004', yy_s,mm_s,day_s,'_new']
        subfolder1 = [folder_PIP,'004', yy_s,mm_s,day_s]
        flist_1 = os.listdir(fullfile(subfolder1,'*_a_v_2.dat'))
        #subfolder2 = [folder_PIP,'004', yy_e,mm_e,day_e,'_new']
        subfolder2 = [folder_PIP,'004', yy_e,mm_e,day_e]
        flist_2 = os.listdir(fullfile(subfolder2,'*_a_v_2.dat'))
    
    indx1 = 0
    indx2 = 0
    # Select files that fall within start and end times
    for (findx in xrange(0:length(flist_1)):
        filename = flist_1(findx,1).name
        file_time.day1(findx) = datenum([int([filename(4:7)]), int(filename(8:9)), int(filename(10:11)), int(filename(12:13)), int(filename(14:15)), 0])
        if (file_time.day1(findx) >= event_start_time and file_time.day1(findx)<=event_end_time):
            indx1      = indx1+1
            file_name.day1(indx1)= {flist_1(findx,1).name}
        
    
    
    PIPD1 =[]
    PIPV1 =[]
    PIPtime1 =[]
    
    for (ii in xrange(0:indx1)):
        DELIMITER = '\t'
        HEADERLINES = 9
        # Import the file
        Data1 = np.loadtxt([str(subfolder1),'\',str(file_name.day1(ii))], DELIMITER, HEADERLINES)
        ID1 = Data1.data(:,2)
        D1 = Data1.data(:,3)
        VV1 = Data1.data(:,5)
        minutes = Data1.data(:,6)
        filename = str(file_name.day1(ii))
        time1 = datenum([int(filename(4:7))*np.ones(size(minutes)), int(filename(8:9))*np.ones(size(minutes)), int(filename(10:11))*np.ones(size(minutes)), int(filename(12:13))*np.ones(size(minutes)),minutes, zeros(size(minutes))])
        PIP(ii).ID1 = []
        PIP(ii).D1 = []
        PIP(ii).VV1 = []
        PIP(ii).time1 = []
        [C,ia,ic] = unique(ID1) 
        for jj = 1:size(ia)
            PIP(ii).ID1 = [PIP(ii).ID1 C(jj)]
            PIP(ii).D1 = [PIP(ii).D1 mean(D1(ID1==C(jj)))]
            PIP(ii).VV1 = [PIP(ii).VV1 mean(VV1(ID1==C(jj)))]
            PIP(ii).time1 = [PIP(ii).time1 time1(ia(jj))]
        
        PIPD1 =[PIPD1 PIP(ii).D1]
        PIPV1 =[PIPV1 PIP(ii).VV1]
        PIPtime1 =[PIPtime1 PIP(ii).time1]
    
    
    for findx2 = 1:length(flist_2)
        filename = flist_2(findx2,1).name
        file_time.day2(findx2) = datenum([int([filename(4:7)]), int(filename(8:9)), int(filename(10:11)), int(filename(12:13)), int(filename(14:15)), 0])
        if file_time.day2(findx2) >= event_start_time and file_time.day2(findx2)<=event_end_time
            indx2      = indx2+1
            file_name.day2(indx2)= {flist_2(findx2,1).name}
        
    
    
    PIPD2 =[]
    PIPV2 =[]
    PIPtime2 =[]
    for ii = 1:indx2
        DELIMITER = '\t'
        HEADERLINES = 9
        # Import the file
        Data2 = np.loadtxt([str(subfolder2),'\',str(file_name.day2(ii))], DELIMITER, HEADERLINES)
        ID2 = Data2.data(:,2)
        D2 = Data2.data(:,3)
        VV2 = Data2.data(:,5)
        minutes = Data2.data(:,6)
        filename = str(file_name.day2(ii))
        time2 = datenum([int(filename(4:7))*np.ones(size(minutes)), int(filename(8:9))*np.ones(size(minutes)), int(filename(10:11))*np.ones(size(minutes)), int(filename(12:13))*np.ones(size(minutes)),minutes, zeros(size(minutes))])
        PIP(ii).ID2 = []
        PIP(ii).D2 = []
        PIP(ii).VV2 = []
        PIP(ii).time2 = []
        [C,ia,ic] = unique(ID2) 
        for jj = 1:size(ia)
            PIP(ii).ID2 = [PIP(ii).ID2 C(jj)]
            PIP(ii).D2 = [PIP(ii).D2 mean(D2(ID2==C(jj)))]
            PIP(ii).VV2 = [PIP(ii).VV2 mean(VV2(ID2==C(jj)))]
            PIP(ii).time2 = [PIP(ii).time2 time2(ia(jj))]
         
         PIPD2 =[PIPD2 PIP(ii).D2]
         PIPV2 =[PIPV2 PIP(ii).VV2]
         PIPtime2 = [PIPtime2 PIP(ii).time2]
    
    
    cc = setdiff(PIPD1, PIPD2) 
    if isempty(cc)
        PIPD = [PIPD1]
        PIPV = [PIPV1]
        PIPtime = [PIPtime1]
    else
        PIPD = [PIPD1PIPD2]
        PIPV = [PIPV1PIPV2]
        PIPtime = [PIPtime1PIPtime2]
    
        
    PIPD(PIPV == -99) = []
    PIPtime(PIPV == -99) = []
    PIPV(PIPV == -99) = []
    PIPV(PIPD == -99) = []
    PIPtime(PIPD == -99) = []
    PIPD(PIPD == -99) = []

