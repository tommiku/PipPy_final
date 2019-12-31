# -*- coding: utf-8 -*-

# This subfunction opens the PIP particle tables and defined the
# strateristics of indivudual particles.

def ReadPIP_par_GPM(event_start_time, event_end_time):

    yy_s = datestr(event_start_time, 'yyyy')
    mm_s = datestr(event_start_time, 'mm')
    day_s = datestr(event_start_time, 'dd')
    yy_e = datestr(event_end_time, 'yyyy')
    mm_e = datestr(event_end_time, 'mm')
    day_e = datestr(event_end_time, 'dd')
    
    
    # the data is derived from the raw data tables
    # File format changed 
    #upgrade_day = datenum(2014,11,23,0,0,0)
    upgrade_day = datenum(2014,1,31,0,0,0)
    if event_start_time < upgrade_day
        D_PIP = 0.125:0.25:26
        D_PIP = [D_PIP 26]
    else
        D_PIP = 0.1:0.2:26
        D_PIP = [D_PIP 26]
    
    # Open PIP files
    # select particle table
    if event_start_time < upgrade_day
        folder_PIP_p = 'M:\GPM_overpass_data\PIP\f_1_2_Particle_Tables_ascii\'
        subfolder_p1 = [folder_PIP_p,'004', yy_s,mm_s,day_s]
        flist_p1 = dir(fullfile(subfolder_p1,'*_a_p.dat'))
        subfolder_p2 = [folder_PIP_p,'004', yy_e,mm_e,day_e]
        flist_p2 = dir(fullfile(subfolder_p2,'*_a_p.dat'))
        folder_PIP_v = 'M:\GPM_overpass_data\PIP\f_2_2_Velocity_Tables\'
        subfolder_v1 = [folder_PIP_v,'004', yy_s,mm_s,day_s]
        flist_v1 = dir(fullfile(subfolder_v1,'*_a_v_2.dat'))
        subfolder_v2 = [folder_PIP_v,'004', yy_e,mm_e,day_e]
        flist_v2 = dir(fullfile(subfolder_v2,'*_a_v_2.dat'))
        indx_v1 = 0
        indx_v2 = 0
    else
        folder_PIP_p = 'M:\GPM_overpass_data\PIP\f_1_2_Particle_Tables_ascii\'
        #subfolder_p1 = [folder_PIP_p,'004', yy_s,mm_s,day_s, '_new']
        subfolder_p1 = [folder_PIP_p,'004', yy_s,mm_s,day_s]
        flist_p1 = dir(fullfile(subfolder_p1,'*_a_p.dat'))
        #subfolder_p2 = [folder_PIP_p,'004', yy_e,mm_e,day_e, '_new']
        subfolder_p2 = [folder_PIP_p,'004', yy_e,mm_e,day_e]
        flist_p2 = dir(fullfile(subfolder_p2,'*_a_p.dat'))
        folder_PIP_v = 'M:\GPM_overpass_data\PIP\f_2_2_Velocity_Tables\'
        #subfolder_v1 = [folder_PIP_v,'004', yy_s,mm_s,day_s,'_new']
        subfolder_v1 = [folder_PIP_v,'004', yy_s,mm_s,day_s]
        flist_v1 = dir(fullfile(subfolder_v1,'*_a_v_2.dat'))
        #subfolder_v2 = [folder_PIP_v,'004', yy_e,mm_e,day_e,'_new']
        subfolder_v2 = [folder_PIP_v,'004', yy_e,mm_e,day_e]
        flist_v2 = dir(fullfile(subfolder_v2,'*_a_v_2.dat'))
        indx_v1 = 0
        indx_v2 = 0
    
    
    # Check, if the start time and end time are on consecutive days and add
    # a flag for calculations
    event2d = 0
    if int(day_s) ~= int(day_e)
        event2d = 1
    
    
    # Select files that fall within start and end times
    for findx = 1:length(flist_v1)
        filename = flist_v1(findx,1).name
        file_time.day1(findx) = datenum([int([filename(4:7)]), int(filename(8:9)), int(filename(10:11)), int(filename(12:13)), int(filename(14:15)), 0])
        if file_time.day1(findx) >= event_start_time && file_time.day1(findx)<=event_end_time
            indx_v1 = indx_v1+1
            file_name.day1(indx_v1)= {flist_v1(findx,1).name}
        
    
    
    PIPD1 =[]
    PIPV1 =[]
    PIPtime1 =[]
    PIPEmaj1 =[]
    PIPEmin1 =[]
    PIPEmaj1_calc =[]
    PIPEmaj1max_calc =[]
    PIPEmin1_calc =[]
    PIPAR1 =[]
    PIPlonX1 =[]
    PIPLen1 =[]
    PIPHig1 =[]
    PIPDia1 =[]
    PIPOR1 =[]
    for ii = 1:indx_v1
        # Open first the velocity table
        DELIMITER = '\t'
        HEADERLINES = 9
        # Import the file
        Data1 = importdata([str(subfolder_v1),'\',str(file_name.day1(ii))], DELIMITER, HEADERLINES)
        REC1 = Data1.data(:,1)
        ID1 = Data1.data(:,2)
        D1 = Data1.data(:,3)
        VV1 = Data1.data(:,5)
        minutes = Data1.data(:,6)
        filename = str(file_name.day1(ii))
        time1 = datenum([int(filename(4:7))*ones(size(minutes)), int(filename(8:9))*ones(size(minutes)), int(filename(10:11))*ones(size(minutes)), int(filename(12:13))*ones(size(minutes)),minutes, zeros(size(minutes))])
        PIP(ii).ID1 = []
        RT(ii).nr1 = []
        RT(ii).Emaj1 = []
        RT(ii).Emin1 = []
        RT(ii).Emaj1_calc = []
        RT(ii).Emin1_calc = []
        RT(ii).lonX1 = []
        RT(ii).Lep1 =[]
        RT(ii).Rp1 =[]
        RT(ii).Up1 =[]
        RT(ii).Lop1 =[]
        RT(ii).AR1 = []
        RT(ii).Dia1 = []
        RT(ii).or1 = []
        PIP(ii).Emaj1 = []
        PIP(ii).Emin1 = []
        PIP(ii).Emaj1_calc = []
        PIP(ii).Emaj1max_calc = []
        PIP(ii).Emin1_calc = []
        PIP(ii).AR1 = []
        PIP(ii).Dia1 = []
        PIP(ii).lonX1 = []
        PIP(ii).Len1 =[]
        PIP(ii).Hig1 =[]
        PIP(ii).or1 = []
        PIP(ii).D1 = []
        PIP(ii).VV1 = []
        PIP(ii).time1 = []
        # combine information from the minutes to the ID1
        ID1_min = ID1 + minutes/100 
        ID1_min = ID1_min(ID1_min>0)
        REC1 = REC1(ID1_min>0)
        [C,ia,ic] = unique(ID1_min)
        # look the particle values from the particle table    
        DELIMITER = '\t'
        HEADERLINES = 10
        xfile = str(file_name.day1(ii))
        # For 18 March 2014 the particle table file is split because of matlab
        # memory problems
        if strcmp(xfile, '0042014032018000_a_v_2.dat')
            xxfile = '0042014032019000_a_v_2.dat'
            # find particle tables that are between the the two velocity tables
            indx_p1 = 0
            for pindx = 1:length(flist_p1)
                filename = flist_p1(pindx,1).name
                particle_table.time(pindx) = datenum([int([filename(4:7)]), int(filename(8:9)), int(filename(10:11)), int(filename(12:13)), int(filename(14:15)), 0])
                if particle_table.time(pindx) >= datenum([int(xfile(4:7)), int(xfile(8:9)), int(xfile(10:11)), int(xfile(12:13)), 0,0])&& particle_table.time(pindx)<datenum([int(xxfile(4:7)), int(xxfile(8:9)), int(xxfile(10:11)), int(xxfile(12:13)), 0,0])
                    indx_p1 = indx_p1+1
                    particle_table.name(indx_p1)= {flist_p1(pindx,1).name}
                
                
            dataset1 = []
            for pp = 1:indx_p1
                xpfile = str(particle_table.name(pp)) 
                Data1 = importdata([str(subfolder_p1),'\',xpfile(1:18),'_p.dat'], DELIMITER, HEADERLINES)
                dataset1 = [dataset1 Data1.data]
            
            dataset1(dataset1(:,6) == -1,:) = []
        else
            if event_start_time < upgrade_day
                Data1 = importdata([str(subfolder_p1),'\',xfile(1:18),'_p.dat'], DELIMITER, HEADERLINES)
            else
                #Data1 = importdata([str(subfolder_p1),'\',xfile(1:18),'_p_60.dat'], DELIMITER, HEADERLINES)
                Data1 = importdata([str(subfolder_p1),'\',xfile(1:18),'_p.dat'], DELIMITER, HEADERLINES)    
            
        
        for jj = 1:size(ia)
            PIP(ii).ID1 = [PIP(ii).ID1 floor(C(jj))]
            RT(ii,jj).nr1 = [REC1((ID1_min==C(jj)))]
            for kk = 1:size(RT(ii,jj).nr1,1)
                if strcmp(xfile, '0042014032018000_a_v_2.dat')
                    test = find(RT(ii,jj).nr1(kk) == dataset1(:,1))
                else
                    test = find(RT(ii,jj).nr1(kk) == Data1.data(:,1))
                
                if  isempty(test) ~= 1
                    if strcmp(xfile, '0042014032018000_a_v_2.dat')
                        RT(ii,jj).Emaj1(kk) = [dataset1(test,18)]
                        RT(ii,jj).Emin1(kk) = [dataset1(test,19)]
                        RT(ii,jj).AR1(kk) = [dataset1(test,20)]
                        RT(ii,jj).lonX1(kk) = [dataset1(test,22)]
                        RT(ii,jj).Dia1(kk) = [dataset1(test,27)]
                        RT(ii,jj).or1(kk) = [dataset1(test,23)]
                        RT(ii,jj).Lep1(kk) = [dataset1(test,28)]
                        RT(ii,jj).Rp1(kk) = [dataset1(test,29)]
                        RT(ii,jj).Up1(kk) = [dataset1(test,30)]
                        RT(ii,jj).Lop1(kk) = [dataset1(test,31)]
                    else
                        RT(ii,jj).Emaj1(kk) = [Data1.data(test,18)]
                        RT(ii,jj).Emin1(kk) = [Data1.data(test,19)]
                        RT(ii,jj).AR1(kk) = [Data1.data(test,20)]
                        RT(ii,jj).lonX1(kk) = [Data1.data(test,22)]
                        RT(ii,jj).Dia1(kk) = [Data1.data(test,27)]
                        RT(ii,jj).or1(kk) = [Data1.data(test,23)]
                        RT(ii,jj).Lep1(kk) = [Data1.data(test,28)]
                        RT(ii,jj).Rp1(kk) = [Data1.data(test,29)]
                        RT(ii,jj).Up1(kk) = [Data1.data(test,30)]
                        RT(ii,jj).Lop1(kk) = [Data1.data(test,31)]
                    
    #                 if ii == 1 && jj == 900 && kk == 2
    #                     disp('Stop here')
    #                 end
                  
                        
                    # Derive the ellipse properties from the bounding box
                    # define bounding box dimensions
                    Wmax = 0.1*nonzeros(RT(ii,jj).Rp1(kk)-RT(ii,jj).Lep1(kk))
                    Hmax = 0.1*nonzeros(RT(ii,jj).Lop1(kk)-RT(ii,jj).Up1(kk))
                    # Exception in data for 24 Dec 2014
    #                 if ii == 5 && jj == 3014 && kk == 2
    #                    RT(ii,jj).Lop1(kk) = 214
    #                    RT(ii,jj).Up1(kk) = 208
    #                    Hmax =0.1*6
    #                 end
                    Dim = [Wmax/2 Hmax/2]
                    Orient = RT(ii,jj).or1(kk)*pi/180
                    ini_Emaj = RT(ii,jj).Emaj1(kk)/2 
                    ini_Emin = RT(ii,jj).Emin1(kk)/2 
                    options = optimset('MaxFunEvals',1000)
                    disp(['ii = ', num2str(ii), 'jj = ', num2str(jj), 'kk = ', num2str(kk)])
                    for hh = 1:length(Orient)
                        if Orient(hh) <= pi/2
                            x0 = [sqrt(ini_Emaj(hh)) sqrt(ini_Emin(hh)) atan((ini_Emin(hh)*cos(Orient(hh)))/(ini_Emaj(hh)*sin(Orient(hh)))) atan((-ini_Emaj(hh)*sin(Orient(hh)))/(ini_Emaj(hh)*cos(Orient(hh))))]
                        else
                            Orient(hh) = pi-Orient(hh)
                            x0 =  [sqrt(ini_Emaj(hh)) sqrt(ini_Emin(hh)) atan((ini_Emin(hh)*cos(Orient(hh)))/(ini_Emaj(hh)*sin(Orient(hh)))) atan((-ini_Emaj(hh)*sin(Orient(hh)))/(ini_Emaj(hh)*cos(Orient(hh))))]
                        
                        if Orient(hh) == 0 || Orient(hh) == pi/2 || Orient(hh) == pi
                            a(hh) = Dim(hh,1)
                            b(hh) = Dim(hh,2)
                        else
                            
                            [x, fval, failed] = fsolve(@(x)ellipsefun(x,Dim(hh,:),Orient(hh)), x0, options)
                            if failed == 0
                                a(hh) = 0
                                b(hh) = 0
                            else
                                a(hh) = x(1)^2
                                b(hh) = x(2)^2
                            
                        
                    
                    RT(ii,jj).Emaj1_calc(kk) = max([a b])
                    RT(ii,jj).Emin1_calc(kk) = min([a b])
                
            
            PIP(ii).Emaj1 = [PIP(ii).Emaj1 mean(nonzeros(RT(ii,jj).Emaj1))]
            PIP(ii).Emin1 = [PIP(ii).Emin1 mean(nonzeros(RT(ii,jj).Emin1))]
            PIP(ii).AR1 = [PIP(ii).AR1 mean(nonzeros(RT(ii,jj).AR1))]
            PIP(ii).lonX1 = [PIP(ii).lonX1 mean(nonzeros(RT(ii,jj).lonX1))]
            PIP(ii).Dia1 = [PIP(ii).Dia1 mean(nonzeros(RT(ii,jj).Dia1))]
            PIP(ii).Len1 = [PIP(ii).Len1 mean(nonzeros(RT(ii,jj).Rp1-RT(ii,jj).Lep1))]
            PIP(ii).Hig1 = [PIP(ii).Hig1 mean(nonzeros(RT(ii,jj).Lop1-RT(ii,jj).Up1))]
            PIP(ii).Emaj1_calc = [PIP(ii).Emaj1_calc mean(nonzeros(RT(ii,jj).Emaj1_calc))]            
            if isempty(find(RT(ii,jj).Emaj1_calc))==1
                PIP(ii).Emaj1max_calc = [PIP(ii).Emaj1max_calc NaN]
            else
                PIP(ii).Emaj1max_calc = [PIP(ii).Emaj1max_calc max(nonzeros(RT(ii,jj).Emaj1_calc))]           
            
            PIP(ii).Emin1_calc = [PIP(ii).Emin1_calc mean(nonzeros(RT(ii,jj).Emin1_calc))]
            # check the quality of the iteration
            PIP(ii).Emaj1_calc(2*PIP(ii).Emaj1_calc > sqrt((PIP(ii).Len1*0.1).^2+ (PIP(ii).Hig1*0.1).^2)) = NaN 
            PIP(ii).Emin1_calc(2*PIP(ii).Emaj1_calc > sqrt((PIP(ii).Len1*0.1).^2+ (PIP(ii).Hig1*0.1).^2)) = NaN
            PIP(ii).Emaj1max_calc(2*PIP(ii).Emaj1_calc > sqrt((PIP(ii).Len1*0.1).^2+ (PIP(ii).Hig1*0.1).^2)) = NaN
            PIP(ii).Emaj1_calc(2*PIP(ii).Emin1_calc > sqrt((PIP(ii).Len1*0.1).^2+ (PIP(ii).Hig1*0.1).^2)) = NaN
            PIP(ii).Emin1_calc(2*PIP(ii).Emin1_calc > sqrt((PIP(ii).Len1*0.1).^2+ (PIP(ii).Hig1*0.1).^2)) = NaN
            PIP(ii).Emaj1max_calc(2*PIP(ii).Emin1_calc > sqrt((PIP(ii).Len1*0.1).^2+ (PIP(ii).Hig1*0.1).^2)) = NaN
            PIP(ii).Emaj1_calc(2*PIP(ii).Emaj1max_calc > sqrt((PIP(ii).Len1*0.1).^2+ (PIP(ii).Hig1*0.1).^2)) = NaN
            PIP(ii).Emin1_calc(2*PIP(ii).Emaj1max_calc > sqrt((PIP(ii).Len1*0.1).^2+ (PIP(ii).Hig1*0.1).^2)) = NaN
            PIP(ii).Emaj1max_calc(2*PIP(ii).Emaj1max_calc > sqrt((PIP(ii).Len1*0.1).^2+ (PIP(ii).Hig1*0.1).^2)) = NaN
            PIP(ii).Emaj1_calc(PIP(ii).Emaj1_calc < 0) = NaN 
            PIP(ii).Emin1_calc(PIP(ii).Emaj1_calc < 0) = NaN 
            PIP(ii).Emaj1max_calc(PIP(ii).Emaj1_calc < 0) = NaN 
            PIP(ii).Emaj1_calc(PIP(ii).Emin1_calc < 0) = NaN 
            PIP(ii).Emin1_calc(PIP(ii).Emin1_calc < 0) = NaN 
            PIP(ii).Emaj1max_calc(PIP(ii).Emin1_calc < 0) = NaN 
            PIP(ii).Emaj1_calc(PIP(ii).Emaj1max_calc < 0) = NaN 
            PIP(ii).Emin1_calc(PIP(ii).Emaj1_calc < 0) = NaN 
            PIP(ii).Emaj1max_calc(PIP(ii).Emaj1max_calc < 0) = NaN 
            # check the variance
            or1_var = var(nonzeros(RT(ii,jj).or1))
            or1_var(or1_var == 0) = NaN
            PIP(ii).or1 = [PIP(ii).or1 sqrt(or1_var)]
            PIP(ii).D1 = [PIP(ii).D1 mean(D1(ID1_min==C(jj)))]
            PIP(ii).VV1 = [PIP(ii).VV1 mean(VV1(ID1_min==C(jj)))]
            PIP(ii).time1 = [PIP(ii).time1 time1(ia(jj))]
        
        PIPD1 =[PIPD1 PIP(ii).D1]
        PIPV1 =[PIPV1 PIP(ii).VV1]
        PIPtime1 =[PIPtime1 PIP(ii).time1]
        PIPEmaj1 =[PIPEmaj1 PIP(ii).Emaj1']
        PIPEmin1 =[PIPEmin1 PIP(ii).Emin1']
        PIPAR1 =[PIPAR1 PIP(ii).AR1']
        PIPlonX1 =[PIPlonX1 PIP(ii).lonX1']
        PIPLen1 =[PIPLen1 PIP(ii).Len1']
        PIPHig1 =[PIPHig1 PIP(ii).Hig1']
        PIPDia1 =[PIPDia1 PIP(ii).Dia1']
        PIPOR1 =[PIPOR1 PIP(ii).or1']
        PIPEmaj1_calc =[PIPEmaj1_calc PIP(ii).Emaj1_calc']
        PIPEmin1_calc =[PIPEmin1_calc PIP(ii).Emin1_calc']
        PIPEmaj1max_calc =[PIPEmaj1max_calc PIP(ii).Emaj1max_calc']
    
    
    clear dataset1
    if event2d == 1
        # Select files that fall within start and end times
        for findx = 1:length(flist_v2)
            filename = flist_v2(findx,1).name
            file_time.day2(findx) = datenum([int([filename(4:7)]), int(filename(8:9)), int(filename(10:11)), int(filename(12:13)), int(filename(14:15)), 0])
            if file_time.day2(findx) >= event_start_time && file_time.day2(findx)<=event_end_time
                indx_v2      = indx_v2+1
                file_name.day2(indx_v2)= {flist_v2(findx,1).name}
            
        
    
        PIPD2 =[]
        PIPV2 =[]
        PIPtime2 =[]
        PIPEmaj2 =[]
        PIPEmin2 =[]
        PIPEmaj2_calc =[]
        PIPEmaj2max_calc =[]
        PIPEmin2_calc =[]
        PIPAR2 =[]
        PIPlonX2 =[]
        PIPLen2 =[]
        PIPHig2 =[]
        PIPDia2 =[]
        PIPOR2 =[]
        for ii = 1:indx_v2
            # Open first the velocity table
            DELIMITER = '\t'
            HEADERLINES = 9
            # Import the file
            Data2 = importdata([str(subfolder_v2),'\',str(file_name.day2(ii))], DELIMITER, HEADERLINES)
            REC2 = Data2.data(:,1)
            ID2 = Data2.data(:,2)
            D2 = Data2.data(:,3)
            VV2 = Data2.data(:,5)
            minutes = Data2.data(:,6)
            filename = str(file_name.day2(ii))
            time2 = datenum([int(filename(4:7))*ones(size(minutes)), int(filename(8:9))*ones(size(minutes)), int(filename(10:11))*ones(size(minutes)), int(filename(12:13))*ones(size(minutes)),minutes, zeros(size(minutes))])
            PIP(ii).ID2 = []
            RT(ii).nr2 = []
            RT(ii).Emaj2 = []
            RT(ii).Emin2 = []
            RT(ii).Emaj2_calc = []
            RT(ii).Emin2_calc = []
            RT(ii).AR2 = []
            RT(ii).Dia2 = []
            RT(ii).lonX2 = []
            RT(ii).Lep2 =[]
            RT(ii).Rp2 =[]
            RT(ii).Up2 =[]
            RT(ii).Lop2 =[]
            RT(ii).or2 = []
            PIP(ii).Emaj2 = []
            PIP(ii).Emin2 = []
            PIP(ii).Emaj2_calc = []
            PIP(ii).Emaj2max_calc = []
            PIP(ii).Emin2_calc = []
            PIP(ii).AR2 = []
            PIP(ii).lonX2 = []
            PIP(ii).Len2 = []
            PIP(ii).Hig2 = []
            PIP(ii).Dia2 = []
            PIP(ii).or2 = []
            PIP(ii).D2 = []
            PIP(ii).VV2 = []
            PIP(ii).time2 = []
            # combine information from the minutes to the ID1
            ID2_min = ID2 + minutes/100 
            ID2_min = ID2_min(ID2_min>0)
            REC2  = REC2(ID2_min>0)
            [C,ia,ic] = unique(ID2_min)
            # look the particle values from the particle table  
            DELIMITER = '\t'
            HEADERLINES = 10
            xfile = str(file_name.day2(ii))
            ##################################
            # For 18 March 2014 the particle table file is split because of matlab
            # memory problems
            if strcmp(xfile, '0042014032018000_a_v_2.dat')
                xxfile = '0042014032019000_a_v_2.dat'
                # find particle tables that are between the the two velocity tables
                indx_p2 = 0
                for pindx = 1:length(flist_p2)
                    filename = flist_p2(pindx,1).name
                    particle_table.time(pindx) = datenum([int([filename(4:7)]), int(filename(8:9)), int(filename(10:11)), int(filename(12:13)), int(filename(14:15)), 0])
                    if particle_table.time(pindx) >= datenum([int(xfile(4:7)), int(xfile(8:9)), int(xfile(10:11)), int(xfile(12:13)), 0,0])&& particle_table.time(pindx)< datenum([int(xxfile(4:7)), int(xxfile(8:9)), int(xxfile(10:11)), int(xxfile(12:13)), 0,0])
                        indx_p2 = indx_p2+1
                        particle_table.name(indx_p2)= {flist_p2(pindx,1).name}
                    
                
                dataset2 = []
                for pp = 1:indx_p2
                    xpfile = str(particle_table.name(pp)) 
                    Data2 = importdata([str(subfolder_p2),'\',xpfile(1:18),'_p.dat'], DELIMITER, HEADERLINES)
                    dataset2 = [dataset2 Data2.data]
                
                dataset2(dataset2(:,6) == -1,:) = []
            else
                if event_start_time < upgrade_day
                    Data2 = importdata([str(subfolder_p2),'\',xfile(1:18),'_p.dat'], DELIMITER, HEADERLINES)
                else
                    Data2 = importdata([str(subfolder_p2),'\',xfile(1:18),'_p.dat'], DELIMITER, HEADERLINES)
                
    
            for jj = 1:size(ia)
                PIP(ii).ID2 = [PIP(ii).ID2 floor(C(jj))]
                RT(ii,jj).nr2 = [REC2((ID2_min == C(jj)))]
                for kk = 1:size(RT(ii,jj).nr2,1)
                    if strcmp(xfile, '0042014032018000_a_v_2.dat')
                        test = find(RT(ii,jj).nr2(kk) == dataset2(:,1))
                    else
                        test = find(RT(ii,jj).nr2(kk) == Data2.data(:,1))
                    
                    if  isempty(test) ~= 1
                        if strcmp(xfile, '0042014032018000_a_v_2.dat')
                            RT(ii,jj).Emaj2(kk) = [dataset2(test,18)]
                            RT(ii,jj).Emin2(kk) = [dataset2(test,19)]
                            RT(ii,jj).AR2(kk) = [dataset2(test,20)]
                            RT(ii,jj).lonX2(kk) = [dataset2.data(test,22)]
                            RT(ii,jj).Dia2(kk) = [dataset2.data(test,27)]
                            RT(ii,jj).or2(kk) = [dataset2.data(test,23)]
                            RT(ii,jj).Lep2(kk) = [dataset2.data(test,28)]
                            RT(ii,jj).Rp2(kk) = [dataset2.data(test,29)]
                            RT(ii,jj).Up2(kk) = [dataset2.data(test,30)]
                            RT(ii,jj).Lop2(kk) = [dataset2.data(test,31)]
                        else
                            RT(ii,jj).Emaj2(kk) = [Data2.data(test,18)]
                            RT(ii,jj).Emin2(kk) = [Data2.data(test,19)]
                            RT(ii,jj).AR2(kk) = [Data2.data(test,20)]
                            RT(ii,jj).lonX2(kk) = [Data2.data(test,22)]
                            RT(ii,jj).Dia2(kk) = [Data2.data(test,27)]
                            RT(ii,jj).or2(kk) = [Data2.data(test,23)]
                            RT(ii,jj).Lep2(kk) = [Data2.data(test,28)]
                            RT(ii,jj).Rp2(kk) = [Data2.data(test,29)]
                            RT(ii,jj).Up2(kk) = [Data2.data(test,30)]
                            RT(ii,jj).Lop2(kk) = [Data2.data(test,31)]
                        
                        # Derive the ellipse properties from the bounding box
                        # define bounding box dimensions
                        Wmax = 0.1*nonzeros(RT(ii,jj).Rp2(kk)-RT(ii,jj).Lep2(kk))
                        Hmax = 0.1*nonzeros(RT(ii,jj).Lop2(kk)-RT(ii,jj).Up2(kk))
                        Dim = [Wmax/2 Hmax/2]
                        Orient = RT(ii,jj).or2(kk)*pi/180
                        ini_Emaj = RT(ii,jj).Emaj2(kk)/2 
                        ini_Emin = RT(ii,jj).Emin2(kk)/2 
                        options = optimset('MaxFunEvals',1000)
                        for hh = 1:length(Orient)
                            if Orient(hh) <= pi/2
                                x0 = [sqrt(ini_Emaj(hh)) sqrt(ini_Emin(hh)) atan((ini_Emin(hh)*cos(Orient(hh)))/(ini_Emaj(hh)*sin(Orient(hh)))) atan((-ini_Emaj(hh)*sin(Orient(hh)))/(ini_Emaj(hh)*cos(Orient(hh))))]
                            else
                                Orient(hh) = pi-Orient(hh)
                                x0 =  [sqrt(ini_Emaj(hh)) sqrt(ini_Emin(hh)) atan((ini_Emin(hh)*cos(Orient(hh)))/(ini_Emaj(hh)*sin(Orient(hh)))) atan((-ini_Emaj(hh)*sin(Orient(hh)))/(ini_Emaj(hh)*cos(Orient(hh))))]
                            
                            if Orient(hh) == 0 || Orient(hh) == pi/2 || Orient(hh) == pi
                                a(hh) = Dim(hh,1)
                                b(hh) = Dim(hh,2)
                            else
                                [x, fval, failed] = fsolve(@(x)ellipsefun(x,Dim(hh,:),Orient(hh)), x0, options)
                                if failed == 0
                                    a(hh) = 0
                                    b(hh) = 0
                                else
                                    a(hh) = x(1)^2
                                    b(hh) = x(2)^2
                                
                            
                        
                        RT(ii,jj).Emaj2_calc(kk) = max([a b])
                        RT(ii,jj).Emin2_calc(kk) = min([a b])
                    
                
                PIP(ii).Emaj2 = [PIP(ii).Emaj2 mean(nonzeros(RT(ii,jj).Emaj2))]
                PIP(ii).Emin2 = [PIP(ii).Emin2 mean(nonzeros(RT(ii,jj).Emin2))]
                PIP(ii).Emaj2_calc = [PIP(ii).Emaj2_calc mean(nonzeros(RT(ii,jj).Emaj2_calc))]
                PIP(ii).Emin2_calc = [PIP(ii).Emin2_calc mean(nonzeros(RT(ii,jj).Emin2_calc))]
                if isempty(find(RT(ii,jj).Emaj2_calc))==1
                    PIP(ii).Emaj2max_calc = [PIP(ii).Emaj2max_calc NaN]
                else
                    PIP(ii).Emaj2max_calc = [PIP(ii).Emaj2max_calc max(nonzeros(RT(ii,jj).Emaj2_calc))]           
                
                PIP(ii).AR2 = [PIP(ii).AR2 mean(nonzeros(RT(ii,jj).AR2))]
                PIP(ii).Len2 = [PIP(ii).Len2 mean(nonzeros(RT(ii,jj).Rp2-RT(ii,jj).Lep2))]
                PIP(ii).Hig2 = [PIP(ii).Hig2 mean(nonzeros(RT(ii,jj).Lop2-RT(ii,jj).Up2))]
                PIP(ii).lonX2 = [PIP(ii).lonX2 mean(nonzeros(RT(ii,jj).lonX2))]
                PIP(ii).Dia2 = [PIP(ii).Dia2 mean(nonzeros(RT(ii,jj).Dia2))]
                # check the quality of the iteration
                PIP(ii).Emaj2_calc(2*PIP(ii).Emaj2_calc > sqrt((PIP(ii).Len2*0.1).^2+ (PIP(ii).Hig2*0.1).^2)) = NaN 
                PIP(ii).Emin2_calc(2*PIP(ii).Emaj2_calc > sqrt((PIP(ii).Len2*0.1).^2+ (PIP(ii).Hig2*0.1).^2)) = NaN
                PIP(ii).Emaj2max_calc(2*PIP(ii).Emaj2_calc > sqrt((PIP(ii).Len2*0.1).^2+ (PIP(ii).Hig2*0.1).^2)) = NaN 
                PIP(ii).Emaj2_calc(2*PIP(ii).Emin2_calc > sqrt((PIP(ii).Len2*0.1).^2+ (PIP(ii).Hig2*0.1).^2)) = NaN
                PIP(ii).Emin2_calc(2*PIP(ii).Emin2_calc > sqrt((PIP(ii).Len2*0.1).^2+ (PIP(ii).Hig2*0.1).^2)) = NaN
                PIP(ii).Emaj2max_calc(2*PIP(ii).Emin2_calc > sqrt((PIP(ii).Len2*0.1).^2+ (PIP(ii).Hig2*0.1).^2)) = NaN
                PIP(ii).Emaj2_calc(2*PIP(ii).Emaj2max_calc > sqrt((PIP(ii).Len2*0.1).^2+ (PIP(ii).Hig2*0.1).^2)) = NaN 
                PIP(ii).Emin2_calc(2*PIP(ii).Emaj2max_calc > sqrt((PIP(ii).Len2*0.1).^2+ (PIP(ii).Hig2*0.1).^2)) = NaN
                PIP(ii).Emaj2max_calc(2*PIP(ii).Emaj2max_calc > sqrt((PIP(ii).Len2*0.1).^2+ (PIP(ii).Hig2*0.1).^2)) = NaN 
                PIP(ii).Emaj2_calc(PIP(ii).Emaj2_calc < 0) = NaN 
                PIP(ii).Emin2_calc(PIP(ii).Emaj2_calc < 0) = NaN 
                PIP(ii).Emaj2max_calc(PIP(ii).Emaj2_calc < 0) = NaN 
                PIP(ii).Emaj2_calc(PIP(ii).Emin2_calc < 0) = NaN 
                PIP(ii).Emin2_calc(PIP(ii).Emin2_calc < 0) = NaN 
                PIP(ii).Emaj2max_calc(PIP(ii).Emin2_calc < 0) = NaN 
                PIP(ii).Emaj2_calc(PIP(ii).Emaj2max_calc < 0) = NaN 
                PIP(ii).Emin2_calc(PIP(ii).Emaj2max_calc < 0) = NaN 
                PIP(ii).Emaj2max_calc(PIP(ii).Emaj2max_calc < 0) = NaN 
                
                # check the variance
                or2_var = var(nonzeros(RT(ii,jj).or2))
                or2_var(or2_var == 0) = NaN
                PIP(ii).or2 = [PIP(ii).or2 sqrt(or2_var)]
                PIP(ii).D2 = [PIP(ii).D2 mean(D2(ID2_min==C(jj)))]
                PIP(ii).VV2 = [PIP(ii).VV2 mean(VV2(ID2_min==C(jj)))]
                PIP(ii).time2 = [PIP(ii).time2 time2(ia(jj))]
            
            PIPD2 =[PIPD2 PIP(ii).D2]
            PIPV2 =[PIPV2 PIP(ii).VV2]
            PIPtime2 =[PIPtime2 PIP(ii).time2]
            PIPEmaj2 =[PIPEmaj2 PIP(ii).Emaj2']
            PIPEmin2 =[PIPEmin2 PIP(ii).Emin2']
            PIPEmaj2_calc =[PIPEmaj2_calc PIP(ii).Emaj2_calc']
            PIPEmin2_calc =[PIPEmin2_calc PIP(ii).Emin2_calc']
            PIPEmaj2max_calc =[PIPEmaj2max_calc PIP(ii).Emaj2max_calc']
            PIPAR2 =[PIPAR2 PIP(ii).AR2']
            PIPLen2 =[PIPLen2 PIP(ii).Len2']
            PIPHig2 =[PIPHig2 PIP(ii).Hig2']
            PIPlonX2 =[PIPlonX2 PIP(ii).lonX2']
            PIPDia2 =[PIPDia2 PIP(ii).Dia2']
            PIPOR2 =[PIPOR2 PIP(ii).or2']
        
        cc = setdiff(PIPD1, PIPD2) 
        if isempty(cc)
            PIPD = [PIPD1]
            PIPV = [PIPV1]
            PIPtime = [PIPtime1]
            PIPEmaj =[PIPEmaj1]
            PIPEmin =[PIPEmin1]
            PIPEmaj_calc =[PIPEmaj1_calc]
            PIPEmin_calc =[PIPEmin1_calc]
            PIPEmajmax_calc =[PIPEmaj1max_calc]
            PIPAR =[PIPAR1]
            PIPlonX =[PIPlonX1]
            PIPDia =[PIPDia1]
            PIPOR =[PIPOR1]
            PIPLen =[PIPLen1]
            PIPHig =[PIPHig1]
        else
            PIPD = [PIPD1PIPD2]
            PIPV = [PIPV1PIPV2]
            PIPtime = [PIPtime1PIPtime2]
            PIPEmaj =[PIPEmaj1PIPEmaj2]
            PIPEmin =[PIPEmin1PIPEmin2]
            PIPEmaj_calc =[PIPEmaj1_calcPIPEmaj2_calc]
            PIPEmin_calc =[PIPEmin1_calcPIPEmin2_calc]
            PIPEmajmax_calc =[PIPEmaj1_calcPIPEmaj2_calc]
            PIPAR =[PIPAR1PIPAR2]
            PIPlonX =[PIPlonX1PIPlonX2]
            PIPDia =[PIPDia1PIPDia2]
            PIPOR =[PIPOR1PIPOR2]
            PIPLen =[PIPLen1PIPLen2]
            PIPHig =[PIPHig1PIPLen2]
        
    else
        PIPD = [PIPD1]
        PIPV = [PIPV1]
        PIPtime = [PIPtime1]
        PIPEmaj =[PIPEmaj1]
        PIPEmin =[PIPEmin1]
        PIPEmaj_calc =[PIPEmaj1_calc]
        PIPEmin_calc =[PIPEmin1_calc]
        PIPEmajmax_calc =[PIPEmaj1max_calc]
        PIPAR =[PIPAR1]
        PIPlonX =[PIPlonX1]
        PIPDia =[PIPDia1]
        PIPOR =[PIPOR1]
        PIPLen =[PIPLen1]
        PIPHig =[PIPHig1]
    
    
    # Velocity sanity - check
    Ur_Atlas = abs(9.65-10.3*exp(-0.6*PIPD))
    PIPD(PIPV>(1.5*Ur_Atlas)) = NaN
    PIPEmaj(PIPV>(1.5*Ur_Atlas)) = NaN
    PIPEmin(PIPV>(1.5*Ur_Atlas)) = NaN
    PIPEmaj_calc(PIPV>(1.5*Ur_Atlas)) = NaN
    PIPEmin_calc(PIPV>(1.5*Ur_Atlas)) = NaN
    PIPEmajmax_calc(PIPV>(1.5*Ur_Atlas)) = NaN
    PIPAR(PIPV>(1.5*Ur_Atlas)) = NaN
    PIPlonX(PIPV>(1.5*Ur_Atlas)) = NaN
    PIPDia(PIPV>(1.5*Ur_Atlas)) = NaN
    PIPOR(PIPV>(1.5*Ur_Atlas)) = NaN
    PIPLen(PIPV>(1.5*Ur_Atlas)) = NaN
    PIPHig(PIPV>(1.5*Ur_Atlas)) = NaN
    PIPV(PIPV>(1.5*Ur_Atlas)) = NaN
    
    # Take velocity values below 0.5 m/s away
    PIPD(PIPV<0.5) = NaN
    PIPEmin_calc(PIPV<0.5) = NaN
    PIPEmaj_calc(PIPV<0.5) = NaN
    PIPEmajmax_calc(PIPV<0.5) = NaN
    PIPEmaj(PIPV<0.5) = NaN
    PIPEmin(PIPV<0.5) = NaN
    PIPAR(PIPV<0.5) = NaN
    PIPlonX(PIPV<0.5) = NaN
    PIPDia(PIPV<0.5) = NaN
    PIPOR(PIPV<0.5) = NaN
    PIPLen(PIPV<0.5) = NaN
    PIPHig(PIPV<0.5) = NaN
    PIPV(PIPV<0.5) = NaN
    
    # Diameter check
    #   PIPV(PIPD<0.375) = NaN
    #   PIPEmaj(PIPD<0.375) = NaN
    #   PIPEmin(PIPD<0.375) = NaN
    #   PIPEmaj_calc(PIPD<0.375) = NaN
    #   PIPEmin_calc(PIPD<0.375) = NaN
    #   PIPAR(PIPD<0.375) = NaN
    #   PIPlonX(PIPD<0.375) = NaN
    #   PIPDia(PIPD<0.375) = NaN
    #   PIPOR(PIPD<0.375) = NaN
    #   PIPD(PIPD<0.375) = NaN







