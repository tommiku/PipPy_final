# -*- coding: utf-8 -*-

# This function reads the PIP PSD data and defines the parameters
# or reads the data from file

def ReadPIP_PSD_GPM(event_start_time,event_end_time,n,foldername):
    yy_s = datestr(event_start_time, 'yyyy')
    mm_s = datestr(event_start_time, 'mm')
    day_s = datestr(event_start_time, 'dd')
    h_s = datestr(event_start_time, 'HH')
    min_s = datestr(event_start_time, 'MM')
    s_s = datestr(event_start_time, 'SS')
    
    yy_e = datestr(event_end_time, 'yyyy')
    mm_e = datestr(event_end_time, 'mm')
    day_e = datestr(event_end_time, 'dd')
    h_e = datestr(event_end_time, 'HH')
    min_e = datestr(event_end_time, 'MM')
    s_e = datestr(event_end_time, 'SS')
    
    # File format changed 
    #upgrade_day = datenum(2014,11,23,0,0,0)
    upgrade_day = datenum(2014,01,31,0,0,0)
    if event_start_time < upgrade_day:
        D_PIP = 0.125:0.25:26
        D_PIP = [D_PIP 26]
    else:
        D_PIP = 0.1:0.2:26
        D_PIP = [D_PIP 26]
    # Open PIP PSD files
    file_name.day1 = {}
    file_name.day2 = {}
    # define which table to read
    if event_start_time < upgrade_day
        #folder_PIP = 'G:\GPM_overpass_data\PIP\f_1_4_DSD_Tables_ascii\'
        folder_PIP = 'M:\GPM_overpass_data\PIP'
        flist_1 = dir(fullfile(folder_PIP,['004',yy_s,mm_s,day_s,'*_a_d.dat']))
        flist_2 = dir(fullfile(folder_PIP,['004',yy_e,mm_e,day_e,'*_a_d.dat']))
        indx1 = 0
        indx2 = 0
        Data1 = []
        Data2 = []
        file_time.day1 = []
        file_time.day2 = []
        # Select files that fall within start and end times
        for findx = 1:length(flist_1)
            filename = flist_1(findx,1).name
            file_time.day1(findx) = datenum([str2num([filename(4:7)]), str2num(filename(8:9)), str2num(filename(10:11)),str2num(h_s), str2num(min_s),str2num(s_s)])
            if file_time.day1(findx) >= event_start_time and file_time.day1(findx)<=event_end_time
                indx1      = indx1+1
                file_name.day1(indx1)= {flist_1(findx,1).name}
            end
        end
            
        PIPtime1 = []
        PIPtime2 = []
        DELIMITER = '\t'
        HEADERLINES = 12
        # Import the file
        Data1 = importdata([folder_PIP,str(file_name.day1(1))], DELIMITER, HEADERLINES)
        hour_PIP = Data1.data(:,2)
        min_PIP = Data1.data(:,3)
        filename = str(file_name.day1(1))
        time1 = datenum([str2num(filename(4:7))*ones(size(min_PIP)), str2num(filename(8:9))*ones(size(min_PIP)), str2num(filename(10:11))*ones(size(min_PIP)), hour_PIP, min_PIP, zeros(size(min_PIP))])
        for jj = 1:size(min_PIP)
            PIPDSD1(jj,:) = 2*Data1.data(jj,6:110)
            PIPtime1(jj) = time1(jj)
        end
    
        # take the average DSD from PIP date
        DELIMITER = '\t'
        HEADERLINES = 8
        # Import the file
        Data_avgdsd = importdata([folder_PIP,str(file_name.day1(1))], DELIMITER, HEADERLINES)
        PIPDSD1_avg = 2*Data_avgdsd.data(1,1:105)
        PIPtime1_avd = time1(1)
    
        for findx2 = 1:length(flist_2)
            filename = flist_2(findx2,1).name
            file_time.day2(findx2) = datenum([str2num([filename(4:7)]), str2num(filename(8:9)), str2num(filename(10:11)),str2num( datestr(event_end_time-0.5*1/24, 'hh')),0,0])
            if file_time.day2(findx2) >= event_start_time and file_time.day2(findx2)<=event_end_time
                indx2      = indx2+1
                file_name.day2(indx2)= {flist_2(findx2,1).name}
            end
        end
    
        # Import the file
        DELIMITER = '\t'
        HEADERLINES = 12
        Data2 = importdata([folder_PIP,str(file_name.day2(1))], DELIMITER, HEADERLINES)
        hour_PIP = Data2.data(:,2)
        min_PIP = Data2.data(:,3)
        filename = str(file_name.day2(1))
        time2 = datenum([str2num(filename(4:7))*ones(size(min_PIP)), str2num(filename(8:9))*ones(size(min_PIP)), str2num(filename(10:11))*ones(size(min_PIP)), hour_PIP, min_PIP, zeros(size(min_PIP))])
        for jj = 1:size(min_PIP)
            PIPDSD2(jj,:) = 2*Data2.data(jj,6:110)
            PIPtime2(jj) = time2(jj)
        end
    
        # take the average DSD from PIP date
        DELIMITER = '\t'
        HEADERLINES = 8
        # Import the file
        Data_avgdsd = importdata([folder_PIP,str(file_name.day2(1))], DELIMITER, HEADERLINES)
        PIPDSD2_avg = 2*Data_avgdsd.data(1,1:105)
        PIPtime2_avd = time1(1)
    else
        folder_PIP = 'M:\GPM_overpass_data\PIP\f_1_4_DSD_Tables_ascii\'
        flist_1 = dir(fullfile([folder_PIP,'004',yy_s, mm_s,day_s,'*_dsd.dat']))
        flist_2 = dir(fullfile([folder_PIP,'004',yy_e, mm_e,day_e,'*_dsd.dat']))
        #flist_1 = dir(fullfile([folder_PIP,'004',yy_s, mm_s,day_s,'*_dsd_new.dat']))
        #flist_2 = dir(fullfile([folder_PIP,'004',yy_e, mm_e,day_e,'*_dsd_new.dat']))
        indx1 = 0
        indx2 = 0
        Data1 = []
        Data2 = []
        file_time.day1 = []
        file_time.day2 = []
        # Select files that fall within start and end times
        for findx = 1:length(flist_1)
            filename = flist_1(findx,1).name
            file_time.day1(findx) = datenum([str2num([filename(4:7)]), str2num(filename(8:9)), str2num(filename(10:11)),str2num(filename(12:13)), str2num(filename(14:15)),str2num(s_s)])
            if file_time.day1(findx) >= event_start_time and file_time.day1(findx)<=ceil(event_end_time)
                indx1      = indx1+1
                file_name.day1(indx1)= {flist_1(findx,1).name}
            end
         end
         PIPtime1 = []
         PIPDSD1 = []
           
         # Import the files
         for kk = 1:indx1
            DELIMITER = '\t'
            HEADERLINES = 12
            Data1 = importdata([folder_PIP,str(file_name.day1(kk))], DELIMITER, HEADERLINES)
            hour_PIP = Data1.data(:,2)
            min_PIP = Data1.data(:,3)
            filename = str(file_name.day1(1))
            time1 = datenum([str2num(filename(4:7))*ones(size(min_PIP)), str2num(filename(8:9))*ones(size(min_PIP)), str2num(filename(10:11))*ones(size(min_PIP)), hour_PIP, min_PIP, zeros(size(min_PIP))])
            for jj = 1:size(min_PIP)
                PIPDSD1 = [PIPDSD1 Data1.data(jj,6:136)]
                PIPtime1 = [PIPtime1 time1(jj)]
            end
                
             # with later software version PIP data average DSD is not
             # meaningful
                
             # Exclude the odd lines
             [jline krow] = find(PIPDSD1(:,1) == -99)
             PIPDSD1(jline,:) = []
             PIPtime1(jline) = []
             # Take out the dublicates
             [C, ia, ic] = unique(PIPtime1)
             PIPDSD1 = PIPDSD1(ia,:)
             PIPtime1 = C
         end
         
         # Take also the last file of the day
         DELIMITER = '\t'
         HEADERLINES = 12
         Data1 = importdata([folder_PIP,str(file_name.day1(end))], DELIMITER, HEADERLINES)
         hour_PIP = Data1.data(:,2)
         min_PIP = Data1.data(:,3)
         filename = str(file_name.day1(end))
         time1 = datenum([str2num(filename(4:7))*ones(size(min_PIP)), str2num(filename(8:9))*ones(size(min_PIP)), str2num(filename(10:11))*ones(size(min_PIP)), hour_PIP, min_PIP, zeros(size(min_PIP))])
         for jj = 1:size(min_PIP)
            PIPDSD1 = [PIPDSD1 Data1.data(jj,6:136)]
            PIPtime1 = [PIPtime1 time1(jj)]
         end
                
         # Exclude the odd lines
         [jline krow] = find(PIPDSD1(:,1) == -99)
         PIPDSD1(jline,:) = []
         PIPtime1(jline) = []
         # Take out the dublicates
         [C, ia, ic] = unique(PIPtime1)
         PIPDSD1 = PIPDSD1(ia,:)
         PIPtime1 = C
         PIPtime1 = PIPtime1'
            
         for findx2 = 1:length(flist_2)
            filename = flist_2(findx2,1).name
            file_time.day2(findx2) = datenum([str2num([filename(4:7)]), str2num(filename(8:9)), str2num(filename(10:11)),str2num(filename(12:13)),str2num(filename(14:15)),0])
            if file_time.day2(findx2) >= event_start_time and file_time.day2(findx2)<=ceil(event_end_time)
                indx2      = indx2+1
                file_name.day2(indx2)= {flist_2(findx2,1).name}
            end
         end
            
         PIPDSD2 = []
         PIPtime2 = []
            
            
         for kk = 1:indx2
            # Import the file
            DELIMITER = '\t'
            HEADERLINES = 12
            Data2 = importdata([folder_PIP,str(file_name.day2(kk))], DELIMITER, HEADERLINES)
            hour_PIP = Data2.data(:,2)
            min_PIP = Data2.data(:,3)
            filename = str(file_name.day2(1))
            time2 = datenum([str2num(filename(4:7))*ones(size(min_PIP)), str2num(filename(8:9))*ones(size(min_PIP)), str2num(filename(10:11))*ones(size(min_PIP)), hour_PIP, min_PIP, zeros(size(min_PIP))])
            for jj = 1:size(min_PIP)
                PIPDSD2 = [PIPDSD2 Data2.data(jj,6:136)]
                PIPtime2 = [PIPtime2 time2(jj)]
            end
    
            # Exclude the odd lines
            [jline krow] = find(PIPDSD2(:,1) == -99)
            PIPDSD2(jline,:) = []
            PIPtime2(jline) = []
            # Take out the dublicates
            [C, ia, ic] = unique(PIPtime2)
            PIPDSD2 = PIPDSD2(ia,:)
            PIPtime2 = C
         end
         # Take also the last file of the day
         DELIMITER = '\t'
         HEADERLINES = 12
         Data2 = importdata([folder_PIP,str(file_name.day2(end))], DELIMITER, HEADERLINES)
         hour_PIP = Data2.data(:,2)
         min_PIP = Data2.data(:,3)
         filename = str(file_name.day2(end))
         time2 = datenum([str2num(filename(4:7))*ones(size(min_PIP)), str2num(filename(8:9))*ones(size(min_PIP)), str2num(filename(10:11))*ones(size(min_PIP)), hour_PIP, min_PIP, zeros(size(min_PIP))])
         for jj = 1:size(min_PIP)
            PIPDSD2 = [PIPDSD2 Data2.data(jj,6:136)]
            PIPtime2 = [PIPtime2 time2(jj)]
         end
                
         # Exclude the odd lines
         [jline krow] = find(PIPDSD2(:,1) == -99)
         PIPDSD2(jline,:) = []
         PIPtime2(jline) = []
         # Take out the dublicates
         [C, ia, ic] = unique(PIPtime2)
         PIPDSD2 = PIPDSD2(ia,:)
         PIPtime2 = C
         PIPtime2 = PIPtime2'
    end
        
    # Select the data inside the event
    PIP_DSD1 = PIPDSD1(PIPtime1>= event_start_time & PIPtime1<=event_end_time,:)
    PIP_DSD2 = PIPDSD2(PIPtime2>= event_start_time & PIPtime2<=event_end_time,:)
    PIPtime1 = PIPtime1(PIPtime1>= event_start_time & PIPtime1<=event_end_time)
    PIPtime2 = PIPtime2(PIPtime2>= event_start_time & PIPtime2<=event_end_time)
    cc = setdiff(PIPtime1, PIPtime2) 
    
    if isempty(cc)
        PIP_DSD = PIP_DSD1
        PIP_time = PIPtime1
    else
        PIP_DSD = [PIP_DSD1 PIP_DSD2]
        PIP_time = [PIPtime1 PIPtime2]
    end
     
    # Define wheather the instrument is PIP or Parsivel
    [N_mean, Nt, Dm, Dmax, N0, lambda, mu, Nw, D02_exp, N02_exp, lambda2_exp, time_vector] = PSDQuantities_GPM(PIP_DSD, PIP_time, D_PIP, n)
    
    PSD_img = figure()
    set(PSD_img, 'Color', 'white')
    set(PSD_img, 'Position', [100 150 800 680]) # left bottom width height
    # Interpolate the time include all minutes before plotting
    days = floor(PIP_time(end)- PIP_time(1))
    hours = floor(((PIP_time(end)- PIP_time(1))-days)*24)
    minutes = floor((((PIP_time(end)- PIP_time(1))-days)-hours/24)*24*60)
    
    PIP_time_plot = (PIP_time(1):1/(24*60):PIP_time(end))'
    PIP_DSD_plot = zeros(size(PIP_time_plot,1),size(PIP_DSD,2)) 
    for bb = 1:size(PIP_DSD_plot,1)
        xx = find(abs((PIP_time_plot(bb)- PIP_time)) < 1e-6, 1,'first')
        if isempty(xx) ==0
            PIP_DSD_plot(bb,:) = PIP_DSD(xx,:)
        end
    end
    pcolor(PIP_time_plot,D_PIP,log10(PIP_DSD_plot'))
    whitebg(gcf,'w'), set(gcf,'Color','w')
    shading interp, lighting phong
    hold on
    set(gca,'FontWeight','normal','FontSize',14,'FontName','Ariel')
    currentMap = colormap
    newMap = [1 1 1 currentMap]
    colormap(newMap)
    axis([PIP_time_plot(1) PIP_time_plot(end) 0 10])
    set(gca, 'XTick',PIP_time_plot(1):(1/24):PIP_time_plot(end))
    datetick('x', 'keeplimits')
    colorbar
    xlabel('Time [UTC]','FontSize', 18, 'FontWeight', 'normal','FontName', 'Arial')
    ylabel('D_{PIP} [mm]','FontSize', 18, 'FontWeight', 'normal','FontName', 'Arial')
    title(['PSD with PIP ', datestr(event_start_time,31), ' - ', datestr(event_end_time,31)],'FontSize', 18, 'FontWeight', 'normal','FontName', 'Arial')
            
    fname = [foldername,'\PSDimage_', datestr(event_start_time,30),'_',datestr(event_end_time,30)]
    savefig(fname)
    close
    