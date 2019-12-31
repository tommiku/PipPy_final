# -*- coding: utf-8 -*-

def PIPVelRel_GPM(event_start_time, event_end_time,D_PIP,PIPtime,PIPD,PIPV,n,foldername):

    # Define the time vector over the averaged time
    min_start = floor(str2num(datestr(event_start_time, 'MM'))/n)
    min_end = floor(str2num(datestr(event_end_time, 'MM'))/n)
    time_start = datenum([str2num(datestr(event_start_time, 'yyyy')), str2num(datestr(event_start_time, 'mm')), str2num(datestr(event_start_time, 'dd')),str2num(datestr(event_start_time,'HH')),(min_start*n), 0])
    time_end = datenum([str2num(datestr(event_end_time, 'yyyy')),str2num(datestr(event_end_time, 'mm')),str2num(datestr(event_end_time, 'dd')),str2num(datestr(event_end_time,'HH')),(min_end*n), 0])
    time_vector = time_start:n/(24*60):time_end
    
    mkdir([foldername,'\VelRel'])
    
    
    # Velocity factors
    avel_PIP_vel = zeros(1,size(time_vector,2)-1)
    bvel_PIP_vel = zeros(1,size(time_vector,2)-1)
    
    diff_d = diff(D_PIP)
    for aa = 1:size(time_vector,2)-1
        PIP_vel(aa).vel = []
        PIP_vel(aa).d_v = []
        PIP_vel(aa).vel_kde = []
        vel = []
        
        d = []
        indx = find(PIPtime > time_vector(aa)+30/(60*60*24) & PIPtime <= time_vector(aa+1)+30/(60*60*24))
        if isempty(indx) || isempty(find(isnan(PIPD(indx))==0))==1 || isempty(find(isnan(PIPV(indx))==0))==1 
            # relation for velocity
            avel_PIP_vel(aa) = NaN 
            bvel_PIP_vel(aa) = NaN  
        else
            vel = PIPV(indx)
            d = PIPD(indx)
            data_vel = [d vel]
            # KDE
            # The data is 1D
            if all(d==d(1)) == 1
                [density_v,xi] = ksdensity(vel)
                vq = interp1(xi,density_v,vel)
                # Search bin-wise the particles that are inside the 50# of
                # the maximum probability for diameter
                for bb = 1:size(D_PIP,2)
                    if bb == 1
                        indx_data = find(d <=(D_PIP(1)+diff_d(1)/2))
                    elseif bb == size(D_PIP,2)
                        indx_data = find(d >(D_PIP(1)-diff_d(end)/2))
                    else
                        indx_data = find(d >(D_PIP(bb)-diff_d(bb-1)/2) & d <=(D_PIP(bb)+diff_d(bb)/2)) 
                    end
                    data_v = vel(indx_data) 
                    data_d = d(indx_data)
                    vq_data = vq(indx_data)/max(vq(indx_data))
                    PIP_vel(aa).vel = [PIP_vel(aa).vel data_v(vq_data>=0.5)]
                    PIP_vel(aa).d_v = [PIP_vel(aa).d_v data_d(vq_data>=0.5)]
                    PIP_vel(aa).vel_kde = [PIP_vel(aa).vel_kde vq_data(vq_data>=0.5)]
                end
            else
                if length(data_vel) > 30 
                    [bandwidth_v,density_v,X_v,Y_v] = kde2d(data_vel)
                    [XI_v,YI_v,ZI_v] = griddata(X_v,Y_v,density_v,data_vel(:,1),data_vel(:,2), 'nearest')
                    # Diameter bins
                    for bb = 1:size(D_PIP,2)
                        if bb == 1
                            indx_data_v = find(XI_v <=(D_PIP(1)+diff_d(1)/2))
                        elseif bb == size(D_PIP,2)
                            indx_data_v = find(XI_v >(D_PIP(1)-diff_d(end)/2))
                        else
                            indx_data_v = find(XI_v >(D_PIP(bb)-diff_d(bb-1)/2) & XI_v <=(D_PIP(bb)+diff_d(bb)/2))
                        end
                        data_v = YI_v(indx_data_v) 
                        data_dv = XI_v(indx_data_v)
                        vq_data = ZI_v(indx_data_v)/max(ZI_v(indx_data_v))
                        PIP_vel(aa).vel = [PIP_vel(aa).vel data_v(vq_data>=0.5)]
                        PIP_vel(aa).d_v = [PIP_vel(aa).d_v data_dv(vq_data>=0.5)]
                        PIP_vel(aa).vel_kde = [PIP_vel(aa).vel_kde vq_data(vq_data>=0.5)]
                    end
                end
            end
            # Define the ar-D and v-D relations
            # nonlinear regression fit 
            if numel(PIP_vel(aa).d_v) > 5
                data_dv_log = log(PIP_vel(aa).d_v)
                data_v_log = log(PIP_vel(aa).vel)
                C1 = robustfit(data_dv_log, data_v_log)
                bvel_PIP_vel(aa) = C1(2)
                avel_PIP_vel(aa) = exp(C1(1))
            end        
        end
        ###################################################################
        # Images of fall velocity, area ratio and axis ratio together
        ###################################################################
        dv_nnan = d
        vel_nnan = vel
        dv_nnan(isnan(vel_nnan)) = []
        vel_nnan(isnan(vel_nnan))= []
        vel_nnan(isnan(dv_nnan))= []
        dv_nnan(isnan(dv_nnan)) = []
        
        insitu_img = figure()
        set(insitu_img, 'Color', 'white')
        set(insitu_img, 'Position', [400 50 800 870]) # left bottom width height
                
        if numel(vel_nnan) > 30
            dscatter(dv_nnan,vel_nnan)
            num_amt = numel(vel_nnan)/5.5
            clb = colorbar
            set(clb, 'TickLabels',{num2str(round(0.1*num_amt)),num2str(round(0.2*num_amt)),num2str(round(0.3*num_amt)),num2str(round(0.4*num_amt)),num2str(round(0.5*num_amt)),num2str(round(0.6*num_amt)),num2str(round(0.7*num_amt)),num2str(round(0.8*num_amt)),num2str(round(0.9*num_amt)),num2str(round(num_amt))})
            hold on
        end
         plot(D_PIP, avel_PIP_vel(aa)*D_PIP.^bvel_PIP_vel(aa), 'k-', 'LineWidth', 4)
         hold on
         grid on
         axis([0 5 0 5])
         set(gca, 'FontSize', 20,'FontName', 'Arial')
         set(gca, 'YTick', [0, 1, 2, 3, 4, 5],'YTickLabel', {'   ','   1', '   2', '   3', '   4' ,'   5'})
         ylabel('Fall velocity [m/s]','FontSize', 20, 'FontName', 'Arial')
         title(['Time ', datestr(time_vector(aa)), ' - ',datestr(time_vector(aa+1))], 'FontSize', 22, 'FontName', 'Arial')
         annotation(insitu_img,'textbox',[0.46 0.78 0.29 0.066],'String',{['V = ',num2str(avel_PIP_vel(aa),3),'D^{', num2str(bvel_PIP_vel(aa),3),'}']},...
            'FontName', 'Arial','FontSize', 20,'FitBoxToText','off','BackgroundColor',[1 1 1])
    
         xlabel('D_{PIP} [mm]','FontSize', 20, 'FontName', 'Arial')
         fname = [foldername,'\VelRel\V_Ax_Ar_DPIP_', datestr(time_vector(aa+1), 30)]
         savefig(fname)
         close
    end
    time_vector_velrel = time_vector(2:end)

