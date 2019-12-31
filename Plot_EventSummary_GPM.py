# -*- coding: utf-8 -*-

fsize_title = 16
fsize_label = 11
fsize_axes = 11
fname = 'Arial'
linew = 2
msize = 6

# Choose the relation
if diacorr == 0.9:
    accum = accum_PIP_fac_MH05_maxmaxD_corrconst08
elif diacorr == 0.82:
    accum = accum_PIP_fac_MH05_maxmaxD_corrconst06
elif diacorr == 0.7:
    accum = accum_PIP_fac_MH05_maxmaxD_corrconst04
elif diacorr == 1 and flag.version == 1:
    accum = accum_PIP_fac_B_maxmaxD
elif diacorr == 1 and flag.version == 2:
    accum = accum_PIP_fac_KC05_maxmaxD
elif diacorr == 1 and flag.version == 3:
    accum = accum_PIP_fac_HW_maxmaxD
else:
    accum = accum_PIP_fac_MH05_maxmaxD


#PL200_rr_mtv = interp1(tv_PL200,PL200_rr,master_time_vector)
# Limit the observations when the precipitation rate is >0.2 mm/h from
#ind_time_den_rr = find(PL200_rr_mtv<0.2)
ind_vel = find(avel_mtv_max<=0 | bvel_mtv_max<=0)

fig_summary = figure('units','normalized','outerposition',[0 0 1 1])
set(fig_summary, 'Color', 'white')

# Create subplot with LWE accumulation
subplot1 = subplot('Position', [0.11 0.82 0.8 0.12],'FontSize',fsize_axes,'FontName',fname,'FontWeight', 'normal')
axis([master_time_vector(1) master_time_vector() 0 8])
box(subplot1,'on')
grid(subplot1,'on')
hold(subplot1,'all')
PL200_acc_plot  = PL200_acc-min(PL200_acc)
PL200_acc_plot(PL200_acc_plot<0) = NaN
accum_plot = nancumsum(accum)
ytl = 0:2:ceil(max([max(PL200_acc_plot) max(accum_plot)]))+2
plot1 = plot(tv_PL200,PL200_acc_plot)
set(plot1,'LineStyle','-', 'Color','k','LineWidth',linew)
ax1 = gca
set(ax1,'YLim', [ytl(1) ytl()],'YTick', ytl, ...
    'YTickLabel', cellstr(num2str(ytl')), 'YMinorTick', 'off', 'XLim', [master_time_vector(1) master_time_vector()], 'XTick',master_time_vector(1):30/(24*60):master_time_vector(), 'XTickLabel', {},...
    'FontName', fname, 'FontSize', fsize_axes,'FontWeight', 'normal')
set(get(ax1,'Ylabel'),'String','LWE [mm]', 'FontName', fname, 'FontSize', fsize_label,'FontWeight', 'normal') 

LWE_e = plot(master_time_vector,  accum_plot, '-', 'Color', [0 0.45 0.74], 'LineWidth',linew)

lg = leg('Pluvio^2 200', 'PIP')
set(lg, 'Location', 'NorthWest','FontName', fname, 'FontSize', fsize_axes,'FontWeight', 'normal')
title(['Event ', datestr(master_time_vector(1),0), ' - ',datestr(master_time_vector(),0)],'FontName', fname, 'FontSize', fsize_title, 'FontWeight', 'normal')

subplot2 = subplot('Position', [0.11 0.68 0.8 0.12],'FontSize',fsize_axes,'FontName',fname,'FontWeight', 'normal')
axis([master_time_vector(1) master_time_vector() 0 0.2])
box(subplot2,'on')
grid(subplot2,'on')
hold(subplot2,'all')

#amass_mtv(ind_time_den_rr) = NaN
amass_mtv(ind_vel) = NaN

plot2 = plot(master_time_vector, amass_mtv)
set(plot2,'LineStyle','-', 'Color','k','LineWidth',linew)
ax2 = gca
set(ax2,'YLim', [0,0.2],'YTick', [0 0.05 0.1 0.15 0.2], ...
    'YTickLabel', {'   ', '0.05', '0.1', '0.15', '0.2'}, 'YMinorTick', 'on', 'XLim', [master_time_vector(1) master_time_vector()], 'XTick',master_time_vector(1):30/(24*60):master_time_vector(), 'XTickLabel', {},...
    'FontName', fname, 'FontSize', fsize_axes,'FontWeight', 'normal')
set(get(ax2,'Ylabel'),'String','a_m [gcm^{b_m}]', 'FontName', fname, 'FontSize', fsize_label,'FontWeight', 'normal') 
lg = leg('a_{m}')
set(lg, 'Location', 'NorthWest','FontName', fname, 'FontSize', fsize_axes,'FontWeight', 'normal')


subplot3 = subplot('Position', [0.11 0.54 0.8 0.12],'FontSize',fsize_axes,'FontName',fname,'FontWeight', 'normal')
axis([master_time_vector(1) master_time_vector() 0 3.5])
box(subplot3,'on')
grid(subplot3,'on')
hold(subplot3,'all')
#bmass_mtv(ind_time_den_rr) = NaN
bmass_mtv(ind_vel) = NaN
plot3 = plot(master_time_vector, bmass_mtv)
set(plot3,'LineStyle','-', 'Color','k','LineWidth',linew)
ax3 = gca
set(ax3,'YLim', [1.5,3.5],'YTick', [1.5 2.0 2.5 3.0 3.5], ...
            'YTickLabel', {'   ', '2.0', '2.5', '3.0', '3.5'}, 'YMinorTick', 'on', ...
           'XLim', [master_time_vector(1) master_time_vector()], 'XTick',master_time_vector(1):60/(24*60):master_time_vector(), 'XTickLabel', {},...
           'FontName', fname, 'FontSize', fsize_axes,'FontWeight', 'normal')        
set(get(ax3,'Ylabel'),'String','b_m [ ]', 'FontName', fname, 'FontSize', fsize_label,'FontWeight', 'normal')
lg = leg('b_{m}')
set(lg, 'Location', 'NorthWest','FontName', fname, 'FontSize', fsize_axes,'FontWeight', 'normal')

subplot4 = subplot('Position', [0.11 0.40 0.8 0.12],'FontSize',fsize_axes,'FontName',fname,'FontWeight', 'normal')
axis([master_time_vector(1) master_time_vector() 0 6])
box(subplot4,'on')
grid(subplot4,'on')
hold(subplot4,'all')

plot(master_time_vector, log10(N02_exp_mtv), 'LineStyle', '-', 'Color','k','LineWidth',linew)
set(gca, 'YLim', [0,8],'YTick', [0 2 4 6 8], ...
    'YTickLabel', {' ', '2', '4', '6', '8'}, 'YMinorTick', 'on', 'XLim', [master_time_vector(1) master_time_vector()], 'XTick',master_time_vector(1):30/(24*60):master_time_vector(), 'XTickLabel', {},...
    'FontName', fname, 'FontSize', fsize_axes,'FontWeight', 'normal')
ylabel({'log_{10}(N_0)', '[log_{10}(mm^{-1}m^{-3})]'}, 'FontName', fname, 'FontSize', fsize_label,'FontWeight', 'normal') 
lg = leg('N_0 (exp PSD)')
set(lg, 'Location', 'SouthWest','FontName', fname, 'FontSize', fsize_axes,'FontWeight', 'normal')    

subplot5 = subplot('Position', [0.11 0.26 0.8 0.12],'FontSize',fsize_axes,'FontName',fname,'FontWeight', 'normal')
axis([master_time_vector(1) master_time_vector() 0 6])
box(subplot5,'on')
grid(subplot5,'on')

hold(subplot5,'all')
#[AX,IMG1,IMG2] = plotyy(PIPtime_n, D02_PIP,PIPtime_n, Dmax_PIP)
[AX,IMG1,IMG2] = plotyy(master_time_vector, D02_exp_mtv,master_time_vector, Dmax_mtv)
set(IMG1,'LineStyle','-', 'Color','k','LineWidth',linew)
set(IMG2,'LineStyle','--', 'Color',[0.3 0.3 0.3],'LineWidth',linew)
set(AX(1),'XColor','k','YColor','k','YLim', [0,6],'YTick', [0 1 2 3 4 5 6], ...
    'YTickLabel', {'    ','1', '2', '3', '4','5', '6'}, 'YMinorTick', 'on', 'XLim', [master_time_vector(1) master_time_vector()], 'XTick',master_time_vector(1):30/(24*60):master_time_vector(), 'XTickLabel', {},...
    'FontName', fname, 'FontSize', fsize_axes,'FontWeight', 'normal')
set(AX(2),'XColor','k','YColor',[0.3 0.3 0.3],'YLim', [0 18],'YTick', [0 3 6 9 12 15 18], ...
            'YTickLabel', {'  ',' 3', ' 6', ' 9', '12', '15', '18'}, 'YMinorTick', 'on', ...
           'XLim', [master_time_vector(1) master_time_vector()], 'XTick',master_time_vector(1):60/(24*60):master_time_vector(),'XTickLabel', {},...
           'FontName', fname, 'FontSize', fsize_axes,'FontWeight', 'normal')        
set(get(AX(1),'Ylabel'),'String','D_0 [mm]', 'FontName', fname, 'FontSize', fsize_label,'FontWeight', 'normal') 
set(get(AX(2),'Ylabel'),'String','D_{max} [mm]', 'FontName', fname, 'FontSize', fsize_label,'FontWeight', 'normal')
lg = leg([IMG1, IMG2], 'D_0 (exp PSD)', 'D_{max} of 5 min')
set(lg, 'Location', 'NorthWest','FontName', fname, 'FontSize', fsize_axes,'FontWeight', 'normal')

datetick('x','HH:MM','keeplimits')
set(get(AX(2),'Xlabel'),'String','Time [UTC]', 'FontName', fname, 'FontSize', fsize_label,'FontWeight', 'normal')

fname = [foldername,'\EventSummary_', datestr(event_start_time,30),'_',datestr(event_end_time,30)]
savefig(fname)
close

