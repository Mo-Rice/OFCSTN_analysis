clear; clc; close all;
load('cleandata.mat');
load('tf_all.mat');

figure(1)

subplot(321)
    contourf(t, f_x, squeeze(nanmean(tf_all(2, OFCchanidx,:,:), 2)), 40, 'linecolor', 'none')
    set(gca,'fontsize',14, 'clim', [-2, 2], 'Yscale', 'log')
    xlabel('Time [ms]'), ylabel(' Frequency [Hz]')
    colormap jet
    cb = colorbar;
    ylabel(cb,'Power [dB]');
    cb.FontSize = 14;
    title('OFC: Go Trials')
    

subplot(322)
    contourf(t, f_x, squeeze(nanmean(tf_all(2, STNchanidx,:,:), 2)), 40, 'linecolor', 'none')
    set(gca,'fontsize',14, 'clim', [-2, 2], 'Yscale', 'log')
    xlabel('Time [ms]'), ylabel(' Frequency [Hz]')
    colormap jet
    cb = colorbar;
    ylabel(cb,'Power [dB]');
    cb.FontSize = 14;
    title('STN: Go Trials')

subplot(323)
    contourf(t, f_x, squeeze(nanmean(tf_all(1, OFCchanidx,:,:), 2)), 40, 'linecolor', 'none')
    set(gca,'fontsize',14, 'clim', [-2, 2], 'Yscale', 'log')
    xlabel('Time [ms]'), ylabel(' Frequency [Hz]')
    colormap jet
    cb = colorbar;
    ylabel(cb,'Power [dB]');
    cb.FontSize = 14;
    title('OFC: Stop Trials')

subplot(324)
    contourf(t, f_x, squeeze(nanmean(tf_all(1, STNchanidx,:,:), 2)), 40, 'linecolor', 'none')
    set(gca,'fontsize',14, 'clim', [-2, 2], 'Yscale', 'log')
    xlabel('Time [ms]'), ylabel(' Frequency [Hz]')
    colormap jet
    cb = colorbar;
    ylabel(cb,'Power [dB]');
    cb.FontSize = 14;
    title('STN: Stop Trials')

subplot(325)
    contourf(t, f_x, (squeeze(nanmean(tf_all(1, OFCchanidx,:,:), 2))-squeeze(nanmean(tf_all(2, OFCchanidx,:,:), 2))), 40, 'linecolor', 'none')
    set(gca,'fontsize',14, 'clim', [-2, 2], 'Yscale', 'log')
    xlabel('Time [ms]'), ylabel(' Frequency [Hz]')
    colormap jet
    cb = colorbar;
    ylabel(cb,'Power [dB]');
    cb.FontSize = 14;
    title('OFC: Difference')


subplot(326)
    contourf(t, f_x, (squeeze(nanmean(tf_all(1, STNchanidx,:,:), 2))-squeeze(nanmean(tf_all(2, STNchanidx,:,:), 2))), 40, 'linecolor', 'none')
    set(gca,'fontsize',14, 'clim', [-2, 2], 'Yscale', 'log')
    xlabel('Time [ms]'), ylabel(' Frequency [Hz]')
    colormap jet
    cb = colorbar;
    ylabel(cb,'Power [dB]');
    cb.FontSize = 14;
    title('STN: Difference')
    
