clear; clc; close all;
load('tf_all.mat');

figure(1)
f_x = logspace(log10(1),log10(100),100);
for fig_i=1:length(OFCchanidx)
   subplot(6, 6, fig_i)
   contourf(t, f_x, squeeze(nanmean(tf_all(1, OFCchanidx(fig_i),:,:), 2)), 40, 'linecolor', 'none')
   set(gca,'fontsize',14, 'clim', [-2, 2], 'Yscale', 'log')
   colormap jet
   title (sprintf('Stop: OFC ch. %d', fig_i))
   
end

figure(2)
for fig_ii=1:length(STNchanidx)
   subplot(6, 6, fig_ii)
   contourf(t, f_x, squeeze(nanmean(tf_all(1, STNchanidx(fig_ii),:,:), 2)), 40, 'linecolor', 'none')
   set(gca,'fontsize',14, 'clim', [-2, 2], 'Yscale', 'log')
   colormap jet
   title (sprintf('Stop: STN ch. %d', fig_ii))
end

figure(3)
f_x = logspace(log10(1),log10(100),100);
for fig_i=1:length(OFCchanidx)
   subplot(6, 6, fig_i)
   contourf(t, f_x, squeeze(nanmean(tf_all(2, OFCchanidx(fig_i),:,:), 2)), 40, 'linecolor', 'none')
   set(gca,'fontsize',14, 'clim', [-2, 2], 'Yscale', 'log')
   colormap jet
   title (sprintf('Go: OFC ch. %d', fig_i))
   
end

figure(4)
for fig_ii=1:length(STNchanidx)
   subplot(6, 6, fig_ii)
   contourf(t, f_x, squeeze(nanmean(tf_all(2, STNchanidx(fig_ii),:,:), 2)), 40, 'linecolor', 'none')
   set(gca,'fontsize',14, 'clim', [-2, 2], 'Yscale', 'log')
   colormap jet
   title (sprintf('Go: STN ch. %d', fig_ii))
end