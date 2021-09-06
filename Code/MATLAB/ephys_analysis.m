%% Load data

clear; close all; clc
[fnames,fpaths] = deal({});
filename = 'k';

while filename~=0
    % load in a new one
    [filename,filepath] = uigetfile('*cleandata.mat');
    fnames = cat(1,fnames,filename);
    fpaths = cat(1,fpaths,filepath);
end
fnames = fnames(1:end-1);
fpaths = fpaths(1:end-1);

%% Time-frequency analysis
for di = 1:length(fnames)

    % load LFP data
    load([fpaths{di} fnames{di}]);
    LFP = EEG;
    
    times2save = -500:50:2000;
    times2saveidx = dsearchn(LFP.times',times2save');
    baseidx = dsearchn(LFP.times',[-500 -200]');

    % specify frequencies
    frex = logspace(log10(1),log10(100),100);
    nCycles = logspace(log10(3),log10(8),100);

    % parameters for complex Morlet wavelets
    srate = 1000;
    wavtime  = -2:1/LFP.srate:2;
    half_wav = floor(length(wavtime)-1)/2;

    % FFT parameters
    nWave = length(wavtime);
    nData = LFP.pnts*LFP.trials;
    nConv = nWave+nData-1;

    % create wavelets
    cmwX = zeros(length(frex),nConv);
    for fi=1:length(frex)
        s          = nCycles(fi) / (2*pi*frex(fi));
        cmw        = exp(1i*2*pi*frex(fi).*wavtime) .* exp( (-wavtime.^2) ./ (2*s^2) );
        tempX      = fft(cmw,nConv);
        cmwX(fi,:) = tempX ./ max(tempX);
    end

    % run convolution to extract time-frequency power
    tf = zeros(LFP.nbchan,length(frex),length(times2save));

    for chani=1:LFP.nbchan

        dataX = fft(reshape(LFP.data(chani,:,:),1,[]) ,nConv);

        for fi=1:length(frex)

            % run convolution
            as = ifft(cmwX(fi,:).*dataX);
            as = as(half_wav+1:end-half_wav);
            as = abs( reshape(as,LFP.pnts,LFP.trials) ).^2;
            
            basepow = mean( mean(as(baseidx(1):baseidx(2),:),2) ,1);

            % power values averaged over trials and baseline normalized
            tf(chani,fi,:) = 10*log10( mean(as(times2saveidx,:),2) ./ basepow );
            % separate into go/no-go trials

        end
    end

    % save the power values
    save([ fpaths{di} 'tf.mat' ],'tf','times2save','frex','OFCchanidx','STNchanidx')

    %% Event-related potential

    % compute erp
    ERP_OFC = mean(LFP.data(OFCchanidx,:,:),3);
    ERP_OFC = eegfilt(ERP_OFC,LFP.srate,0,20);
    ERP_STN = mean(LFP.data(STNchanidx,:,:),3);
    ERP_STN = eegfilt(ERP_STN,LFP.srate,0,20);

    % find time indices for the relevant part of the ERP
    ERP_bounds = [-200 2000];
    time_idx = dsearchn(LFP.times',ERP_bounds');
    timeforERP = LFP.times(time_idx(1):time_idx(2));

    % plot the ERPs for each channel from OFC
    figure(1); 
    subplot(211);
    plot(timeforERP,ERP_OFC(:,time_idx(1):time_idx(2))); hold on;
    set(gca,'FontName','Times New Roman','fontsize',12)
    xlim(ERP_bounds)
    xlabel('time (ms) relative to go cue')
    ylabel('amplitude (mV)')
    title('ERPs for all OFC channels')
    ylim([-1 1]*max(max(abs(ERP_OFC(:,time_idx(1):time_idx(2)))))*1.1)

    % plot the ERPs for each channel from STN
    subplot(212);
    plot(timeforERP,ERP_STN(:,time_idx(1):time_idx(2))); hold on;
    set(gca,'FontName','Times New Roman','fontsize',12)
    xlim(ERP_bounds)
    xlabel('time (ms) relative to go cue')
    ylabel('amplitude (mV)')
    title('ERPs for all STN channels')
    ylim([-1 1]*max(max(abs(ERP_STN(:,time_idx(1):time_idx(2)))))*1.1)

    % save figure
    cd (fpaths{di})
    saveas(gcf,'ERP.jpg');

    %% Make time-frequency plots

    % plot time-frequency plot for OFC channels 
    figure(2);
    subplot(121);
    contourf(times2save,frex,squeeze(nanmean(tf(OFCchanidx,:,:),1)),40,'linecolor','none')
    set(gca,'FontName','Times New Roman','fontsize',14)
    xlabel('time (ms) relative to go cue'), ylabel('frequency (Hz)')
    title('Total power OFC')
    axis square
    colormap jet
    cb = colorbar;
    ylabel(cb,'power (dB)');
    cb.FontSize = 16;
    caxis([-1 1]*max(max(max(abs(tf(OFCchanidx,:,:)))))*.5)

    % plot time-frequency plot for STN channels 
    subplot(122);
    contourf(times2save,frex,squeeze(nanmean(tf(STNchanidx,:,:),1)),40,'linecolor','none')
    set(gca,'FontName','Times New Roman','fontsize',14)
    xlabel('time (ms) relative to go cue'), ylabel('frequency (Hz)')
    title('Total power STN')
    axis square
    colormap jet
    cb = colorbar;
    ylabel(cb,'power (dB)');
    cb.FontSize = 16;
    caxis([-1 1]*max(max(max(abs(tf(STNchanidx,:,:)))))*.5)

    % save figure
    cd (fpaths{di})
    saveas(gcf,'POWER.jpg')

end

%% End