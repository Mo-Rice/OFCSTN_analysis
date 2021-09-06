clear; clc; close all;

load('cleandata.mat');

% trial times
t =  -500:50:2000;
t_i = dsearchn(EEG.times', t');
baseidx = dsearchn(EEG.times',[-500 -200]');

% morlet frequency and cycle space
f_x = logspace(log10(1),log10(100),100); 
nCycles = logspace(log10(3),log10(8),100); 

% morlet
s_rate = 1000; % Hz
wavetime  = -2:1/EEG.srate:2; % s
half_wave = floor(length(wavetime)-1)/2;

nWave = length(wavetime);
nData = EEG.pnts*EEG.trials;
nConv = nWave+nData-1;

cmwX = zeros(length(f_x), nConv);
for fi=1:length(f_x)
    s          = nCycles(fi) / (2*pi*f_x(fi));
    cmw        = exp(1i*2*pi*f_x(fi).*wavetime) .* exp( (-wavetime.^2) ./ (2*s^2) );
    tempX      = fft(cmw,nConv);
    cmwX(fi,:) = tempX ./ max(tempX);
end

cond = EEG.conditions < 3; 

tf_all  = zeros(2, EEG.nbchan, length(f_x),length(t));

 for chani=1:EEG.nbchan
 
     dataX = fft(reshape(EEG.data(chani,:,:),1,[]) ,nConv);
 
     for fi=1:length(f_x)
         % wavelet convolution
         as = ifft(cmwX(fi,:).*dataX);
         as = as(half_wave+1:end-half_wave);
         as = abs(reshape(as,EEG.pnts, EEG.trials) ).^2;
         basepow = mean( mean(as(baseidx(1):baseidx(2),:),2) ,1);
         
         for coni=1:2
             % power values averaged over trials and baseline normalized
             % separated into all/go/no-go trials
             tf_all(coni, chani, fi, :) = 10*log10( mean(as(t_i, cond==(coni-1)), 2) ./ basepow );
         end
     end
 end
 
 save('tf_all.mat', 'tf_all','cond', 't','f_x', 'OFCchanidx','STNchanidx')
 
 