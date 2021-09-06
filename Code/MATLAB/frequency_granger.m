% bivariate spectral granger
clear; clc;
load('../../Data/cleandata.mat')
EEG.data = double(EEG.data);

% downsample data to 250 Hz -> maintain spectral resolution
data = downsample(EEG.data, 4);
srate = EEG.srate/4;

% frequencies of interest
lo_freq = 2;
hi_freq = 30;
numfrex = 40;
frex = linspace(lo_freq, hi_freq, numfrex);

% window and time shift
window = 200;
t_windowed = size(data, 2)-window;
% order -> we want 250 ms~ sampling rate - 250 Hz -> 1/250 = 4 ms 250/4 = 62.5 -> order of 62? ask mike
order = 25;

xGy_f = zeros(length(frex), size(data,2) - window);
yGx_f = zeros(length(frex), size(data,2) - window);

