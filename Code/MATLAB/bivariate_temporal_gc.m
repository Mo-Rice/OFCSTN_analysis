% bivariate granger
% number of trials is currently set to 30 for speed purposes
clear; clc;
load('../../Data/cleandata.mat')
EEG.data = double(EEG.data);

% parameters
optimal_order = 5;
window = 200;

% randomly select channels from OFC and STN data (for now)
OFC = OFCchanidx(randi(length(OFCchanidx), 1));
STN = OFCchanidx(randi(length(OFCchanidx), 1));
data = [EEG.data(OFC, :, 1:10); EEG.data(STN, :, 1:10)];
t_windowed = size(data, 2)-window;
xGy = zeros(1, size(data, 2) - window);
yGx = zeros(1, size(data, 2) - window);
GC_time = zeros(1, size(data, 2) - window);

% For each window of data
for time_i = 1:t_windowed
    temp        = reshape(data(:,time_i:time_i+window-1,:),2,window*size(data,3));
    [A_x, E_x]  = armorf(temp(1, :), size(data, 3), window, optimal_order);
    [A_y, E_y]  = armorf(temp(2, :), size(data, 3), window, optimal_order);
    [A, E]      = armorf(temp, size(data, 3), window, optimal_order);
    xGy(time_i) = log(E_x/E(1, 1));
    yGx(time_i) = log(E_y/E(2, 2));
    GC_time(time_i) = mean(EEG.times(time_i:time_i+window-1));
end

figure
hold on
plot(GC_time, xGy);
plot(GC_time, yGx);
xlabel('Time (ms)')
ylabel('Granger Causality')
legend("OFC ch. " + num2str(OFC)+ " \rightarrow " + "STN ch. " + num2str(STN),  "STN ch. " + num2str(STN) + " \rightarrow " + "OFC ch. " + num2str(OFC))

% Frequency Domain
