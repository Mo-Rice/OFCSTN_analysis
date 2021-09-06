clear; clc;
load('../../Data/cleandata.mat')

runs = 500;
min_order = zeros(runs, 1);
min_BIC = zeros(runs, 1);
EEG.data = double(EEG.data);
t_idx = dsearchn(EEG.times',[0 800]');

for r_idx=1:runs
    test = [EEG.data(OFCchanidx(randi(length(OFCchanidx), 1)), t_idx(1):t_idx(2), :); EEG.data(STNchanidx(randi(length(STNchanidx), 1)), t_idx(1):t_idx(2), :)]; % random pair between OFC/STN
    test2 = test(:, :, 1:10);

    window = 800; % ms
    orders = [1:4 5:5:30]; % time lag -> sr = 1 kHz
    BIC_full = zeros(1, length(orders));
    
    for o_idx = 1:length(orders)
        [~, E] = armorf(test2, size(test2, 3), window, orders(o_idx));
        BIC_full(o_idx) = log(det(E)) + (2^2*orders(o_idx)*log(window))/window;
    end
    
    [min_o, argmin] = min(BIC_full);
    min_BIC(r_idx) = min_o;
    min_order(r_idx) = argmin;
end

figure
subplot(121)
histogram(min_order, 15)
xlabel('k (ms)')
title('Time Lag Distribution: R = 500')

subplot(122)
histogram(min_BIC, 15)
xlabel('BIC')
title('BIC Distribution: R = 500')

% timelag of 3 seems to be the best choice

% sliding window doesn't gives poor results -> flat line
% w_s = 200
% for t_idx = 1:size(test2, 2) - window
%     temp = test2(:, t_idx:t_idx+window-1, :);
%     BIC_temp = zeros(1, length(orders));
%     
%     for o_idx = 1:length(orders)
%         [~, E]  = armorf(temp, size(temp, 3), window, orders(o_idx));
%         BIC_temp(o_idx) = log(det(E)) + (2^2*orders(o_idx)*log(w_s))/w_s;
%     end
%     
%     [~, argmin] = min(BIC_temp);
%     BIC_win(t_idx) = argmin;
% end


% figure
% % subplot(211)
% plot(orders, BIC_full, 's-', 'color', 'black', 'linewidth', 2);
% xlabel('Order [ms]')
% ylabel('BIC')

% subplot(212)
% plot(1:size(test2, 2) - window, BIC_win, 'color', 'black', 'linewidth', 2)
% xlabel('Time [ms]')
% ylabel('Optimal Order')

