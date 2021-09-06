% Taking a closer look at the eigenspectrum for the frequencies in band
% 1. need to choose an FWHM for the full spectrum
% 2. need to choose one for the frequency bands
%   * use bands from "optimal" macro
%   * then calculate aNMI for in band clustering
%   * I suppose I can do this for the large bands as well...? just calculate the aNMI for the macro scopic bands? <- to this first
%% 
data = '../../Data/cleandata.mat';

f_lo = 2;
f_hi = 100;
N_f = 100;
frex = linspace(f_lo, f_hi, N_f);
FWHMs = linspace(0.5, 4, 8);

GED_evecs = [];
freqbands_f = [];
location = "OFC";
aNMI = zeros(1, size(FWHMs, 2));
sil = [];

for FWHM=1:size(FWHMs, 2)
    [~, evecs] = freq_GED(data, FWHMs(FWHM), location, f_lo, f_hi, N_f);
    GED_evecs = [GED_evecs; evecs];
    [E_CorMat, freqbands, ~] = GEDBounds(evecs, f_lo, f_hi, N_f);
    freqbands_f = [freqbands_f; freqbands];
    figure
    [S, H] = silhouette(E_CorMat, freqbands);
    title("FWHM = " + num2str(FWHMs(FWHM)))
    sil = [sil; S];
end

M = max(max(freqbands_f));
max_sil = zeros(M, size(FWHMs, 2));

for i=1:M
    for j=1:size(FWHMs, 2)
        max_sil(i, j) = (freqbands_f(j, :) == i)*sil(((j-1)*100+1):j*100)/nnz(freqbands_f(j, :) == i);
    end
end

figure
imagesc(max_sil)
colorbar
xticklabels(FWHMs)
xlabel('FWHM (Hz)')
ylabel('Cluster Label')
title('Average Cluster Silhoutte: ' + location)


% figure
% plot(FWHMs, aNMI, '-o', 'linewidth', 2)
% xlabel('FWHM (Hz)');
% ylabel('aNMI');

% I'm not quite sure what this is telling me -> the aNMI is always going to be worse for clusters that are not produced from the same
% batch.... 
% just do permutation testing -> randomize the data and perform the analysis 
% STD over covariance matrix (Z-scoring like)
% look up a paper on spectral analysis :']
% write an abstract for the master's thesis to guide the analysis steps
% abstract to define a goal