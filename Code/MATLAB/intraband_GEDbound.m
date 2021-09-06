% do it on the first band...
data = '../../Data/cleandata.mat';
bands = STN_freqbands;
no_bands = max(bands);
%intrabands = cell(4, no_bands);

% Q: Threshold mini clusters?
for i=1:no_bands
    min_freq = min(frex(bands==i));
    max_freq = max(frex(bands==i));
    n_frex = 50;

    [evals, evecs] = freq_GED(data, 3, 'STN', min_freq, max_freq, n_frex);
    [E_CorMat, freqbands, avecorcoef] = GEDBounds(evecs, min_freq, max_freq, n_frex);
end
