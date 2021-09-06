function [evals, evecs] = freq_GED(data, FWHM, location, lo_freq, hi_freq, nfrex)
% More general GED function, for doing intraband analysis

load(data);
EEG.data = double(EEG.data);

tidx = dsearchn(EEG.times',[0 800]');
frex = logspace(log10(lo_freq), log10(hi_freq), nfrex);
ch_idx = nan;

if location == 'STN'
    ch_idx = STNchanidx;
    bbcov = zeros(length(ch_idx));

elseif location == 'OFC'
    ch_idx = OFCchanidx;
    bbcov = zeros(length(ch_idx));
    
else
    disp('Area not recognized, enter either STN or OFC')
    return;
end


for triali=1:EEG.trials
    
    % get a data snippet
    tmpdat = EEG.data(ch_idx,tidx(1):tidx(2),triali);
    tmpdat = tmpdat - mean(tmpdat, 2);
    
    % add to covariance matrix
    bbcov = bbcov + tmpdat*tmpdat' / diff(tidx);

end

evals = zeros(nfrex, 1);
evecs = zeros(nfrex, length(ch_idx));

for fi=1:nfrex
    
    % filter data
    fdat = filterFGx(EEG.data, EEG.srate, frex(fi), FWHM);
    % covariance
    Scov = zeros(length(ch_idx));
    
    for triali=1:EEG.trials
        
        % get a data snippet
        tmpdat = fdat(ch_idx, tidx(1):tidx(2), triali);
        tmpdat = tmpdat - mean(tmpdat,2);     
      
        % add to covariance matrix
        Scov = Scov + tmpdat*tmpdat' / diff(tidx);
    end
    
    %% GED
    
    % eig and sort
    [V, D] = eig(Scov, bbcov);
    
    [d, sidx] = sort(diag(D), 'descend');
    
    V = V(:, sidx);
    
    evals(fi, :) = d(1);
    evecs(fi, :) = V(:, 1);
    
    %% get component time series
    
    %compts_STN = V_STN(:,1)' * reshape(fdat(STNchanidx,:,:),length(STNchanidx),[]);
    %compts_STN = reshape(compts_STN,[EEG.pnts EEG.trials]);
    
    %compts_OFC = V_OFC(:,1)' * reshape(fdat(OFCchanidx,:,:),length(OFCchanidx),[]);
    %compts_OFC = reshape(compts_OFC,[EEG.pnts EEG.trials]);
    
    % you can do analyses on this time series
    
end