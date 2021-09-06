function [evals_STN, evecs_STN, evals_OFC, evecs_OFC] = OFCSTN_GED(default, e_val, lo_freq, hi_freq, nfrex)


load('../../Data/cleandata.mat')
EEG.data = double(EEG.data);
plot_tog = false;

%%
% TODO: pre stimulus and stimulus time windows
% time window for covariance matrices

tidx = dsearchn(EEG.times',[0 800]');

if default==true
    % Frequency parameters
    lo_freq = 2;
    hi_freq = 100;
    nfrex = 100;
end

frex = logspace(log10(lo_freq), log10(hi_freq), nfrex);

bbcov_STN = zeros(length(STNchanidx));
bbcov_OFC = zeros(length(OFCchanidx));

for triali=1:EEG.trials
    % licking onset window?
    
    % get a data snippet
    tmpdat_STN = EEG.data(STNchanidx,tidx(1):tidx(2),triali);
    tmpdat_STN = tmpdat_STN - mean(tmpdat_STN,2);
    
    tmpdat_OFC = EEG.data(OFCchanidx,tidx(1):tidx(2),triali);
    tmpdat_OFC = tmpdat_OFC - mean(tmpdat_OFC,2);
    
    % add to covariance matrix
    bbcov_STN = bbcov_STN + tmpdat_STN*tmpdat_STN' / diff(tidx);
    bbcov_OFC = bbcov_OFC + tmpdat_OFC*tmpdat_OFC' / diff(tidx);
end

evals_STN = zeros(nfrex);
evals_OFC = zeros(nfrex);
evecs_STN = zeros(nfrex, 28);
evecs_OFC = zeros(nfrex, 32); 

for fi=1:nfrex
    
    % filter data
    fdat = filterFGx(EEG.data, EEG.srate, frex(fi), 3);
    % vector selectivity for increasing frequency?
    
    % covariance
    Scov_STN = zeros(length(STNchanidx));
    Scov_OFC = zeros(length(OFCchanidx));
    
    for triali=1:EEG.trials
        
        % get a data snippet
        tmpdat_STN = fdat(STNchanidx,tidx(1):tidx(2),triali);
        tmpdat_STN = tmpdat_STN - mean(tmpdat_STN,2);
        
        tmpdat_OFC = fdat(OFCchanidx,tidx(1):tidx(2),triali);
        tmddat_OFC = tmpdat_OFC - mean(tmpdat_OFC, 2);
        
        % add to covariance matrix
        Scov_STN = Scov_STN + tmpdat_STN*tmpdat_STN' / diff(tidx);
        Scov_OFC = Scov_OFC + tmpdat_OFC*tmpdat_OFC' / diff(tidx);
    end
    
    %% GED
    
    % eig and sort
    [V_STN,D_STN] = eig(Scov_STN,bbcov_STN);
    [V_OFC,D_OFC] = eig(Scov_OFC,bbcov_OFC);
    
    [d_STN, sidx_STN] = sort(diag(D_STN),'descend');
    [d_OFC, sidx_OFC] = sort(diag(D_OFC), 'descend');
    
    V_STN = V_STN(:, sidx_STN);
    V_OFC = V_OFC(:, sidx_OFC);
    
    evals_STN(fi, :) = d_STN(e_val);
    evals_OFC(fi, :) = d_OFC(e_val);
    evecs_STN(fi, :) = V_STN(:, e_val);
    evecs_OFC(fi, :) = V_OFC(:, e_val);
    
    %% get component time series
    
    %compts_STN = V_STN(:,1)' * reshape(fdat(STNchanidx,:,:),length(STNchanidx),[]);
    %compts_STN = reshape(compts_STN,[EEG.pnts EEG.trials]);
    
    %compts_OFC = V_OFC(:,1)' * reshape(fdat(OFCchanidx,:,:),length(OFCchanidx),[]);
    %compts_OFC = reshape(compts_OFC,[EEG.pnts EEG.trials]);
    
    % you can do analyses on this time series
    
end