
load('../../Data/cleandata.mat')
EEG.data = double(EEG.data);
plot = false;

%%

% TODO: pre stimulus and stimulus time windows
% time window for covariance matrices
tidx = dsearchn(EEG.times',[0 800]');

% Frequency parameters
lo_freq = 2;
hi_freq = 30;
numfrex = 40;

% frequencies in hz
frex = linspace(lo_freq,hi_freq,numfrex);

%% broadband covariance

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
    
%% narrowband GED

% store the maximum eigenvalues
evals_STN = zeros(numfrex,1);
evals_OFC = zeros(numfrex,1);
evecs_STN = zeros(numfrex, 28);
evecs_OFC = zeros(numfrex, 32);


for fi=1:numfrex
    
    % filter data
    fdat = filterFGx(EEG.data,EEG.srate,frex(fi),3);
    
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
    
    evals_STN(fi, :) = d_STN(1);
    evals_OFC(fi, :) = d_OFC(1);
    evecs_STN(fi, :) = V_STN(:, 1);
    evecs_OFC(fi, :) = V_OFC(:, 1);
    
    %% get component time series
    
    compts_STN = V_STN(:,1)' * reshape(fdat(STNchanidx,:,:),length(STNchanidx),[]);
    compts_STN = reshape(compts_STN,[EEG.pnts EEG.trials]);
    
    compts_OFC = V_OFC(:,1)' * reshape(fdat(OFCchanidx,:,:),length(OFCchanidx),[]);
    compts_OFC = reshape(compts_OFC,[EEG.pnts EEG.trials]);
    
    % you can do analyses on this time series
    
end


%%
% Frequency bounds
E_OFC = zscore(evecs_OFC, [], 2);
E_CorMat_OFC = (E*E'/(32-1).^2;

E_STN = zscore(evecs_STN, [], 2);
E_CorMat_STN = (E*E'/(28-1).^2;



%%


if plot == true
    %%

    figure(1), clf
    plot(frex,evals_STN(:,1),'s-', 'linewidth', 2)
    hold on
    plot(frex, evals_OFC(:,1), 's-', 'linewidth', 2)
    xlabel('Frequency (Hz)')
    ylabel('\lambda')
    legend("STN 1", "OFC")
    hold off

    %%

    figure(2), clf
    hold on

    for i=1:20
        plot(evals_STN(i, :), 's-', 'linewidth', 2, 'DisplayName', num2str(i))
    end
    ylabel('\lambda')
    xlabel('Component')
    title('STN Low Frequency')
    legend()
    hold off

    figure(3), clf
    hold on
    for i=21:40
        plot(evals_STN(i, :), 's-', 'linewidth', 2, 'DisplayName', num2str(i))
    end
    ylabel('\lambda')
    xlabel('Component')
    title('STN High Frequency')
    legend()
    hold off


    figure(4), clf
    hold on
    for i=1:20
        plot(evals_OFC(i, :), 's-', 'linewidth', 2, 'DisplayName', num2str(i))
    end
    ylabel('\lambda')
    xlabel('Component')
    title('OFC Low Frequency')
    legend()
    hold off

    figure(5), clf
    hold on

    for i=21:40
        plot(evals_OFC(i, :), 's-', 'linewidth', 2, 'DisplayName', num2str(i))
    end
    ylabel('\lambda')
    xlabel('Component')
    title('OFC High Frequency')
    legend()
    hold off

    figure(6), clf
       subplot(211)
       imagesc(evals_OFC)
       set(gca,'fontsize',14, 'clim', [0, 1])
       xlabel('\lambda')
       ylabel('Frequency [Hz]')
       title('OFC GED')
       colorbar()

       subplot(212)
       imagesc(evals_STN)
       set(gca,'fontsize',14, 'clim', [0, 1])
       xlabel('\lambda')
       ylabel('Frequency [Hz]')
       title('STN GED')
       colorbar()

    % to do
    % TF condition average, condition difference - done
    % maximum eigen value spectrum (two comps.) - done 
    % permutation testing
    % ar fitting
    % looking into phase mod delay
end