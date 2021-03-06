function [E_CorMat, freqbands, avecorcoef] = GEDBounds(evecs, f_lo, f_high, n_frex)

frex = logspace(log10(f_lo), log10(f_high), n_frex); %linspace(f_lo, f_high, n_frex);
E = zscore(evecs, [], 2);
E_CorMat = (E*E'/(size(evecs, 2)-1)).^2;
nepsis = 50;
epsis  = linspace(.001, .05, nepsis);
qvec   = nan(nepsis, 1);

% find optimal epsilon
for epi=1:length(epsis)
    % scan
    % k in clustering?
    freqbands = dbscan(E_CorMat, epsis(epi), 3, 'Distance', 'Correlation');
    if max(freqbands) < 4
        continue; 
    end
    
    % compute q
    qtmp = zeros(max(freqbands), 1);
    MA = false(size(E_CorMat));
    
    for i=1:max(freqbands)
        M = false(size(E_CorMat));
        M(freqbands==i, freqbands==i) = 1;
        qtmp(i) = mean(mean(E_CorMat(M))) / mean(mean(E_CorMat(~M)));
        MA = MA+M;
    end
    qvec(epi) = mean(qtmp) + log(mean(MA(:)));
end

% run it again on the best epsilon value
[~,epsiidx] = findpeaks(qvec, 'NPeaks', 1, 'SortStr', 'descend');

if isempty(epsiidx)
    epsiidx = round(nepsis/2);
end


freqbands = dbscan(E_CorMat, epsis(epsiidx), 3, 'Distance', 'Correlation');

newc = cell(4,1); n=1;
for i=1:max(freqbands)
cc = bwconncomp(freqbands==i);
    for ci=1:cc.NumObjects
        if length(cc.PixelIdxList{ci})>2
            newc{n} = cc.PixelIdxList{ci};
            n = n+1;
        end
    end
end

freqbands = -ones(size(frex));
for ni=1:n-1
    freqbands(newc{ni}) = ni;
end

%% average correlation coefficient within each cluster

avecorcoef = zeros(max(freqbands),2);
for i=1:max(freqbands)
    submat = E_CorMat(freqbands==i,freqbands==i);
    avecorcoef(i,1) = mean(nonzeros(tril(submat,-1)));
    avecorcoef(i,2) = mean(frex(freqbands==i));
end
bw = zeros(1, max(freqbands));

figure(1), clf, colormap bone
imagesc(1-E_CorMat), hold on
f2u = round(linspace(1, length(frex), 10));
set(gca, 'clim', [.2 1], 'xscale', 'linear', 'yscale', 'linear', 'xtick', f2u, 'xticklabel', round(frex(f2u), 1), 'ytick', f2u, 'yticklabel', round(frex(f2u), 1))
axis square, axis xy

for i=1:max(freqbands)
    
    tbnds = frex(freqbands==i);
    M = max(tbnds);
    m = min(tbnds);
    bw(i) = M-m;
    tbnds = dsearchn(frex',tbnds([1 end])');
    
    % box
    plot(tbnds,[1 1]*tbnds(1), 'm', 'linew', 2)
    plot(tbnds,[1 1]*tbnds(2), 'm', 'linew', 2)
    plot([1 1]*tbnds(1),tbnds, 'm', 'linew', 2)
    plot([1 1]*tbnds(2),tbnds, 'm', 'linew', 2)

end
 
figure(2)
plot(bw, 'o-', 'linewidth', 2)
xlabel('Cluster Label')
ylabel('Bandwidth [Hz]')
% normalize within band spectrum by its components?

% fold validation split test and training
% split the bins
% MC(jackknifing?) -> random replacement
% add regularization (1%)