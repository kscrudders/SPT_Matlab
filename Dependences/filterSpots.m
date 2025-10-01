function [spotsTable, Detector] = filterSpots(img_stack, spotsTable, Detector, plotFlag, save_dir, options)
    % WIP 20250728
    if ~isfield(Detector,'minPeak') % amplitude
        Detector.minPeak = 0;      
    end   
    if ~isfield(Detector,'minIntegrated') % net integrated intensity
        Detector.minIntegrated = 0; 
    end   
    if ~isfield(Detector,'minSNR') % signal / bg
        Detector.minSNR = 0;        
    end   
    if ~isfield(Detector,'maxSNR') % signal / bg
        Detector.maxSNR = inf;        
    end  

    culling_flag = 2; % 1 = PCA, 2 = manual stats thresholds

    switch culling_flag
        case 1 % PCA
            results = pca_cluster_spots(spotsTable);
           
            cluster1_idx = results.cluster.idx == 1;
            cluster2_idx = results.cluster.idx == 2;
    
            % Assume there are more valid spots than there are invalid ones
            % keep the larger PCA cluster
            if sum(cluster1_idx) < sum(cluster2_idx)
                spotsTable = spotsTable(cluster2_idx,:);
            else
                spotsTable = spotsTable(cluster1_idx,:);
            end
        case 2 % manual
            %{
            %------------------------------ Filter based on Sigma ------------------------------%
                data_2_threshold = double(spotsTable{:,'Fit_sigma'});
                [Y, binEdges] = histcounts(data_2_threshold, 'BinMethod', 'fd','Normalization','PDF'); % The Freedman-Diaconis rule is less sensitive to outliers in the data, and might be more suitable for data with heavy-tailed distributions. It uses a bin width of 2*iqr(X(:))*numel(X)^(-1/3).
            
                X = movmean(binEdges,2);
                X = X(2:end);
            
                sig = [5 5];
                [fitresult, lower_thres, upper_thres] = KLS_fit_1guass(X, Y, sig,0);
            
                if plotFlag == 1
                    figure()
                        histogram(data_2_threshold,'BinMethod','fd','Normalization','PDF')
                        xlabel('Guassian Sigma (px)')
                    
                        x_range = min(X):range(X)/1000:max(X);
                        hold on
                            plot(x_range, fitresult(x_range),'--','LineWidth',2,'Color','k');
                        hold off
        
                        lbl = sprintf('Mean + %.1f StdDev', sig(2));
                        xline(upper_thres, '--', lbl, 'LabelVerticalAlignment','middle', 'HandleVisibility','on');
                        lbl = sprintf('Mean - %.1f StdDev', sig(1));
                        xline(lower_thres, '--', lbl, 'LabelVerticalAlignment','middle', 'HandleVisibility','on');

                        lbl = 'Manual Threshold';
                        xline(Detector.maxSigma, '--', lbl, 'LabelVerticalAlignment','top', 'HandleVisibility','on');
        
                        lbl = 'Manual Threshold';
                        xline(Detector.minSigma, '--', lbl, 'LabelVerticalAlignment','top', 'HandleVisibility','on');
                        box off
                end
        
                %valid_idx = data_2_threshold <= upper_thres;
                %spotsTable = spotsTable(valid_idx,:);

                valid_idx = data_2_threshold > Detector.minSigma;
                spotsTable = spotsTable(valid_idx,:);

                data_2_threshold = double(spotsTable{:,'Fit_sigma'});

                valid_idx = data_2_threshold <= Detector.maxSigma;
                spotsTable = spotsTable(valid_idx,:);
        
            %}
            %------------------------------ Filter based on SNR ------------------------------%
                spot_sig = spotsTable{:,'Integrate_fit_Signal'};
                spot_bkgd = spotsTable{:,'Integrate_fit_Bkgd'};
                spot_bkgd(spot_bkgd < 1) = 1; % background < 1 is invalid
                SNR =  spot_sig ./ spot_bkgd;
            
                data_2_threshold = SNR;
                [Y, binEdges] = histcounts(data_2_threshold, 'BinMethod', 'fd','Normalization','PDF'); % The Freedman-Diaconis rule is less sensitive to outliers in the data, and might be more suitable for data with heavy-tailed distributions. It uses a bin width of 2*iqr(X(:))*numel(X)^(-1/3).
            
                X = movmean(binEdges,2);
                X = X(2:end);
            
                if plotFlag == 1 && options.Show.filterSpots_SBRHist == 1
                    figure()
                        histogram(data_2_threshold,'BinMethod','fd','Normalization','PDF')
                        xlabel('Peak/Flat Bkgd (SBR)')
                        ylabel('PDF')
                    
                        %[fitresult, lower_thres, upper_thres] = KLS_fit_1guass(X, Y, sig,0);
                        box off
        
                        %lbl = 'Manual Threshold';
                        %xline(Detector.maxSNR, '--', lbl, 'LabelVerticalAlignment','top', 'HandleVisibility','on');
        
                        lbl = 'Manual Threshold';
                        xline(Detector.minSNR, '--', lbl, 'LabelVerticalAlignment','top', 'HandleVisibility','on');

                        savefig(fullfile(save_dir, 'SBR_Threshold'))
                end
            
                valid_idx = SNR > Detector.minSNR;
                spotsTable = spotsTable(valid_idx,:);
        
                spot_sig = spotsTable{:,'Peak'};
                spot_bkgd = spotsTable{:,'Background_px'};
                spot_bkgd(spot_bkgd < 1) = 1; % background < 1 is invalid
                SNR =  spot_sig ./ spot_bkgd;
        
                %valid_idx = SNR <= Detector.maxSNR;
                %spotsTable = spotsTable(valid_idx,:);

            %}
            %{
            %------------------------------ Filter based on bgkd > signal ------------------------------%
                valid_idx = spotsTable.Integrate_fit_Bkgd < spotsTable.Integrate_fit_Signal;
                spotsTable = spotsTable(valid_idx,:);
            %}
        
    end

    % Extra manual filters that are currently not in use
    %{
    %------------------------------ Filter based on fit peak ------------------------------%
        data_2_threshold = spotsTable{:,'Peak'};
        [Y, binEdges] = histcounts(data_2_threshold, 'BinMethod', 'fd','Normalization','PDF'); % The Freedman-Diaconis rule is less sensitive to outliers in the data, and might be more suitable for data with heavy-tailed distributions. It uses a bin width of 2*iqr(X(:))*numel(X)^(-1/3).
    
        X = movmean(binEdges,2);
        X = X(2:end);
    
        sig = [3 3];
        [fitresult, lower_thres, upper_thres] = KLS_fit_1guass(X, Y, sig,0);
    
        if plotFlag == 1 && options.Show.filterSpots_peakHist == 1
            figure()
                histogram(data_2_threshold,'BinMethod','fd','Normalization','PDF')
                xlabel('Peak Intensity (units set by input image)')
            
                x_range = min(X):range(X)/1000:max(X);
                hold on
                    plot(x_range, fitresult(x_range),'--','LineWidth',2,'Color','k');
                hold off

                lbl = sprintf('Mean + %.1f StdDev', sig(2));
                xline(upper_thres, '--', lbl, 'LabelVerticalAlignment','top', 'HandleVisibility','on');
                lbl = sprintf('Mean - %.1f StdDev', sig(1));
                xline(lower_thres, '--', lbl, 'LabelVerticalAlignment','top', 'HandleVisibility','on');
                box off
        end

        valid_idx = data_2_threshold <= upper_thres;
        spotsTable = spotsTable(valid_idx,:);
    %}

    %{
    %------------------------------ Filter based on LoG Peak ------------------------------%
        data_2_threshold = spotsTable{:,'Quality_LoG'};
        [Y, binEdges] = histcounts(data_2_threshold, 'BinMethod', 'fd','Normalization','PDF'); % The Freedman-Diaconis rule is less sensitive to outliers in the data, and might be more suitable for data with heavy-tailed distributions. It uses a bin width of 2*iqr(X(:))*numel(X)^(-1/3).
    
        X = movmean(binEdges,2);
        X = X(2:end);
    
        sig = [3 3];
        [fitresult, lower_thres, upper_thres] = KLS_fit_1guass(X, Y, sig,0);
    
        if plotFlag == 1 && options.Show.filterSpots_LoGHist == 1
            figure()
                histogram(data_2_threshold,'BinMethod','fd','Normalization','PDF')
                xlabel('Spot LoG Quality')
            
                x_range = min(X):range(X)/1000:max(X);
                hold on
                    plot(x_range, fitresult(x_range),'--','LineWidth',2,'Color','k');
                hold off

                lbl = sprintf('Mean + %.1f StdDev', sig(2));
                xline(upper_thres, '--', lbl, 'LabelVerticalAlignment','top', 'HandleVisibility','on');
                lbl = sprintf('Mean - %.1f StdDev', sig(1));
                xline(lower_thres, '--', lbl, 'LabelVerticalAlignment','top', 'HandleVisibility','on');
                box off
        end

        valid_idx = data_2_threshold >= upper_thres;
        spotsTable = spotsTable(valid_idx,:);
    %}

    %{
    Extra figures if you want to visualize this intermediate results
    if plotFlag == 1
        figure()
        histogram(spotsTable{:,'Fit_sigma'}.*2.35,'BinMethod','fd')
        xlabel('Guassian FWHM (px)')
    end

    if plotFlag == 1
        figure()
        histogram(spotsTable{:,'Fit_sigma'}.*2.35 .* Detector.pixelsize .* 1000,'BinMethod','fd')
        xlabel('Guassian FWHM (nm)')
    end
    %}

% -------------------------------------------------------------------------
% 6) [Plot] Spot Overlay
if plotFlag && options.Show.detectionOverlay == 1
    figure()
    img = img_stack(:,:,options.FrameToShow);
    imshow(img,[]);
    
    hold on
        row_idx = spotsTable{:,'Frame'} == options.FrameToShow;
        % Pull out coordinates and orient as N by T
        x = spotsTable{row_idx,'x_px'};   % N by T
        y = spotsTable{row_idx,'y_px'};   % N by T
        
        % Drop tracks that are entirely NaN (never detected)
        valid = any(~isnan(x),2);            % Nx1 logical
        x = x(valid,:);
        y = y(valid,:);
        
        x = x.'; % reshape the matix so its T by N so plot reads the points correctly
        y = y.';
        
        % plot() accepts a matrix: each column (now becuase we reshaped the vector this is N) becomes one line
        scatter(x, y, 'o', 'LineWidth',2, 'Color', [0 0 1 0.4]); % one color
        %plot(x, y, 'LineWidth',2); % different colors per track
    hold off
    
    axis tight; % Remove extra whitespace

    title('Post-Localization Filter')

    savefig(fullfile(save_dir, 'Localization Check - Post-Localization Filter'))
end
end


function results = pca_cluster_spots(spotsTable, opts)
%PCA_CLUSTER_SPOTS  PCA + clustering for spot features.
%
% results = pca_cluster_spots(spotsTable, opts)
%
% Required input:
%   spotsTable : MATLAB table with columns:
%       {'SpotID','Frame','x_px','y_px','Quality_LoG','Peak','Fit_sigma','Net_Signal','Integrate_Signal','Integrate_bkgd'}
%
% Optional opts (struct, all fields optional):
%   .features        (cellstr) default = {'Quality_LoG','Peak','Fit_sigma','Background_px','Net_Signal','Integrate_fit_Signal','Integrate_fit_Bkgd',Integrate_Signal_ROI'}
%   .standardize     (char) 'robust' (default) | 'zscore' | 'none'
%   .impute          (char) 'median' (default) | 'remove'
%   .varianceTarget  (double) cumulative variance % to retain PCs, default 95
%   .cluster.method  (char) 'kmeans' (default) | 'gmm' | 'dbscan' | 'none'
%   .cluster.k       (scalar) K for kmeans/GMM (if empty, auto-select 2–6 by silhouette)
%   .cluster.rangeK  (1x2) [kmin kmax] for auto-K, default [2 6]
%   .cluster.replicates (int) default 10
%   .cluster.dbscanEps  (double) default [], auto via k-dist heuristic if empty
%   .cluster.dbscanMinPts (int) default max(5, round(0.01*N))
%   .rngSeed         (scalar) default 42
%   .plotFlag        (logical) default true
%
% Output:
%   results : struct with fields
%       .features, .X, .Xscaled, .scaleStats
%       .coeff, .score, .latent, .explained, .cumExplained, .numPCs
%       .cluster.method, .cluster.idx, .cluster.model, .cluster.metrics
%       .rowMeta (table: SpotID, Frame)
%       .figures (handles if plotFlag)
%
% Notes:
% - Clustering is done in the retained PC subspace (results.score(:,1:numPCs)).
% - Robust scaling uses median/MAD (consistent MAD = 1.4826*mad).
%
% Kevin-ready. No toolboxes beyond Statistics & Machine Learning.
%
% Author: Kevin L. Scrudders w/ ChatGPT (GPT‑5 Thinking) 20250812
% Tested: Kevin L. Scrudders 20250812 on SM data

% -------------------------- Defaults & Validate --------------------------
arguments
    spotsTable table
    % <--- use only your six features by default
    opts.features cell = {'Quality_LoG','Peak','Fit_sigma','Net_Signal', ...
                          'Integrate_fit_Signal','Integrate_fit_Bkgd'}
    opts.standardize char {mustBeMember(opts.standardize,{'robust','zscore','none'})} = 'robust'
    opts.impute char {mustBeMember(opts.impute,{'median','remove'})} = 'median'
    opts.varianceTarget (1,1) double {mustBeGreaterThan(opts.varianceTarget,0), mustBeLessThanOrEqual(opts.varianceTarget,100)} = 95
    opts.cluster struct = struct()
    opts.rngSeed (1,1) double = 42
    opts.plotFlag (1,1) logical = true
end

reqCols = {'SpotID','Frame','x_px','y_px','Quality_LoG','Peak','Fit_sigma', ...
           'Net_Signal','Integrate_fit_Signal','Integrate_fit_Bkgd'};

missingCols = setdiff(reqCols, spotsTable.Properties.VariableNames);

if ~isempty(missingCols)
    error('pca_cluster_spots:MissingColumns', ...
        'spotsTable is missing columns: %s', strjoin(missingCols, ', '));
end

% Feature presence
featCols = opts.features(:)';
featMissing = setdiff(featCols, spotsTable.Properties.VariableNames);
if ~isempty(featMissing)
    error('pca_cluster_spots:MissingFeatures', ...
        'Requested feature(s) not found: %s', strjoin(featMissing, ', '));
end


% Cluster defaults
copts = struct( ...
    'method','kmeans', ...
    'k',[], ...
    'rangeK',[2 6], ...              % (fix to match docstring)
    'replicates',10, ...
    'dbscanEps',[], ...
    'dbscanMinPts',max(5, round(0.01*height(spotsTable))), ...
    'space','whiten', ...               % NEW: 'pca' (default) | 'whiten' | 'scaled'
    'pcWeights',[] ...               % NEW: optional weights per PC (or scalar)
);

% override from user struct (only known fields)
if ~isempty(fieldnames(opts.cluster))
    userF = fieldnames(opts.cluster);
    for i=1:numel(userF)
        f = userF{i};
        if isfield(copts, f)
            copts.(f) = opts.cluster.(f);
        else
            warning('pca_cluster_spots:UnknownClusterOpt','Ignoring unknown cluster option "%s".', f);
        end
    end
end

rng(opts.rngSeed); %#ok<RAND>

% -------------------------- Extract & Clean ------------------------------
rowMeta = spotsTable(:, {'SpotID','Frame'});
X = spotsTable{:, featCols};
X = double(X);

% Handle missing
nanRows = any(isnan(X) | isinf(X), 2);
if any(nanRows)
    switch opts.impute
        case 'median'
            med = median(X, 'omitnan');
            for j=1:size(X,2)
                col = X(:,j);
                col(nanRows) = med(j);
                X(:,j) = col;
            end
        case 'remove'
            X(nanRows,:) = [];
            rowMeta(nanRows,:) = [];
    end
end

% -------------------------- Scale ---------------------------------------
switch opts.standardize
    case 'none'
        Xs = X;
        scaleStats = struct('center',zeros(1,size(X,2)), 'scale', ones(1,size(X,2)), ...
                            'mode','none', 'features',{featCols});
    case 'zscore'
        mu = mean(X,1);
        sd = std(X,0,1);
        sd(sd==0) = 1;
        Xs = (X - mu)./sd;
        scaleStats = struct('center',mu, 'scale', sd, 'mode','zscore', 'features',{featCols});
    case 'robust'
        med = median(X,1);
        madc = 1.4826*mad(X,1,1); % consistent MAD
        madc(madc==0) = 1;
        Xs = (X - med)./madc;
        scaleStats = struct('center',med, 'scale', madc, 'mode','robust', 'features',{featCols});
end

% -------------------------- PCA -----------------------------------------
% Use SVD-based PCA; determine #PCs from cumulative variance
[coeff, score, latent, ~, explained, mu] = pca(Xs, 'Algorithm','svd', 'Centered', true);
cumExpl = cumsum(explained);
numPCs  = find(cumExpl >= opts.varianceTarget, 1, 'first');
if isempty(numPCs), numPCs = size(score,2); end
scoreR  = score(:,1:numPCs);

% ------------------ Choose space for clustering  --------------------
Z = scoreR; % default: raw PCs
switch lower(copts.space)
    case 'pca'
        % nothing
    case 'whiten'
        % equalize variance across PCs so PC2 matters as much as PC1
        L = sqrt(latent(1:numPCs)).';  L(L==0) = 1;
        Z = bsxfun(@rdivide, scoreR, L);
    case 'scaled'
        % cluster in standardized feature space (all selected features)
        Z = Xs;
    otherwise
        error('pca_cluster_spots:BadSpace','Unknown cluster.space "%s".', copts.space);
end

% -------------------------- Clustering ----------------------------------
clusterIdx = []; clusterModel = []; clusterMetrics = struct();
switch lower(copts.method)
    case 'none'
        % no-op
    case 'kmeans'
        if isempty(copts.k)
            [bestK, bestIdx, bestSil] = autoKmeans(Z, copts.rangeK, copts.replicates);
            clusterIdx = bestIdx; clusterMetrics.silhouette = bestSil; clusterMetrics.k = bestK;
        else
            [clusterIdx, C] = kmeans(Z, copts.k, 'Replicates',copts.replicates, ...
                'MaxIter',1000, 'Display','off', 'OnlinePhase','on');
            clusterMetrics.k = copts.k; clusterMetrics.centroids = C;
        end
    case 'gmm'
        if isempty(copts.k)
            [bestK, GM, idx, crit] = autoGMM(Z, copts.rangeK);
            clusterIdx = idx; clusterModel = GM; clusterMetrics = crit; clusterMetrics.k = bestK;
        else
            GM = fitgmdist(Z, copts.k, 'RegularizationValue',1e-6, ...
                'CovarianceType','full', 'SharedCovariance',false, 'Replicates',copts.replicates);
            [~, clusterIdx] = max(posterior(GM, Z), [], 2);
            clusterModel = GM; clusterMetrics.k = copts.k;
        end
    case 'dbscan'
        minpts = copts.dbscanMinPts; eps = copts.dbscanEps;
        if isempty(eps), eps = heuristicEps(Z, minpts); end
        clusterIdx = dbscan(Z, eps, minpts);
        clusterMetrics.eps = eps; clusterMetrics.minpts = minpts;
    otherwise
        error('pca_cluster_spots:BadMethod','Unknown cluster method "%s".', copts.method);
end

% -------------------------- Plots ---------------------------------------
figs = struct();
if opts.plotFlag  && opts.Show.filterSpots_PCAstats == 1
    figs = makePlots(score, explained, coeff, featCols, clusterIdx, numPCs);
end

% -------------------------- Pack Results --------------------------------
results = struct( ...
    'features',{featCols}, ...
    'X', X, ...
    'Xscaled', Xs, ...
    'scaleStats', scaleStats, ...
    'coeff', coeff, ...
    'score', score, ...
    'latent', latent, ...
    'explained', explained, ...
    'cumExplained', cumExpl, ...
    'numPCs', numPCs, ...
    'cluster', struct('method', lower(copts.method), 'idx', clusterIdx, ...
                      'model', clusterModel, 'metrics', clusterMetrics), ...
    'rowMeta', rowMeta, ...
    'figures', figs ...
);

end

% ============================ Helpers =====================================

function [bestK, bestIdx, bestSil] = autoKmeans(scoreR, rangeK, reps)
bestSil = -Inf; bestK = rangeK(1); bestIdx = [];
for k = rangeK(1):rangeK(2)
    idx = kmeans(scoreR, k, 'Replicates',reps, 'MaxIter',1000, ...
                 'Display','off', 'OnlinePhase','on');
    s = silhouette(scoreR, idx);
    ms = mean(s, 'omitnan');
    if ms > bestSil
        bestSil = ms; bestK = k; bestIdx = idx;
    end
end
end

function [bestK, GMbest, idxBest, crit] = autoGMM(scoreR, rangeK)
% Use BIC to choose K
opts = statset('MaxIter',2000);
BIC = inf(1, diff(rangeK)+1);
GM = cell(size(BIC));
klist = rangeK(1):rangeK(2);
for i = 1:numel(klist)
    k = klist(i);
    GM{i} = fitgmdist(scoreR, k, 'RegularizationValue',1e-6, ...
        'CovarianceType','full', 'SharedCovariance',false, 'Options',opts, 'Replicates',3);
    BIC(i) = GM{i}.BIC;
end
[~, j] = min(BIC);
GMbest = GM{j};
[~, idxBest] = max(posterior(GMbest, scoreR),[],2);
bestK = klist(j);
crit = struct('kList',klist, 'BIC',BIC);
end

function eps = heuristicEps(scoreR, minpts)
% k-distance elbow heuristic (k = minpts)
D = pdist2(scoreR, scoreR, 'euclidean','Smallest',minpts);
% Take the k-th NN distances (row-wise), use 90th percentile as conservative eps
kdist = D(end, :).';
eps = prctile(kdist, 90);
if eps<=0 || ~isfinite(eps)
    eps = median(kdist(kdist>0));
end
end

function figs = makePlots(score, explained, coeff, featCols, idx, numPCs)
figs = struct();
% Scatter in first two PCs
figs.pcScatter = figure('Name','PCA: PC1 vs PC2','Color','w'); %#ok<LFIG>
ax = axes(figs.pcScatter); hold(ax,'on');
if isempty(idx)
    scatter(ax, score(:,1), score(:,2), 18, 'filled', 'MarkerFaceAlpha',0.7);
    title(ax, sprintf('PC1 vs PC2 (no clustering)')); 
else
    cmap = lines(max(2, numel(unique(idx(idx>0)))));
    uc = unique(idx);
    for i = 1:numel(uc)
        m = idx==uc(i);
        if uc(i) == -1 % DBSCAN noise
            scatter(ax, score(m,1), score(m,2), 12, '.', 'MarkerEdgeAlpha',0.4);
        else
            ci = uc(i);
            scatter(ax, score(m,1), score(m,2), 18, cmap(1+mod(ci-1,size(cmap,1)),:), 'filled', 'MarkerFaceAlpha',0.7);
        end
    end
    title(ax, sprintf('PC1 vs PC2 (clusters), PCs kept = %d', numPCs));
end
xlabel(ax, sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(ax, sprintf('PC2 (%.1f%%)', explained(2)));
grid(ax,'on'); box(ax,'on');

% Scree
figs.scree = figure('Name','PCA Scree','Color','w'); %#ok<LFIG>
bar(explained); ylabel('Variance Explained (%)'); xlabel('Principal Component');
yline(100,':'); grid on; title('Scree Plot');

% Loadings (PC1/PC2)
figs.loadings = figure('Name','PC1/PC2 Loadings','Color','w'); %#ok<LFIG>
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
nexttile; bar(coeff(:,1)); set(gca,'XTickLabel',featCols,'XTick',1:numel(featCols)); xtickangle(30);
ylabel('Loading'); title('PC1 Loadings'); grid on;
nexttile; bar(coeff(:,2)); set(gca,'XTickLabel',featCols,'XTick',1:numel(featCols)); xtickangle(30);
ylabel('Loading'); title('PC2 Loadings'); grid on;
end