function [spotXY, spotAmplitude, fitInfo, keep_all_positive] = refineGaussianGPU(img_stack, spotXY, spotFrame, sigma, win)
% GPU-accelerated 2D Gaussian MLE refinement (amp & bg implicit, x,y,sigma explicit)
% KLS 20250923, bsxfun is likely unnessary since Matlab v2016b as implicit
%   expansion was impliment so I could directly apply the operation instead
%   of also bsxfun to expand the varibles fist

%{
% Future Ideas/Suggestions
2. Implement Per-Spot Convergence
The current termination conditions (if all(...)) check if all fits have 
converged simultaneously. This is inefficient, as the loop continues to 
process spots that have already found a good fit. A much better approach 
is to use a logical mask (active_fits) to track which spots are still 
iterating. Once a spot converges, it's removed from the active set, and 
no further calculations are performed for it. This can lead to substantial 
speedups when fitting a large number of spots with varying difficulty.
%}
    % ---- params ----
    maxIter   = 200;
    lambda0   = 0.01;
    minDelta  = 1e-4;
    minImprove= 1e-6;
    bgMin = gpuArray(single(0.1));  % bkgd can not be less than 0.1

    useFixedSigma = false;
    %singlePrec = true;

    % ---- prep ----
    if isempty(win)
        win = 7;
    end
    half = floor(win/2);
    %[H,W,~] = size(img_stack);
    N  = size(spotXY,1);
    %P  = win*win;

    % pad once
    padImg = padarray(img_stack,[half half],'replicate','both');
    Hp = size(padImg,1); 
    Wp = size(padImg,2);

    % integer center positions in padded coords
    xc = round(spotXY(:,1)) + half;         % integer centers for patch pull
    yc = round(spotXY(:,2)) + half;
    tf = spotFrame(:);

    % offsets for ROI
    [dy,dx] = ndgrid(-half:half,-half:half);

    baseIdx = (tf-1)*Hp*Wp + sub2ind([Hp Wp], yc, xc);

    % move to GPU
    gImg = gpuArray(single(padImg));

    % --- coords & indexing ---
    [Yw, Xw] = ndgrid(1:win, 1:win);        % y first (rows), x second (cols)
    Yv = gpuArray(single(Yw(:)));
    Xv = gpuArray(single(Xw(:)));

    
    offLin  = dy(:) + dx(:)*Hp;             % correct offset
    idxMat  = gpuArray(offLin + baseIdx.'); % keep on GPU
    Z       = gImg(idxMat);                 % P x N

    %keep_all_positive = true([N 1]); % no removeal of localizations near
    %the edge of the mask

    %
    % ---- Cull ROIs that contain 20% or less non-positive pixel (≤ 0) ----------------
    %badROI_gpu = any(Z <= 0, 1);                  % 1×N logical on GPU
    badROI_gpu = mean(Z == 0, 1) > 0.20;   % 1×N logical (gpuArray)
    if any(gather(badROI_gpu))                    % gather once for the if
        keep_all_positive = ~badROI_gpu;                   % for indexing GPU arrays
        %keep_all_positive = ones([size(Z,2) 1]);
        keep     = gather(keep_all_positive);              % for indexing CPU arrays
    
        % Subset GPU data first
        Z      = Z(:, keep_all_positive);
        idxMat = idxMat(:, keep_all_positive);
    
        % Subset CPU-side metadata to stay aligned
        spotXY    = spotXY(keep, :);
        spotFrame = spotFrame(keep);
        xc        = xc(keep);
        yc        = yc(keep);
        baseIdx   = baseIdx(keep);
    
        % Update N
        N = nnz(keep);
    else
        keep_all_positive = ~badROI_gpu;
    end
    %}
    % Early exit if nothing left after culling
    if N == 0
        spotAmplitude = zeros(0,1,'single');
        fitInfo = struct('amp',[],'bg',[],'sigma',[],'chi2',[],'iters',[], ...
                         'lambda',[],'acceptedMask',[],'integrated_signal_roi',[], ...
                         'spotSignal',[],'spotBkgd',[],'spot_NetSignal',[]);
        return;
    end

    clear padImg gImg

    % initial params in window coords
    x0 = gpuArray(single(half+1)) * ones(1,N,'single');
    y0 = gpuArray(single(half+1)) * ones(1,N,'single');
    s0 = gpuArray(single(sigma))   * ones(1,N,'single');

    % --- derivatives ---
    invs2 = 1./(s0.^2);
    dxp   = bsxfun(@minus, Xv, x0);
    dyp   = bsxfun(@minus, Yv, y0);
    r2    = dxp.^2 + dyp.^2;
    f     = exp(-0.5 * bsxfun(@times, r2, invs2));

    dfdx = bsxfun(@times, dxp .* invs2, f);
    dfdy = bsxfun(@times, dyp .* invs2, f);
    if ~useFixedSigma
        dfds = bsxfun(@times, r2 .* (1./(s0.^3)), f);   % df/dσ
    end

    % ---------- LM iterations ----------
    lambda = lambda0*ones(1,N,'single');
    chiPrev = inf(1,N,'single');

    fitChi   = chiPrev;
    fitIter  = zeros(1,N,'uint8');
    everAcc  = false(1,N,'logical');

    % ---------- NEW: bounds for central 3x3 ----------
    c      = gpuArray(single(half+1)); % window center in [1..win]
    xmin_b = c - 1; % lower bound (inclusive)
    xmax_b = c + 1; % upper bound (inclusive)
    ymin_b = xmin_b;
    ymax_b = xmax_b;

    % solve (JTJ*(1+λ) + off-diag) δ = -g
    for it = 1:maxIter
        % --- implicit alpha, beta with guard ---
        P  = size(Z,1);
        F  = sum(f,1);  
        G = sum(Z,1);
        FF = sum(f.^2,1);  
        FG = sum(f.*Z,1);

        den = P*FF - F.^2;
        bad = den <= eps('single')*P.*FF;

        alpha = zeros(1,N,'single'); 
        beta = zeros(1,N,'single');

        good = ~bad;
        % Unconstrained normal-equation solution
        alpha(good) = (P*FG(good) - F(good).*G(good)) ./ den(good);
        beta(good)  = (G(good).*FF(good) - F(good).*FG(good)) ./ den(good);
        % Fallback for ill-conditioned cases
        beta(bad)   = G(bad)/P; % flat bkgd model
        % alpha(bad)=0;
        
        % === enforce β >= bgMin and re-solve α for fixed β ===
        beta        = max(beta, bgMin);
        hasFF       = FF > 0;
        alpha(hasFF)= (FG(hasFF) - beta(hasFF).*F(hasFF)) ./ FF(hasFF);
        alpha(~hasFF)= 0;   % safety

        % (Optional) enforce non-negative peaks
        alpha = max(alpha, 0);

        % --- residual ---
        model = bsxfun(@times,f,alpha) + beta;
        R     = model - Z;
        chi   = sum(R.^2,1);
        
        % --- gradients & JTJ ---
        dfdx = bsxfun(@times, dxp .* invs2, f);
        dfdy = bsxfun(@times, dyp .* invs2, f);
        if ~useFixedSigma
            dfds = bsxfun(@times, r2 .* (1./(s0.^3)), f);
        end
        g1 = alpha .* sum(R.*dfdx,1);
        g2 = alpha .* sum(R.*dfdy,1);
        if ~useFixedSigma, g3 = alpha .* sum(R.*dfds,1); end
        
        JTJ11 = sum((alpha.*dfdx).^2,1);
        JTJ22 = sum((alpha.*dfdy).^2,1);
        JTJ12 = sum((alpha.*dfdx).*(alpha.*dfdy),1);
        if ~useFixedSigma
            JTJ33 = sum((alpha.*dfds).^2,1);
            JTJ13 = sum((alpha.*dfdx).*(alpha.*dfds),1);
            JTJ23 = sum((alpha.*dfdy).*(alpha.*dfds),1);
        end
        
        % --- solve (JTJ + λ·diag) δ = -g ---
        if useFixedSigma
            A11 = JTJ11 + lambda .* JTJ11;
            A22 = JTJ22 + lambda .* JTJ22;
            A12 = JTJ12;
            detA = A11.*A22 - A12.^2;
            dx = ( -g1.*A22 +  g2.*A12) ./ detA;
            dy = ( -g2.*A11 +  g1.*A12) ./ detA;
            ds = zeros(1,N,'single');
        else
            A11 = JTJ11 + lambda .* JTJ11; A22 = JTJ22 + lambda .* JTJ22; A33 = JTJ33 + lambda .* JTJ33;
            A12 = JTJ12; A13 = JTJ13; A23 = JTJ23;
            c11 =  A22.*A33 - A23.*A23;
            c12 = -(A12.*A33 - A13.*A23);
            c13 =  A12.*A23 - A13.*A22;
            c21 =  c12; c22 =  A11.*A33 - A13.*A13;
            c23 = -(A11.*A23 - A13.*A12);
            c31 =  c13; c32 =  c23; c33 =  A11.*A22 - A12.^2;
            detA = A11.*c11 + A12.*c12 + A13.*c13;
            rhs1 = -g1; rhs2 = -g2; rhs3 = -g3;
            dx = (c11.*rhs1 + c12.*rhs2 + c13.*rhs3) ./ detA;
            dy = (c21.*rhs1 + c22.*rhs2 + c23.*rhs3) ./ detA;
            ds = (c31.*rhs1 + c32.*rhs2 + c33.*rhs3) ./ detA;
        end
        
        % --- proposed step and clamps ---
        xNew = x0 + dx;
        yNew = y0 + dy;
        sNew = s0 + ds;

        xNew = min(max(xNew, xmin_b), xmax_b);
        yNew = min(max(yNew, ymin_b), ymax_b);

        sMin = single(0.7);
        sMax = single(max(1.5, 0.45*win));
        sNew = min(max(sNew, sMin), sMax);

        % evaluate proposal (no re-linearization)
        invs2n = 1./(sNew.^2);
        r2n = bsxfun(@minus,Xv,xNew).^2 + bsxfun(@minus,Yv,yNew).^2;
        fn  = exp(-0.5*bsxfun(@times,r2n,invs2n));
        Fn  = sum(fn,1); 
        FFn = sum(fn.^2,1);
        FGn = sum(fn.*Z,1); Gn = G;
        denN= size(Z,1)*FFn - Fn.^2;
        alphaN = (size(Z,1)*FGn - Fn.*Gn)./denN;
        betaN  = (Gn.*FFn - Fn.*FGn)./denN;
        Rn     = bsxfun(@times,fn,alphaN) + betaN - Z;
        chiNew = sum(Rn.^2,1);

        accept = chiNew < chiPrev;
        lambda = min(lambda .* (accept*0.3 + (~accept)*4.0), 1e4);

        P    = size(Z,1);
        denN = P*FFn - Fn.^2;
        badN = denN <= eps('single')*P.*FFn;
        
        alphaN = zeros(1,N,'single');
        betaN  = zeros(1,N,'single');
        
        goodN = ~badN;
        % Unconstrained
        alphaN(goodN) = (P*FGn(goodN) - Fn(goodN).*Gn(goodN)) ./ denN(goodN);
        betaN(goodN)  = (Gn(goodN).*FFn(goodN) - Fn(goodN).*FGn(goodN)) ./ denN(goodN);
        % Fallback
        betaN(badN)   = Gn(badN)/P;
        
        % === enforce β >= bgMin and re-solve α for fixed β ===
        betaN          = max(betaN, bgMin);
        hasFFn         = FFn > 0;
        alphaN(hasFFn) = (FGn(hasFFn) - betaN(hasFFn).*Fn(hasFFn)) ./ FFn(hasFFn);
        alphaN(~hasFFn)= 0;


        % evaluate proposal (no re-linearization)
        invs2n = 1./(sNew.^2);
        r2n = bsxfun(@minus,Xv,xNew).^2 + bsxfun(@minus,Yv,yNew).^2;
        fn  = exp(-0.5*bsxfun(@times,r2n,invs2n));
        Fn  = sum(fn,1); 
        FFn = sum(fn.^2,1);
        FGn = sum(fn.*Z,1); Gn = G;
        denN= size(Z,1)*FFn - Fn.^2;
        alphaN = (size(Z,1)*FGn - Fn.*Gn)./denN;
        betaN  = (Gn.*FFn - Fn.*FGn)./denN;
        Rn     = bsxfun(@times,fn,alphaN) + betaN - Z;
        chiNew = sum(Rn.^2,1);

        improve = (chiPrev - chiNew)./chiPrev; % absolute value of the improvement?
        accept = chiNew < chiPrev;
        lambda = min(lambda .* (accept*0.3 + (~accept)*4.0), 1e4);

        % LM damping
        lambda = lambda .* (~accept*10 + accept*0.1);
        lambda = min(lambda,1e4);

        % commit accepted steps (already clamped)
        x0(accept) = xNew(accept);
        y0(accept) = yNew(accept);
        if ~useFixedSigma
            s0(accept) = sNew(accept);
        end

        alpha(accept) = alphaN(accept);
        beta(accept)  = betaN(accept);
        chiPrev(accept)= chiNew(accept);

        fitChi(accept)  = chiNew(accept);
        fitIter(accept) = uint8(it);
        everAcc         = everAcc | accept;

        if all((abs(improve) < minImprove & it > 25) | improve < 0)
            disp(['improvement (' num2str(max(improve(:))) ') less than minImprove (' num2str(minImprove) '), itration ' num2str(it)])
            break; 
        end
        if all(abs([dx dy ds])./max([abs(x0) abs(y0) abs(s0)],1e-6) < minDelta,'all')
            disp('all fits less than minDelta')
            break; 
        end
    end
    
    % ----- integrated areas (discrete, over the ROI window) -----
    invs2 = 1./(s0.^2);
    r2    = bsxfun(@minus,Xv,x0).^2 + bsxfun(@minus,Yv,y0).^2;
    f     = exp(-0.5*bsxfun(@times,r2,invs2));
    Fsum  = sum(f,1);
    %P     = size(Z,1);

    areaSignal_patch        = alpha .* Fsum;                 % α·Σf (ROI discrete)
    %areaBkgd_patch_true     = beta  .* P;                    % β·P   (ROI discrete)
    %areaNet_patch_model     = areaSignal_patch;              % model net in ROI
    %areaNet_patch_data      = sum(Z,1) - areaBkgd_patch_true;% data net in ROI
    Aeq = 2*pi*(s0.^2);                            % PSF-equivalent area
    spotSignal_cont         = alpha .* Aeq;       % α·2πσ² (continuous)
    spotBkgd_cont           = beta .* Aeq;       % ß·2πσ² (continuous)
    spotNet_cont            = spotSignal_cont - spotBkgd_cont; 


    % back to global coords
    spotXY(:,1) = gather(x0.' + (spotXY(:,1)- (half+1)));
    spotXY(:,2) = gather(y0.' + (spotXY(:,2)- (half+1)));
    spotAmplitude    = gather(alpha.');

    fitInfo = struct( ...
        'amp',    gather(alpha.'), ...
        'bg',     gather(beta.'), ...
        'sigma',  gather(s0.'), ...
        'chi2',   gather(fitChi.'), ...
        'iters',  gather(fitIter.'), ...
        'lambda', gather(lambda.'), ...
        'acceptedMask', gather(everAcc.'), ...
        'integrated_signal_roi',       gather(areaSignal_patch.'), ...
        'spotSignal',   gather(spotSignal_cont.'), ...
        'spotBkgd',    gather(spotBkgd_cont.'), ...
        'spot_NetSignal',     gather(spotNet_cont.'));

    %{
    % Pause here if you want to visualize what is going on with the raw
    % data and it's fits?
    %%
    % 2D data review
    k = 6704;% pick a spot
    
    % --- gather needed vars
    Zk     = gather(Z(:,k));
    x0k    = gather(x0(k));
    y0k    = gather(y0(k));
    sigk   = gather(s0(k));
    alphak = gather(alpha(k));
    betak  = gather(beta(k));
    
    % coords
    [Yw,Xw] = ndgrid(1:win, 1:win);
    Xv = Xw(:); 
    Yv = Yw(:);
    
    % model
    fk   = exp(-0.5 * ((Xv-x0k).^2 + (Yv-y0k).^2) / (sigk^2));
    Zfit = alphak*fk + betak;
    
    % reshape for display
    patch   = reshape(Zk,   win, win); % Image data
    fitImg  = reshape(Zfit, win, win); % Fit of that data
    resid   = patch - fitImg; % Difference between the data and it's fit
    
    figure;
    subplot(1,3,1); imagesc(patch); 
        axis image; title('Data'); colorbar; colormap('gray')
    subplot(1,3,2); imagesc(fitImg); 
        axis image; title('Fit');  colorbar
    subplot(1,3,3); imagesc(resid);  
        axis image; title('Residual'); colorbar
    hold on; 
        plot(x0k, y0k, 'r+', 'MarkerSize',10,'LineWidth',1.5);
    hold off;

    %%
    % ===== 3-D review: data vs fit (αf + β) with background plane =====
    k = 1;

    % --- gather needed vars (CPU scalars for plotting) ---
    Zk     = gather(Z(:,k));
    x0k    = gather(x0(k));
    y0k    = gather(y0(k));
    sigk   = gather(s0(k));
    alphak = gather(alpha(k));
    betak  = gather(beta(k));
    chi2k  = gather(fitChi(k));
    
    % --- coords (consistent with main code: [Yw,Xw]) ---
    [Yw, Xw] = ndgrid(1:win, 1:win);
    X = Xw;  Y = Yw;
    
    % --- reshape data and build model surfaces ---
    patch   = reshape(Zk,   win, win);                     % data
    fk      = exp(-0.5 * ((X-x0k).^2 + (Y-y0k).^2) / (sigk^2));
    fitSurf = alphak .* fk + betak;                        % αf + β (model)
    bgPlane = betak  .* ones(win, win, 'like', fitSurf);   % β
    zPeak   = alphak + betak;                              % model peak height
    
    % --- optional diagnostics ---
    Fsum    = sum(fk(:));                  % Σ f in ROI
    Aeq     = 2*pi*sigk^2;                 % PSF-equivalent area
    capture = Fsum / Aeq;                  % fraction of flux captured by ROI
    
    % --- plot ---
    figure('Name','Gaussian Fit (3D)','Color','w');
    %s1 = surf(X, Y, patch,  'EdgeColor','none','FaceAlpha',0.60); hold on;
    s1 = plot3(X(:), Y(:), patch(:), 'k.', 'MarkerSize',12, 'LineWidth',1.5); hold on;
    s2 = surf(X, Y, fitSurf,'EdgeColor','none','FaceAlpha',0.25);
    s3 = surf(X, Y, bgPlane,'EdgeColor','none','FaceAlpha',0.25,'FaceColor',[0.6 0.6 0.6]);
    
    plot3(x0k, y0k, zPeak, 'r+', 'MarkerSize',12, 'LineWidth',1.5);
    
    % styling
    colormap(turbo);
    shading interp; 
    view(45, 35); axis tight; 
    axis vis3d; 
    box off; 
    grid on;
    xlabel('x (px)'); 
    ylabel('y (px)'); 
    zlabel('Intensity (units provide by input image)');
    
    legend([s1 s2 s3], {'Data','Fit (\alpha f + \beta)','Background \beta'}, ...
           'Location','northeast');
    
    title(sprintf('Spot k=%d | \\alpha=%.3g, \\beta=%.3g, \\sigma=%.3g px, \\chi^2=%.3g | capture=%.2f%%', ...
          k, alphak, betak, sigk, chi2k, 100*capture));
    %}
end