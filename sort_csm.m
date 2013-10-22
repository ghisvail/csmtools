function [chan_raw_sorted, body_raw_sorted] = sort_csm(chan_lab, chan_raw, ...
        body_lab, body_raw)
%
%
% Ghislain Vaillant <ghislain.vaillant@kcl.ac.uk>

    %% extract geometry information
    xDim = size(chan_raw, 1);
    yMax = max(chan_lab(:, 1), [], 1);
    yMin = min(chan_lab(:, 1), [], 1);
    yDim = 1 + yMax - yMin;
    zMax = max(chan_lab(:, 2), [], 1);
    zMin = min(chan_lab(:, 2), [], 1);
    zDim = 1 + zMax - zMin;
    
    if yMin > 0
        warning('offset in y detected');
    end

    if zMin > 0
        warning('offset in z detected');
    end    
    
    nCoils = size(chan_raw, 2);
    nsa = nnz(chan_lab(:, 1) == yMax) / zDim; 
    nPE = yDim * zDim;
    
    % re-order raw data matrices
    chan_raw = reshape(chan_raw,  xDim, nCoils, yDim*zDim, nsa);
    chan_raw_sorted = zeros(xDim, nCoils, yDim, zDim, nsa); 
    for iPE = 1:nPE
        iY = 1 + chan_lab(iPE, 1);
        iZ = 1 + chan_lab(iPE, 2);
        chan_raw_sorted(:, :, iY, iZ, :) = chan_raw(:, :, iPE, :);
    end
    
    body_raw = reshape(body_raw,  xDim, 1, yDim*zDim, nsa);
    body_raw_sorted = zeros(xDim, 1, yDim, zDim, nsa);
    for iPE = 1:nPE
        iY = 1 + body_lab(iPE, 1);
        iZ = 1 + body_lab(iPE, 2);
        body_raw_sorted(:, 1, iY, iZ, :) = body_raw(:, 1, iPE, :);
    end    
    
    % make coil dimension last
    chan_raw_sorted = permute(chan_raw_sorted, [1 3 4 5 2]);
    body_raw_sorted = permute(body_raw_sorted, [1 3 4 5 2]);
