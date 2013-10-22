function [ csm_chan, csm_body ] = recon_csm( chan_raw, ...
        body_raw, acq_res, rec_res, rec_matrix, varargin )
%
%
% Ghislain Vaillant <ghislain.vaillant@kcl.ac.uk>

    threshold = [];
    
    if ~isempty(varargin)
        for iArg=1:length(varargin)
            if isscalar(varargin{iArg}) || ischar(varargin{iArg})
                switch varargin{iArg}                   
                    case {'threshold'}
                        % specify the number of coils to read
                        threshold = varargin{iArg+1};
                end
            end
        end
    end
    
    % default threshold value
    if isempty(threshold)
        threshold = 0.05;
    end

    if any([isempty(acq_res), isempty(rec_res)])
        error('ValueError: all geometry informations must be provided')
    end

   %% polynomial fit
    fprintf('performing fit...\n')
    
    % malloc
    xDim = size(body_raw, 1);
    yDim = size(body_raw, 2);
    zDim = size(body_raw, 3);
    nsa = size(chan_raw, 4);
    nCoils = size(chan_raw, 5);
    
    chan_all = zeros(xDim, yDim, zDim, nCoils);
    body_all = zeros(xDim, yDim, zDim);
    
    for iz = 1:zDim
        % debug
        fprintf('fitting slice %d of %d\n', iz, zDim)
        
        % malloc
        chan_z = zeros(xDim, yDim, nsa, nCoils); 
        body_z = zeros(xDim, yDim, nsa);
        
        % reconstruct current slice
        for isa = 1:nsa       
             % multi-channel data
            for iCoil = 1:nCoils
               chan_z(:, :, isa, iCoil) = (...
                    ktoi(chan_raw(:, :, iz, isa, iCoil)));
               chan_z(:, :, isa, iCoil) = (...
                    ifftshift(chan_z(:, :, isa, iCoil), 2)); 
            end
            % body coil data
            body_z(:, :, isa) = (...
                ktoi(body_raw(:, :, iz, isa)));
            body_z(:, :, isa) = (...
                ifftshift(body_z(:, :, isa), 2));
        end
        
        % mean of NSA images
        chan_z = squeeze(mean(chan_z, 3));
        body_z = squeeze(mean(body_z, 3));
        
        % calculate sensitivity maps
        [s, rms] = coil_smap(chan_z, body_z);
        
        % create mask for low values
        mean_body_z = abs(body_z);
        max_body_z = max(mean_body_z(:));
        mask_threshold = max_body_z * threshold;
        mask = (mean_body_z > mask_threshold);
        
        % post-process mask to remove holes
        dilation_kernel = strel('disk', 6);
        mask = imdilate(mask, dilation_kernel);
        erode_kernel = strel('disk', 3);
        mask = imerode(mask, erode_kernel);
        
        % perform fitting
        fitS_real = zeros(xDim, yDim);
        fitS_imag = zeros(xDim, yDim);
        [yTraj, xTraj] = meshgrid(1:yDim, 1:xDim);
        [yTrajFit, xTrajFit] = meshgrid(1:yDim, 1:xDim);
        masked_out = find(~mask);
        xTraj(masked_out(end:-1:1)) = [];
        yTraj(masked_out(end:-1:1)) = [];
        for iCoil = 1:nCoils
            % discard data outside mask
            s_coil = s(:, :, iCoil);
            s_coil = s_coil(:);
            s_coil(masked_out(end:-1:1)) = [];
            % fit real part
            poly = polyfitn([xTraj(:), yTraj(:)], real(s_coil(:)), 7);
            fitS_real_tmp = polyvaln(poly, [xTrajFit(:), yTrajFit(:)]);
            fitS_real_tmp(~mask) = 0;
            fitS_real(:, :, iCoil) = (...
                    reshape(fitS_real_tmp, xDim, yDim));
            % fit imaginary part
            poly = polyfitn([xTraj(:), yTraj(:)], imag(s_coil(:)), 7);
            fitS_imag_tmp = polyvaln(poly, [xTrajFit(:), yTrajFit(:)]);
            fitS_imag_tmp(~mask) = 0;
            fitS_imag(:,:,iCoil) = (...
                    reshape(fitS_imag_tmp, xDim, yDim));            
        end
        fitS = fitS_real + 1i * fitS_imag; 
        
        % store and loop
        chan_all(:, :, iz, :) = fitS;
        body_all(:, :, iz) = body_z;
    end


    %% interpolation to reconstructed resolution
    fprintf('perform interpolation...\n')
    
    csm_chan = chan_all;
    csm_body = body_all;    
    nDims = length(acq_res);
    
    for iDim = 1:nDims
        old_res = acq_res(iDim);
        new_res = rec_res(iDim);       
        if new_res ~= old_res
            ratio_res = new_res / old_res;
            if ratio_res < 1
                old_shape = size(csm_body);
                new_shape = [ceil(old_shape(1)/ratio_res) old_shape(2:end)];
 
                ll = 1 + floor(new_shape(1) / 2 .* (1 - ratio_res));
                rr = ll + old_shape(1) - 1;
                                
                zpad_body = zeros(new_shape);
                zpad_body(ll:rr, :, :) = itok(csm_body, 1);
                csm_body = ktoi(zpad_body, 1);
                
                zpad_chan = zeros([new_shape, nCoils]);
                for iCoil = 1:nCoils
                    zpad_chan(ll:rr, :, :, iCoil) = (...
                        itok(csm_chan(:, :, :, iCoil), 1));    
                end
                csm_chan = ktoi(zpad_chan, 1);
            else
                old_shape = size(csm_body);
                new_shape = [ceil(old_shape(1)/ratio_res) old_shape(2:end)];
                
                ll = 1 + floor(old_shape(1) / 2 .* (1 - 1/ratio_res));
                rr = ll + floor(new_shape(1)) - 1;                
                
                csm_body = itok(csm_body(ll:rr, :, :), 1);
                csm_body = ktoi(csm_body, 1);               
                
                csm_chan = itok(csm_chan(ll:rr, :, :, :), 1);
                csm_chan = ktoi(csm_chan, 1);             
            end
        end
 
        csm_body = permute(csm_body, [2 3 1]);
        csm_chan = permute(csm_chan, [2 3 1 4]);       
    end

    
    %if ~all(acq_res == rec_res)
         %inc_res = acq_res ./ rec_res;
         %old_shape = size(body_all);
         %new_shape = floor(old_shape .* inc_res);
         
         %% malloc
         %csm_chan = zeros([new_shape nCoils]);
         %csm_body = zeros(new_shape);
         
         %% perform Fourier interpolation on body coil data
         %ll = 1 + floor(new_shape / 2 .* (1 - 1 ./ inc_res));
         %rr = ll + old_shape - 1;
         %csm_body(ll(1):rr(1), ll(2):rr(2), ll(3):rr(3)) = itok(body_all);
         %csm_body = ktoi(csm_body);
         
         %% perform Fourier interpolation on each channel
         %for iCoil = 1:nCoils
             %csm_chan(ll(1):rr(1), ll(2):rr(2), ll(3):rr(3), iCoil) = (...
                     %itok(chan_all(:, :, :, iCoil)));
             %csm_chan(:, :, :, iCoil) = (...
                     %ktoi(csm_chan(:, :, :, iCoil)));             
         %end
    %else
         %csm_body = body_all;
         %csm_chan = chan_all;
    %end


    %% extract reconstructed FOV
    fprintf('extract ROI...\n')
    
    if ~isempty(rec_matrix)
        old_shape = size(csm_body);
        new_shape = rec_matrix;
        inc_matrix = old_shape ./ new_shape;
        % calculate indexes
        ll = 1 + floor(old_shape / 2 .* (1 - 1 ./ inc_matrix));
        rr = ll + new_shape - 1;       
        % extract desired region
        csm_body = csm_body(ll(1):rr(1), ll(2):rr(2), ll(3):rr(3));
        csm_chan = csm_chan(ll(1):rr(1), ll(2):rr(2), ll(3):rr(3), :);       
    end
end
