function [ chan_lab, chan_raw, body_lab, body_raw ] = read_csm(...
        chan_file, body_file)
%
%
% Ghislain Vaillant <ghislain.vaillant@kcl.ac.uk>

    file_ = body_file;
    if exist(file_, 'file') ~= 2
        error('ValueError: file %s not found', file_)
    end

    file_ = chan_file;
    if exist(file_, 'file') ~= 2
        error('ValueError: file %s not found', file_)
    end

    chan_lab_file = strcat(chan_file(1:end-4), '.lab');

    if exist(chan_lab_file, 'file') ~= 2
        error('channel scan label file not found at %s', chan_lab_file)
    end

    body_lab_file = strcat(body_file(1:end-4), '.lab');

    if exist(body_lab_file, 'file') ~= 2
        error('channel scan label file not found at %s', body_lab_file)
    end    
    
    for iCoil = 1:32
        [chan_lab, chan_raw(:, iCoil, :)] = (...
                philipsmr_get_raw_data(chan_lab_file, chan_file, ...
                'no_bar', 'verbose', 'coil', iCoil));   
        
    end
    
    [body_lab, body_raw] = (...
            philipsmr_get_raw_data(body_lab_file, body_file, ...
            'no_bar', 'verbose', 'coil', 1));   
