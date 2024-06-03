function [ data] = convertRawToFrame( data, Acqp, varargin)
% function frame = convertRawToFrame(data, Acqp, ['specified_NRs', NRarray])
% 
% Input:
%   data: the raw data (fid/rawdata.job*), as read by ReadBrukerRaw 
%
%   Acqp: An acqp struct as generated by the function readBrukerParamFile('path/acqp')
%
% Optional Input:   
%   'specified_NRs', NRarray: A list of NRs to be read, NR starting with 1 
%                             'specified_NRs',[2 5 7] -> only NR 2, 5 and 7 are converted
%
% Output:
%   frame: sorted 5D-Matrix with dimensions 
%          (Scansize, ScansPerFrame, NumberOfReceiveChannels, NumberOfObjects, NumberOfRepetitions)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2013
% Bruker BioSpin MRI GmbH
% D-76275 Ettlingen, Germany
%
% All Rights Reserved
%
% $Id: convertRawToFrame.m,v 1.3 2013/07/05 12:56:02 haas Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% Define default-value if necesseary:
    
    % Check arguments
    [trash,specified_NRs]=bruker_addParamValue(varargin, 'specified_NRs','@(x) isnumeric(x)',[]);
    clear trash;

    % Check data input
    if iscell(data)
        if length(data) == 1
            data = data{1};
        else
            error('Raw data is a cell array with multiple elements. Specify which job you would like to convert, e.g. data{1}');
        end
    end
    
    % Check for missing variables in structs:
    cellstruct{1}=Acqp;
    all_here = bruker_requires(cellstruct, {{'Acqp','NI','NR','ACQ_size','ACQ_phase_factor','ACQ_obj_order', 'ACQ_dim'}});
    clear cellstruct;
    if ~all_here
        error('Some parameters are missing');
    end
    
    %Init:
    complexfid=~isreal(data);
    NI=Acqp.NI;
    ACQ_size=Acqp.ACQ_size;
    ACQ_phase_factor=Acqp.ACQ_phase_factor;
    ACQ_obj_order=Acqp.ACQ_obj_order;
    ACQ_dim=Acqp.ACQ_dim;
    if ~isempty(specified_NRs)
        NR=length(specified_NRs);
    else
        NR=Acqp.NR; 
    end
    
    % read precision:
    temp=whos('data');
    precision=temp.class;
    clear temp;
    % Convert precision-string to boolean-variable:
    if(strcmpi(precision, 'single'))
        memsave=true;
    else
        memsave=false;
    end
    clear precision;   

%-----------------------------------------------
    %% Calculate additional Parameters
    
    % calculate numSelectedReceivers
    numSelectedReceivers=size(data,1);
    
    % Calculating number of elements in higher dimensions
    numDataHighDim=prod(ACQ_size(2:end));    

    % Convert if complex: to blockSize of a complex Matrix and change
    % ACQ_size(1)
    if complexfid
        scanSize(1)=ACQ_size(1)/2;
    else
        scanSize(1)=ACQ_size(1);
    end  

    %% Start resort:
    
    
        
    if ACQ_dim>1

        % Permuting receiver and read direction
        % data entering: (scanSize, numDataHighDim*NI*NR, numSelectedReceivers)
        
        data=reshape(data, numSelectedReceivers, scanSize, ACQ_phase_factor, NI, numDataHighDim/ACQ_phase_factor,NR);
        data=permute(data, [2 3 5 1 4 6]); % => scansize, ACQ_phase_factor, numDataHighDim/ACQ_phase_factor, numSelectedReceivers, NI, NR
        data=reshape(data, scanSize, numDataHighDim, numSelectedReceivers, NI, NR);

        data(:,:,:,ACQ_obj_order+1,:)=data;
    else
        % dim=1:
        data=reshape(data,numSelectedReceivers, scanSize,1,NI,NR);
        data=permute(data, [2 3 1 4 5]);
    end
end
