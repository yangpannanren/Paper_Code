function [ estimAOA ] = beamforming( sensorsLocs, snapshot, carrierFreq)

% Coded by Nil Garcia

%   Inputs:
%   - sensorsLocs [m]
%       Cell vector array of length (numBS)
%       The m-th cell contains l (numSensors(m))x2 matrix which stacks the
%       positions of the sensors within such array.
%   - snapshot
%       Cell array of dimensions (numBS)x1
%       The m-th cell contains the (numSensors(m))x1 samples from the
%       single snapshot.
%   - carrierFreq [Hz]
%       The carrier frequency of the emitted signal.
%
%   Outputs:
%   - estimAOA
%       Vector of estimated AOA's.
    

% -----------------------------------------------------------
%   Check that the inputs comply with the correct format.
% -----------------------------------------------------------

    % Check 'sensorsLocs'
    if ~iscell(sensorsLocs)
        error('Beamforming error: The sensors locations must be represented by cell array.')
    elseif ~isvector(sensorsLocs)
        error('Beamforming error: The sensors locations must be a vector cell array.')
    elseif any(~cellfun(@ismatrix,sensorsLocs))
        error('Beamforming error: The sensors locations for each array must be represented by a matrix.')
    elseif any(cellfun(@(x)size(x,2),sensorsLocs)~=2)
        error('Beamforming error: The number of coordinates of the sensors locations must be 2.')
    elseif any(~cellfun(@isreal,sensorsLocs))
        error('Beamforming error: The sensors locations coordinates must be real numbers.')
    end
    numBS = length(sensorsLocs);
    numSensors = cellfun(@(x)size(x,1),sensorsLocs);
    
    % Check 'snapshot'
    if ~iscell(snapshot)
        error('Beamforming error: The spatial samples must be represented by cell array.')
    elseif ~isvector(snapshot)
        error('Beamforming error: The spatial samples must be a vector cell array.')
    elseif length(snapshot)~=numBS
        error(['Beamforming error: The spatial samples must be defined for ', num2str(numBS), ' base stations.'])
    elseif any(~cellfun(@isvector,snapshot))
        error('Beamforming error: The spatial samples for each array must be represented by a vector.')
    elseif any(cellfun(@(x,y)length(x)~=y,snapshot,num2cell(numSensors)))
        error('Beamforming error: The spatial samples are not defined for all antennas.')
    elseif any(~cellfun(@isnumeric,snapshot))
        error('Beamforming error: The spatial samples must be numbers.')
    end
    
    % Check 'carrierFreq'
    if ~isscalar(carrierFreq)
        error('Beamforming error: The carrier frequeny is not a scalar.')
    elseif ~isreal(carrierFreq)
        error('Beamforming error: The carrier frequeny is not a real number.')
    elseif ~isfinite(carrierFreq)
        error('Beamforming error: Infinite is not a valid value for the carrier frequeny.')
    elseif carrierFreq<=0
        error('Beamforming error: The carrier frequeny must be positive.')
    end
    
    
% -----------------------------------------------------------
%   Beamforming
% -----------------------------------------------------------

    angles = 2*pi/1e3:2*pi/1e3:2*pi;
    estimAOA = nan(numBS,1);
    
    for l=1:numBS
    
        [~,ind] = max(abs(angleToSteering(sensorsLocs{l},angles,carrierFreq)'*snapshot{l}));

%         % TEST
%         figure
%         semilogy(180/pi*angles,abs(angleToSteering(sensorsLocs{l},angles,carrierFreq)'*snapshot{l}))
%         figure
%         plot(180/pi*angles,abs(angleToSteering(sensorsLocs{l},angles,carrierFreq)'*snapshot{l}))
%         % TEST

        estimAOA(l) = angles(ind);
    
    end

end

