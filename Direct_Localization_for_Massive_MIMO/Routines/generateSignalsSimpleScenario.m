function [ snapshot, AOAs ] = generateSignalsSimpleScenario( mobileUserLoc, sensorsLocs, receivedPower, carrierFreq, mapWidth, reflectorLoc, BSnlos )

    
% Generate the received signals at multiple widely distributed arrays of sensors.
% The received signals are simulated following the model from (Spencer-Jeffs-Jensen-Swindlehurst'00)
% The emitted signals are narrow-band Gaussian pulses.

% Coded by Nil Garcia

%   Inputs:
%   - mobileUserLoc [m]
%       1x2 matrix
%   - sensorsLocs [m]
%       Cell vector array of length (numBS)
%       The m-th cell contains l (numSensors(m))x2 matrix which stacks the
%       positions of the sensors within such array.
%   - receivedPower
%       (numBS)x1 matrix
%       Fixes the average power of the LOS arrival to each BS.
%   - carrierFreq [Hz]
%       The carrier frequency of the emitted signal.
%   - mapWidth
%       The width and heigh of the map where the source is assumed to be.
%   - reflectorLoc [m]
%       1x2 matrix
%   - BSnlos
%       Indices of the BSs receiving a NLOS path bounced from the
%       reflector.
%
%   Outputs:
%   - receivedSignal
%       Cell array of dimensions (numBS)x1
%       The m-th cell contains the (numSensors(m))x1 samples from the
%       single snapshot.
%   - startIndex
%       Sample index of 'receivedSignal' corresponding to delay 0.
    

% -----------------------------------------------------------
%   Check that the inputs comply with the correct format.
% -----------------------------------------------------------
    
    % Check 'mobileUserLoc'
    if ~ismatrix(mobileUserLoc)
        error('Signal generation error: The mobile user''s location must be represented by a matrix.')
    elseif any(size(mobileUserLoc)~=[1,2])
        error('Signal generation error: The mobile user''s location must be specified by a vector of size 1x2.')
    elseif ~isreal(mobileUserLoc)
        error('Signal generation error: The mobile user''s location coordinates must be real numbers.')
    elseif any(~isfinite(mobileUserLoc))
        error('Signal generation error: The mobile user''s location coordinates include infinite values.')
    end
    
    % Check 'sensorsLocs'
    if ~iscell(sensorsLocs)
        error('Signal generation error: The sensors locations must be represented by cell array.')
    elseif ~isvector(sensorsLocs)
        error('Signal generation error: The sensors locations must be a vector cell array.')
    elseif any(~cellfun(@ismatrix,sensorsLocs))
        error('Signal generation error: The sensors locations for each array must be represented by a matrix.')
    elseif any(cellfun(@(x)size(x,2),sensorsLocs)~=2)
        error('Signal generation error: The number of coordinates of the sensors locations must be 2.')
    elseif any(~cellfun(@isreal,sensorsLocs))
        error('Signal generation error: The sensors locations coordinates must be real numbers.')
    end
    numBS = length(sensorsLocs);
    numSensors = cellfun(@(x)size(x,1),sensorsLocs);
    
    % Check 'carrierFreq'
    if ~isscalar(carrierFreq)
        error('Signal generation error: The carrier frequeny is not a scalar.')
    elseif ~isreal(carrierFreq)
        error('Signal generation error: The carrier frequeny is not a real number.')
    elseif ~isfinite(carrierFreq)
        error('Signal generation error: Infinite is not a valid value for the carrier frequeny.')
    elseif carrierFreq<=0
        error('Signal generation error: The carrier frequeny must be positive.')
    end
    
    % Check 'mapWidth'
    if ~isscalar(mapWidth)
        error('Signal generation error: The width and heigh of the map is not a scalar.')
    elseif ~isreal(mapWidth)
        error('Signal generation error: The width and heigh of the map is not a real number.')
    elseif ~isfinite(mapWidth)
        error('Signal generation error: Infinite is not a valid value for the width and heigh of the map.')
    elseif mapWidth<0
        error('Signal generation error: The width and heigh of the map cannot be negative.')
    end
    
    % Check 'reflectorLoc'
    if ~ismatrix(reflectorLoc)
        error('Signal generation error: The reflectors''s location must be represented by a matrix.')
    elseif any(size(reflectorLoc)~=[1,2])
        error('Signal generation error: The reflectors''s location must be specified by a vector of size 1x2.')
    elseif ~isreal(reflectorLoc)
        error('Signal generation error: The reflectors''s location coordinates must be real numbers.')
    elseif any(~isfinite(reflectorLoc))
        error('Signal generation error: The reflectors''s location coordinates include infinite values.')
    end

    
% -----------------------------------------------------------
%   Signal simulation
% -----------------------------------------------------------

    centerOfGravity = cell2mat(cellfun(@(x)mean(x,1),sensorsLocs,'UniformOutput',false));

    temp = ones(numBS,1)*mobileUserLoc -centerOfGravity;
    anglesLOS = atan(temp(:,2)./temp(:,1)) +(temp(:,1)<0)*pi;
    
    temp = ones(numBS,1)*reflectorLoc -centerOfGravity;
    anglesNLOS = atan(temp(:,2)./temp(:,1)) +(temp(:,1)<0)*pi;   
    
    snapshot = cell(numBS,1);

    AOAs = cell(numBS,1);
    
    % Aggregate LOS path components to signal
    for l=1:numBS

        % Compute Rayleigh-distributed random signal strength
        amplitude = sqrt(receivedPower/2)*sqrt(randn^2+randn^2);

        % Compute random phase
        phase = 2*pi*rand;

        % Aggregate path component to signal
        snapshot{l} = amplitude*exp(1i*phase)*angleToSteering(sensorsLocs{l},anglesLOS(l),carrierFreq);
        
        AOAs{l} = mod(anglesLOS(l),2*pi);
    end
    
    % Aggregate NLOS path components to signal
    for i=1:length(BSnlos)
        
        % Compute Rayleigh-distributed random signal strength
        amplitude = sqrt(receivedPower/2)*sqrt(randn^2+randn^2);

        % Compute random phase
        phase = 2*pi*rand;

        % Aggregate path component to signal
        snapshot{BSnlos(i)} = snapshot{BSnlos(i)} +amplitude*exp(1i*phase)...
            *angleToSteering(sensorsLocs{BSnlos(i)},anglesNLOS(BSnlos(i)),carrierFreq);
        
        AOAs{BSnlos(i)} = [AOAs{BSnlos(i)},mod(anglesNLOS(BSnlos(i)),2*pi)];
    end
end
