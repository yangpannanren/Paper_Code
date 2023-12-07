function [ TOAs, AOAs, compAmpli ] = generateChannels( mobileUserLoc, sensorsLocs, receivedPowers, clusterDecayTime, rayDecayTime, clusterInterTime, rayInterTime, angleStdDev )

    
% Generate the multipath channel at all base stations according to the model from (Spencer-Jeffs-Jensen-Swindlehurst'00)
% The emitted signals are narrow-band Gaussian pulses.

% Coded by Nil Garcia

%   Inputs:
%   - mobileUserLoc [m]
%       1x2 matrix
%   - sensorsLocs [m]
%       Cell vector array of length (numBS)
%       The m-th cell contains l (numSensors(m))x2 matrix which stacks the
%       positions of the sensors within such array.
%   - receivedPowers
%       (numBS)x1 matrix
%       Fixes the average power of the LOS arrival to each BS.
%   - clusterDecayTime
%       Decides how quickly decreases the power of each cluster of arrivals.
%   - rayDecayTime
%       Decides how quickly decreases the power of each arrival within a cluster.
%   - clusterInterTime
%       Mean time difference between two clusters.
%   - rayInterTime
%       Mean time difference between two arrivals.
%   - angleStdDev
%       Standard deviation of the Laplacian random variable controlling the
%       angles of arrival within a cluster.
%
%   Outputs:
%   - TOAs
%       Cell vector array of length (numBS)
%       The m-th cell contains a vector of all TOAs at the m-th BS.
%   - AOAs
%       Cell vector array of length (numBS)
%       The m-th cell contains a vector of all AOAs at the m-th BS.
%   - compAmpli
%       Cell vector array of length (numBS)
%       The m-th cell contains a vector of all complex amplitudes at the m-th BS.
    

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
    
    % Check 'clusterDecayTime'
    if ~isscalar(clusterDecayTime)
        error('Signal generation error: The cluster decay time is not a scalar.')
    elseif ~isreal(clusterDecayTime)
        error('Signal generation error: The cluster decay time is not a real number.')
    elseif ~isfinite(clusterDecayTime)
        error('Signal generation error: Infinite is not a valid value for the cluster decay time.')
    elseif clusterDecayTime<=0
        error('Signal generation error: The cluster decay time must be positive.')
    end
    
    % Check 'rayDecayTime'
    if ~isscalar(rayDecayTime)
        error('Signal generation error: The ray decay time is not a scalar.')
    elseif ~isreal(rayDecayTime)
        error('Signal generation error: The ray decay time is not a real number.')
    elseif ~isfinite(rayDecayTime)
        error('Signal generation error: Infinite is not a valid value for the ray decay time.')
    elseif rayDecayTime<=0
        error('Signal generation error: The ray decay time must be positive.')
    end
    
    % Check 'clusterInterTime'
    if ~isscalar(clusterInterTime)
        error('Signal generation error: The mean cluster inter arrival time is not a scalar.')
    elseif ~isreal(clusterInterTime)
        error('Signal generation error: The mean cluster inter arrival time is not a real number.')
    elseif ~isfinite(clusterInterTime)
        error('Signal generation error: Infinite is not a valid value for the mean cluster inter arrival time.')
    elseif clusterInterTime<=0
        error('Signal generation error: The mean cluster inter arrival time must be positive.')
    end
    
    % Check 'rayInterTime'
    if ~isscalar(rayInterTime)
        error('Signal generation error: The mean ray inter arrival time is not a scalar.')
    elseif ~isreal(rayInterTime)
        error('Signal generation error: The mean ray inter arrival time is not a real number.')
    elseif ~isfinite(rayInterTime)
        error('Signal generation error: Infinite is not a valid value for the mean ray inter arrival time.')
    elseif rayInterTime<=0
        error('Signal generation error: The mean ray inter arrival time must be positive.')
    end
    
    % Check 'angleStdDev'
    if ~isscalar(angleStdDev)
        error('Signal generation error: The angle standard deviation is not a scalar.')
    elseif ~isreal(angleStdDev)
        error('Signal generation error: The angle standard deviation is not a real number.')
    elseif ~isfinite(angleStdDev)
        error('Signal generation error: Infinite is not a valid value for the angle standard deviation.')
    elseif angleStdDev<=0
        error('Signal generation error: The angle standard deviation must be positive.')
    end

    
% -----------------------------------------------------------
%   Signal simulation
% -----------------------------------------------------------

    centerOfGravity = cell2mat(cellfun(@(x)mean(x,1),sensorsLocs,'UniformOutput',false));
    temp = ones(numBS,1)*mobileUserLoc -centerOfGravity;
    anglesLOS = atan(temp(:,2)./temp(:,1)) +(temp(:,1)<0)*pi;
    rangeLOS = sqrt(sum(temp.^2,2));
    
    TOAs = cell(numBS,1);
    [TOAs{:}] = deal(nan(1e3,1));
    AOAs = cell(numBS,1);
    [AOAs{:}] = deal(nan(1e3,1));
    compAmpli = cell(numBS,1);
    [compAmpli{:}] = deal(nan(1e3,1));
    
    for l=1:numBS

        % Compute parameters first cluster
        clusterTime = 0;
        clusterAngle =  anglesLOS(l);

        TOAcount = 0;
        while exp(-clusterTime/clusterDecayTime)>1e-3 % Loops over the clusters

            % Compute parameters first ray in the cluster
            rayTime = 0;
            if clusterTime == 0
                % If it's l LOS path, then the angle is fixed.
                rayAngle =  0;
            else
                rayAngle =  angleStdDev*laprnd(1);
            end

            while exp(-rayTime/rayDecayTime)>1e-3 % Loops over the rays                
                % Compute average power
                powerArrival = receivedPowers(l)*exp(-clusterTime/clusterDecayTime-rayTime/rayDecayTime);

                % Compute Rayleigh-distributed random signal strength
                % amplitude = sqrt(powerArrival/2)*sqrt(randn^2+randn^2);
                amplitude = sqrt(powerArrival/2)*sqrt(randn^2+randn^2);
                
                % Compute random phase
                phase = 2*pi*rand;
                
                % Store TOA time
                TOAs{l}(TOAcount+1) = rangeLOS(l)/3e8+clusterTime+rayTime;
                AOAs{l}(TOAcount+1) = mod(clusterAngle+rayAngle,2*pi);
                compAmpli{l}(TOAcount+1) = amplitude*exp(1i*phase);
                TOAcount = TOAcount +1;
                
                % ----------------------------------------
                % Compute random ray time for next arrival
                rayTime = rayTime +exprnd(rayInterTime);
                
                rayAngle =  angleStdDev*laprnd(1);
            end

            % --------------------------------------------
            % Compute random cluster time for next arrival
            clusterTime = clusterTime +exprnd(clusterInterTime);

            % Compute random cluster angle
            clusterAngle =  2*pi*rand;
        end
        
        TOAs{l} = TOAs{l}(1:TOAcount);
        AOAs{l} = AOAs{l}(1:TOAcount);
        compAmpli{l} = compAmpli{l}(1:TOAcount);
    end
    
end
