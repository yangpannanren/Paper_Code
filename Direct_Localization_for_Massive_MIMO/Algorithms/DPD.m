function [ mobileUserLocEstimate ] = DPD( sensorsLocs, receivedSignals, carrierFreq, bandwidth, oversamplingFactor, diamaterSearchArea, resolution, numIter )

% Direct Positiong Determination (DPD) from paper
% <A. J. Weiss and A. Amar, "Direct position determination of multiple radio signals," EURASIP Journal on Applied Signal Processing, 2005.>.
% Coded by Nil Garcia

%   Inputs:
%   - sensorsLocs [m]
%       Cell vector array of length (numBS)
%       The m-th cell contains l (numSensors(m))x2 matrix which stacks the
%       positions of the sensors within such array.
%   - receivedSignals
%       Cell array of dimensions (numBS)x1
%       The m-th cell contains a matrix of size (numSensors(m))x(numSamples)
%   - carrierFreq [Hz]
%       The carrier frequency.
%   - bandwidth [Hz]
%       The time-domain sampling frequency.
%   - oversamplingFactor
%       The oversampling factor at reception
%   - resolution [m]
%       The starting resolution of grid employed by DPD.
%   - diamaterSearchArea [m]
%   - numIter
%       # of grid refinements
%
%   Outputs:
%   - mobileUserLocEstimate [m]
%       Matrix of size 1x2 with the coordinates of the location of the M
    

% -----------------------------------------------------------
%   Check that the inputs comply with the correct format.
% -----------------------------------------------------------

    % Check 'sensorsLocs'
    if ~iscell(sensorsLocs)
        error('Localization error: The sensors locations must be represented by cell array.')
    elseif ~isvector(sensorsLocs)
        error('Localization error: The sensors locations must be a vector cell array.')
    elseif any(~cellfun(@ismatrix,sensorsLocs))
        error('Localization error: The sensors locations for each array must be represented by a matrix.')
    elseif any(cellfun(@(x)size(x,2),sensorsLocs)~=2)
        error('Localization error: The number of coordinates of the sensors locations must be 2.')
    elseif any(~cellfun(@isreal,sensorsLocs))
        error('Localization error: The sensors locations coordinates must be real numbers.')
    end
    numBS = length(sensorsLocs);
    numSensors = cellfun(@(x)size(x,1),sensorsLocs);
    
    % Check 'receivedSignals'
    if ~iscell(receivedSignals)
        error('Localization error: The received signals must be represented by cell array.')
    elseif ~isvector(receivedSignals)
        error('Localization error: The received signals must be a vector cell array.')
    elseif length(receivedSignals)~=numBS
        error(['Localization error: The received signals must be defined for ', num2str(numBS), ' base stations.'])
    elseif any(~cellfun(@ismatrix,receivedSignals))
        error('Localization error: The received signals for each array must be represented by a vector.')
    elseif any(cellfun(@(x,y)size(x,1)~=y,receivedSignals,num2cell(numSensors)))
        error('Localization error: The received signals are not defined for all antennas.')
    end
    numSamples = size(receivedSignals{1},2);
    if any(cellfun(@(x)size(x,2)~=numSamples,receivedSignals))
        error('Localization error: There must be an equal amount of observations at each antenna.')
    elseif any(~cellfun(@isnumeric,receivedSignals))
        error('Localization error: The received signals must be numbers.')
    end
    
    % Check 'bandwidth'
    if ~isscalar(bandwidth)
        error('Localization error: The sampling frequency is not a scalar.')
    elseif ~isreal(bandwidth)
        error('Localization error: The sampling frequency is not a real number.')
    elseif ~isfinite(bandwidth)
        error('Localization error: Infinite is not valid for the sampling frequency.')
    elseif bandwidth<=0
        error('Localization error: The sampling frequency must be positive.')
    end
    
    % Check 'carrierFreq'
    if ~isscalar(carrierFreq)
        error('Localization error: The sampling frequency is not a scalar.')
    elseif ~isreal(carrierFreq)
        error('Localization error: The sampling frequency is not a real number.')
    elseif ~isfinite(carrierFreq)
        error('Localization error: Infinite is not valid for the sampling frequency.')
    elseif carrierFreq<=0
        error('Localization error: The sampling frequency must be positive.')
    end
    
    % Check 'oversamplingFactor'
    if ~isscalar(oversamplingFactor)
        error('Localization error: The oversampling factor is not a scalar.')
    elseif ~isreal(oversamplingFactor) || oversamplingFactor-round(oversamplingFactor)~=0
        error('Localization error: The oversampling factor is not an integer.')
    elseif ~isfinite(oversamplingFactor)
        error('Localization error: Infinite is not valid for the oversampling factor.')
    elseif oversamplingFactor<=0
        error('Localization error: The oversampling factor must be positive.')
    end
    
    % Check 'resolution'
    if ~isscalar(resolution)
        error('Localization error: The resolution is not a scalar.')
    elseif ~isreal(resolution)
        error('Localization error: The resolution is not a real number.')
    elseif ~isfinite(resolution)
        error('Localization error: Infinite is not valid for the resolution.')
    elseif resolution<=0
        error('Localization error: The resolution must be positive.')
    end
    
    % Check 'diamaterSearchArea'
    if ~isscalar(diamaterSearchArea)
        error('Localization error: The diamater of the search area is not a scalar.')
    elseif ~isreal(diamaterSearchArea)
        error('Localization error: The cdiamater of the search area is not a real number.')
    elseif ~isfinite(diamaterSearchArea)
        error('Localization error: Infinite is not a valid value for the diamater of the search area.')
    elseif diamaterSearchArea<=0
        error('Localization error: The diamater of the search area must be positive.')
    end
    
    % Check 'numIter'
    if ~isscalar(numIter)
        error('Localization error: The number grid refinement iterations is not a scalar.')
    elseif ~isreal(numIter) || numIter-round(numIter)~=0
        error('Localization error: The number grid refinement iterations is not an integer.')
    elseif ~isfinite(numIter)
        error('Localization error: Infinite is not valid for the number grid refinement iterations.')
    elseif numIter<=0
        error('Localization error: The number grid refinement iterations must be positive.')
    end
    
    
% -----------------------------------------------------------
%   Algorithm
% -----------------------------------------------------------
    
    % Generate grid of locations
    [x,y] = meshgrid(-diamaterSearchArea/2:resolution:diamaterSearchArea/2);
    gridLocs = [x(:),y(:)];
    numLocs = size(gridLocs,1);
    
    for i=1:numIter
        
        if i>1
            resolution = resolution/2;
            gridLocs = refineGrid(resolution,mobileUserLocEstimate);
            numLocs = size(gridLocs,1);
        end
    
        % Compute angles and delays from all grid locations
        baseStationsLocs = cell2mat(cellfun(@(x)mean(x,1),sensorsLocs,'UniformOutput',false));
        temp = repmat(permute(gridLocs,[3,2,1]),[numBS,1,1])-repmat(baseStationsLocs,[1,1,numLocs]);
        distBS2locs = sqrt(squeeze(sum(temp.^2,2)));
        angleBS2locs = squeeze(atan(temp(:,2,:)./temp(:,1,:)) +(temp(:,1,:)<0)*pi);

        % Global match filter
        offsetSamples = ceil(3*oversamplingFactor/pi*sqrt(2*log(2)));
        likelihood = zeros(1,numLocs);
        for l=1:numBS
            atom = oversamplingFactor^(-1/2)*(pi/log(2))^(1/4)...
                *exp(-((0:numSamples-1).'*ones(1,numLocs)./oversamplingFactor-bandwidth*ones(numSamples,1)*distBS2locs(l,:)./3e8-offsetSamples/oversamplingFactor).^2.*pi.^2./2./log(2));

            spatialData = receivedSignals{l}*conj(atom);

            likelihood = likelihood +abs(sum(conj(angleToSteering(sensorsLocs{l},angleBS2locs(l,:),carrierFreq)).*spatialData,1)).^2/numSensors(l);
        end
        
        [~,ind] = max(likelihood);
        mobileUserLocEstimate = gridLocs(ind,:);
    
    end
end

function [activeLocs] = refineGrid(locRes,location)
    
    % Refine grid of locations
    [extraX,extraY] = meshgrid(-2*locRes:locRes:2*locRes,-2*locRes:locRes:2*locRes);
    activeLocs = ones(25,1)*location +[extraX(:),extraY(:)];
end