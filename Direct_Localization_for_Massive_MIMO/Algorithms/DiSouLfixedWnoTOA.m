function [ mobileUserLocEstimate, detectedLOSbs ] = DiSouLfixedWnoTOA( sensorsLocs, snapshot, noisePower, locGridRes, carrierFreq, diamaterSearchArea, wSquare )
%DIRECT Summary of this function goes here
%   Detailed explanation goes here

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
%   - noisePower
%       The noise power.
%   - locGridRes [m]
%       The starting cell size for the grid refinement procedure.
%   - carrierFreq [Hz]
%       The carrier frequency of the emitted signal.
%   - diamaterSearchArea [m]
%
%   Outputs:
%   - mobileUserLocEstimate [m]
%       Matrix of size 1x2 with the coordinates of the location of the MU.
%   - detectedLOSbs
%       Indexes of the BS's which have been detected to include a LOS component.
    

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
    if numBS<2
        mobileUserLocEstimate = nan;
        detectedLOSbs = [];
        return;
    end
    numSensors = cellfun(@(x)size(x,1),sensorsLocs);
    
    % Check 'snapshot'
    if ~iscell(snapshot)
        error('Localization error: The spatial samples must be represented by cell array.')
    elseif ~isvector(snapshot)
        error('Localization error: The spatial samples must be a vector cell array.')
    elseif length(snapshot)~=numBS
        error(['Localization error: The spatial samples must be defined for ', num2str(numBS), ' base stations.'])
    elseif any(~cellfun(@isvector,snapshot))
        error('Localization error: The spatial samples for each array must be represented by a vector.')
    elseif any(cellfun(@(x,y)length(x)~=y,snapshot,num2cell(numSensors)))
        error('Localization error: The spatial samples are not defined for all antennas.')
    elseif any(~cellfun(@isnumeric,snapshot))
        error('Localization error: The spatial samples must be numbers.')
    end
    
    % Check 'noisePower'
    if ~isscalar(noisePower)
        error('Localization error: The noise power is not a scalar.')
    elseif ~isreal(noisePower)
        error('Localization error: The noise power is not a real number.')
    elseif ~isfinite(noisePower)
        error('Localization error: Infinite is not a valid value for the noise power.')
    elseif noisePower<0
        error('Localization error: The noise power cannot be negative.')
    end
    
    % Check 'locGridRes'
    if ~isscalar(locGridRes)
        error('Localization error: The grid resolution is not a scalar.')
    elseif ~isreal(locGridRes)
        error('Localization error: The grid resolution is not a real number.')
    elseif ~isfinite(locGridRes)
        error('Localization error: Infinite is not a valid value for the grid resolution.')
    elseif locGridRes<=0
        error('Localization error: The grid resolution must be positive.')
    end
    
    % Check 'carrierFreq'
    if ~isscalar(carrierFreq)
        error('Localization error: The carrier frequeny is not a scalar.')
    elseif ~isreal(carrierFreq)
        error('Localization error: The carrier frequeny is not a real number.')
    elseif ~isfinite(carrierFreq)
        error('Localization error: Infinite is not a valid value for the carrier frequeny.')
    elseif carrierFreq<=0
        error('Localization error: The carrier frequeny must be positive.')
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
    
    
% -----------------------------------------------------------
%   Localization
% -----------------------------------------------------------
    
    maxError = sqrt(noisePower/2*chi2inv(.99,2*sum(numSensors)));
    
    mobileUserLocEstimate = [];
    detectedLOSbs = [];

    if  sqrt(sum(cellfun(@(x)norm(x)^2,snapshot))) > maxError
        
        baseStationsLocs = cell2mat(cellfun(@(x)mean(x,1),sensorsLocs,'UniformOutput',false));

        angleGridRes = 2*locGridRes/diamaterSearchArea;
        angleGridRes = 2*pi/ceil(2*pi/angleGridRes);

        [x,y] = meshgrid(-diamaterSearchArea/2:locGridRes:diamaterSearchArea/2);
        gridLocs = [x(:),y(:)];
        gridAngles = cell(numBS,1);
        [gridAngles{:}] = deal(angleGridRes:angleGridRes:2*pi);
        
        
% -----------------------------------------------------------
%       % Grid refinement loop
% -----------------------------------------------------------
        
        oldCost = inf;
        newCost = inf;

        activeLocs = gridLocs;
        activeAngles = gridAngles;

        locRes = locGridRes;
        angleRes = angleGridRes;

        while newCost <= .999*oldCost % Stop when practical convergence is achieved

            if newCost<inf % 'newCost' is equal to inf at the 1st iteration only
                locRes = locRes/2;
                angleRes = angleRes/2;
                [activeLocs,activeAngles] = refineGrids(locRes,angleRes,activeLocs,activeAngles,angleBS2locs);
            end

            % Compute AOA and TOA from all grid locations
            temp = repmat(permute(activeLocs,[3,2,1]),[numBS,1,1]) -repmat(baseStationsLocs,[1,1,size(activeLocs,1)]);
            angleBS2locs = permute(atan(temp(:,2,:)./temp(:,1,:)) +(temp(:,1,:)<0)*pi,[1,3,2]);

            [optVal,locVar,angleVar] = solveOptProblem(sensorsLocs,snapshot,activeLocs,activeAngles,angleBS2locs,maxError,wSquare,carrierFreq);

            oldCost = newCost;
            newCost = optVal;

            maxVal = max([abs(locVar(:).'),abs(cell2mat(angleVar))]);

            for l=1:numBS
                select = abs(angleVar{l}) > maxVal*1e-7;
                activeAngles{l} = activeAngles{l}(select);
            end

            select = abs(locVar)>maxVal*1e-7;
            if sum(select(:))>0
                [~,ind] = max(sum(abs(locVar).^2,1));
                mobileUserLocEstimate = activeLocs(ind,:);
                detectedLOSbs = find(select(:,ind));
                activeLocs = activeLocs(any(select,1),:);
                angleBS2locs = angleBS2locs(:,any(select,1));
            else
                break
            end
        end
    end

end


% -----------------------------------------------------------
%   Local Functions
% -----------------------------------------------------------

function [optVal,locVar,angleVar] = solveOptProblem(sensorsLocs,snapshot,activeLocs,activeAngles,angleBS2locs,maxError,wSquare,carrierFreq)
    
    numBS = length(sensorsLocs);
    numSensors = cellfun(@(x)size(x,1),sensorsLocs);
    numLocs = size(activeLocs,1);
    numAngles = cellfun(@length,activeAngles);

    numNonZeros = 4*sum(numSensors.*numAngles) +4*sum(numSensors)*numLocs +2*sum(numSensors);
    
    rowIndices = zeros(numNonZeros,1);
    columnIndices = zeros(numNonZeros,1);
    values = zeros(numNonZeros,1);
    offset = 0;

    % Angle variables
    for l=1:numBS
        select = offset+(1:4*numSensors(l)*numAngles(l));
        rowIndices(select) = 2*sum(numSensors(1:l-1)) +(1:2*numSensors(l)).'*ones(1,2*numAngles(l));
        temp =  reshape(1:3*numAngles(l),[3,numAngles(l)]);
        temp = temp([2,3],:);
        columnIndices(select) = 3*sum(numAngles(1:l-1)) +ones(2*numSensors(l),1)*temp(:).';
        atom = angleToSteering(sensorsLocs{l},activeAngles{l},carrierFreq);
        values(select) = [real(atom);imag(atom);-imag(atom);real(atom)];
        offset = offset +4*numSensors(l)*numAngles(l);
    end
    
    % Location variables
    for l=1:numBS
        select = offset+(1:4*numSensors(l)*numLocs);
        rowIndices(select) = 2*sum(numSensors(1:l-1)) +(1:2*numSensors(l)).'*ones(1,2*numLocs);
        temp =  reshape(1:(1+2*numBS)*numLocs,[1+2*numBS,numLocs]);
        temp = temp([1+l,1+numBS+l],:);
        columnIndices(select) = 3*sum(numAngles) +ones(2*numSensors(l),1)*temp(:).';
        atom = angleToSteering(sensorsLocs{l},angleBS2locs(l,:),carrierFreq);
        values(select) = [real(atom);imag(atom);-imag(atom);real(atom)];
        offset = offset +4*numSensors(l)*numLocs;
    end
    
    select = offset+(1:2*sum(numSensors));
    rowIndices(select) = 1:2*sum(numSensors);
    columnIndices(select) = 3*sum(numAngles) +(1+2*numBS)*numLocs +1 +(1:2*sum(numSensors));
    values(select) = 1;
    offset = offset +2*sum(numSensors);
    
    % Data
    rhsVector = zeros(2*sum(numSensors),1);
    for l=1:numBS
        rhsVector(2*sum(numSensors(1:l-1)) +(1:numSensors(l))) = real(snapshot{l});
        rhsVector(2*sum(numSensors(1:l-1)) +numSensors(l) +(1:numSensors(l))) = imag(snapshot{l});
    end
    
    prob.c = zeros(3*sum(numAngles)+(1+2*numBS)*numLocs+1+2*sum(numSensors),1);
    prob.c(1:3:3*sum(numAngles)) = 1;
    prob.c(3*sum(numAngles) +(1:1+2*numBS:(1+2*numBS)*numLocs)) = sqrt(wSquare);
    
    prob.a = sparse(rowIndices,columnIndices,values,2*sum(numSensors),3*sum(numAngles)+(1+2*numBS)*numLocs+1+2*sum(numSensors),numNonZeros);
    
    prob.blc = rhsVector;
    prob.buc = rhsVector;

    prob.blx = -inf(3*sum(numAngles)+(1+2*numBS)*numLocs+1+2*sum(numSensors),1);
    prob.bux = inf(3*sum(numAngles)+(1+2*numBS)*numLocs+1+2*sum(numSensors),1); 
    prob.bux(3*sum(numAngles)+(1+2*numBS)*numLocs+1) = maxError;

    % Specify the cones.
    prob.cones.type = zeros(sum(numAngles)+numLocs+1,1);
    prob.cones.sub = 1:3*sum(numAngles) +(1+2*numBS)*numLocs +1+2*sum(numSensors);
    prob.cones.subptr = [1:3:3*sum(numAngles),3*sum(numAngles)+(1:1+2*numBS:(1+2*numBS)*numLocs),3*sum(numAngles)+(1+2*numBS)*numLocs+1];
    
    [~,res]=mosekopt('minimize echo(0)',prob);
    
    optVal = res.sol.itr.pobjval;

    angleVar = reshape(res.sol.itr.xx(1:3*sum(numAngles)),[3,sum(numAngles)]);
    angleVar = angleVar(2,:) +1i*angleVar(3,:);
    angleVar = mat2cell(angleVar,1,numAngles);
    
    locVar = reshape(res.sol.itr.xx(3*sum(numAngles)+(1:(1+2*numBS)*numLocs)),[1+2*numBS,numLocs]);
    locVar = locVar(1+(1:numBS),:) +1i*locVar(1+numBS+(1:numBS),:);
end

function [activeLocs,activeAngles] = refineGrids(locRes,angleRes,activeLocs,activeAngles,angleBS2locs)

    numBS = length(activeAngles);
    
    % Refine grid of angles
    for l=1:numBS
        newAngles = zeros(length(activeAngles{l}),5);
        for i=1:length(activeAngles{l})
            newAngles(i,:) = round(activeAngles{l}(i)/angleRes)*angleRes +(-2*angleRes:angleRes:2*angleRes);
        end
        activeAngles{l} = unique(mod(newAngles(:),2*pi));
        activeAngles{l} = unique([activeAngles{l};round(mod(angleBS2locs(l,:).',2*pi)/angleRes)*angleRes;]);
    end
    
    % Refine grid of locations
    [extraX,extraY] = meshgrid(-2*locRes:locRes:2*locRes,-2*locRes:locRes:2*locRes);
    if ~isempty(activeLocs)
        newLocs = zeros(25,2,size(activeLocs,1));
        for i=1:size(activeLocs,1)
            newLocs(:,:,i) = ones(25,1)*activeLocs(i,:) +[extraX(:),extraY(:)];
        end
        activeLocs = unique(reshape(permute(newLocs,[1,3,2]),[25*size(activeLocs,1),2]),'rows');
    end
end
