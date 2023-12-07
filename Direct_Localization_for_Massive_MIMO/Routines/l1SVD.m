function [ estimAOA ] = l1SVD( sensorsLocs, snapshot, carrierFreq, noisePower)

% Coded by Nil Garcia

%   Inputs:
%   - sensorsLocs [m]
%       Matrix of size (numSensors)x2 matrix.
%   - snapshot
%       Vector of samples of length (numSensors).
%   - carrierFreq [Hz]
%       The carrier frequency of the emitted signal.
%   - noisePower
%       The noise power.
%
%   Outputs:
%   - estimAOA
%       The estimated angle of the strongest arrival.
    

% -----------------------------------------------------------
%   Check that the inputs comply with the correct format.
% -----------------------------------------------------------

    % Check 'sensorsLocs'
    if ~ismatrix(sensorsLocs)
        error('l1-SVD error: The sensors locations must be represented by a matrix.')
    elseif size(sensorsLocs,2)~=2
        error('l1-SVD error: The number of coordinates of the sensors locations must be 2.')
    elseif ~isreal(sensorsLocs)
        error('l1-SVD error: The sensors locations coordinates must be real numbers.')
    end
    numSensors = size(sensorsLocs,1);
    
    % Check 'snapshot'
    if ~isvector(snapshot)
        error('l1-SVD error: The spatial samples must be represented by a vector.')
    elseif length(snapshot)~=numSensors
        error('l1-SVD error: The spatial samples are not defined for all antennas.')
    elseif ~isnumeric(snapshot)
        error('l1-SVD error: The spatial samples coordinates must be numbers.')
    end
    
    % Check 'carrierFreq'
    if ~isscalar(carrierFreq)
        error('l1-SVD error: The carrier frequeny is not a scalar.')
    elseif ~isreal(carrierFreq)
        error('l1-SVD error: The carrier frequeny is not a real number.')
    elseif ~isfinite(carrierFreq)
        error('l1-SVD error: Infinite is not a valid value for the carrier frequeny.')
    elseif carrierFreq<=0
        error('l1-SVD error: The carrier frequeny must be positive.')
    end
    
    % Check 'noisePower'
    if ~isscalar(noisePower)
        error('l1-SVD error: The noise power is not a scalar.')
    elseif ~isreal(noisePower)
        error('l1-SVD error: The noise power is not a real number.')
    elseif ~isfinite(noisePower)
        error('l1-SVD error: Infinite is not a valid value for the noise power.')
    elseif noisePower<0
        error('l1-SVD error: The noise power cannot be negative.')
    end
    
    
% -----------------------------------------------------------
%   l1-SVD
% -----------------------------------------------------------
    
    estimAOA = [];
    
    maxError = sqrt(noisePower/2*chi2inv(.99,2*numSensors));
    
    if norm(snapshot)>maxError

        angles = 2*pi/1e3:2*pi/1e3:2*pi;
        numAngles = length(angles);

        numNonZeros = 4*numSensors*numAngles +2*numSensors;

        rowIndices = zeros(numNonZeros,1);
        columnIndices = zeros(numNonZeros,1);
        values = zeros(numNonZeros,1);
        offset = 0;

        select = offset+(1:4*numSensors*numAngles);
        rowIndices(select) = (1:2*numSensors).'*ones(1,2*numAngles);
        temp =  reshape(1:3*numAngles,[3,numAngles]);
        temp = temp([2,3],:);
        columnIndices(select) = ones(2*numSensors,1)*temp(:).';
        atoms = angleToSteering( sensorsLocs, angles, carrierFreq );
        values(select) = [real(atoms);imag(atoms);-imag(atoms);real(atoms)];
        offset = offset +4*numSensors*numAngles;

        select = offset+(1:2*numSensors);
        rowIndices(select) = 1:2*numSensors;
        columnIndices(select) = 3*numAngles +1 +(1:2*numSensors);
        values(select) = 1;
        offset = offset +2*numSensors;

        % Data
        rhsVector = zeros(2*numSensors,1);
        rhsVector(:) = [real(snapshot(:));imag(snapshot(:))];

        % Specify the non-conic part of the problem.
        prob.c = zeros(3*numAngles+1+2*numSensors,1);
        prob.c(1:3:3*numAngles) = 1;
        prob.a = sparse(rowIndices,columnIndices,values,2*numSensors,3*numAngles+1+2*numSensors,numNonZeros);

        prob.blc = rhsVector;
        prob.buc = rhsVector;

        prob.blx = -inf(3*numAngles+1+2*numSensors,1);
        prob.bux = inf(3*numAngles+1+2*numSensors,1); 
        prob.bux(3*numAngles+1) = maxError;

        % Specify the cones.
        prob.cones.type = zeros(numAngles+1,1);
        prob.cones.sub = 1:3*numAngles+1+2*numSensors;
        prob.cones.subptr = [1:3:3*numAngles,3*numAngles+1];

        [~,res]=mosekopt('minimize echo(0)',prob);

        if  strcmpi(res.sol.itr.prosta,'PRIMAL_AND_DUAL_FEASIBLE')
            temp = res.sol.itr.xx(1:3:3*numAngles);
%             % TEST
%             figure
%             semilogy(180/pi*angles,temp)
%             % TEST
            [~,ind] = findpeaks(temp,'MINPEAKHEIGHT',max(temp)*1e-7,'SORTSTR','none');
            amplitudes = pinv(atoms(:,ind))*snapshot;
            [~,ix] = max(abs(amplitudes));
            estimAOA = angles(ind(ix));
        end
    end

end

