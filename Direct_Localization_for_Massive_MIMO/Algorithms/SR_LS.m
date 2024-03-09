function [ mobileUserLocEstimate ] = SR_LS( baseStationsLocs, TOAs )
%DIRECT Summary of this function goes here
%   Detailed explanation goes here

% Coded by Nil Garcia

%   Inputs:
%   - baseStationsLocs [m]
%       Array of size (numBS)x2 matrix.
%   - TOAs [s]
%       The estimated TOA at all BS's.
%
%   Outputs:
%   - mobileUserLocEstimate [m]
%       Matrix of size 1x2 with the coordinates of the location of the MU.
    

% -----------------------------------------------------------
%   Check that the inputs comply with the correct format.
% -----------------------------------------------------------

    % Check 'baseStationsLocs'
    if ~ismatrix(baseStationsLocs)
        error('Localization error: The locations of the base stations must be represented by a matrix.')
    elseif size(baseStationsLocs,2)~=2
        error('Localization error: The number of coordinates of the base stations must be two.')
    elseif any(~isreal(baseStationsLocs))
        error('Localization error: The base stations locations coordinates must be real numbers.')
    end
    numBS = size(baseStationsLocs,1);
    
    % Check 'TOAs'
    if ~isvector(TOAs)
        error('Localization error: The TOA''s must be represented by a vector.')
    elseif any(length(TOAs)~=numBS)
        error('Localization error: The TOA''s must be specified by a vector of size 1x2.')
    elseif ~isreal(TOAs)
        error('Localization error: The TOA''s must be real numbers.')
    elseif any(~isfinite(TOAs))
        error('Localization error: The TOA''s cannot be infinite.')
    end
    
    
% -----------------------------------------------------------
%   Localization
% -----------------------------------------------------------

    select = ~isnan(TOAs);

    TOAs = TOAs(select);
    baseStationsLocs = baseStationsLocs(select,:);
    numBS = sum(select);
    
    mobileUserLocEstimate = [];
    if numBS>2
        matrixA = [-2*baseStationsLocs, ones(numBS,1)];
        vectorB = (TOAs(:)*3e8).^2 -sum(baseStationsLocs.^2,2);
        matrixD = [1,0,0; 0,1,0; 0,0,0];
        vectorF = [0;0;-.5];

        % Compute lamba by bisection
        LBlambda = -1/eigs(matrixD,matrixA.'*matrixA,1);
        UBlambda = inf;
        lambda = nan;
        d = eig(matrixA.'*matrixA);
        precision = min(d)/1e3;
        ite = 0;
        while isnan(lambda)
            if UBlambda==inf
                testLambda = LBlambda+exp(ite);
                yHat = (matrixA.'*matrixA+testLambda*matrixD)\(matrixA.'*vectorB-testLambda*vectorF);
                func = yHat.'*matrixD*yHat +2*vectorF.'*yHat;
                if func==0
                    lambda = testLambda;
                elseif func<0
                    UBlambda = testLambda;
                end
            else
                testLambda = (LBlambda+UBlambda)/2;
                yHat = (matrixA.'*matrixA+testLambda*matrixD)\(matrixA.'*vectorB-testLambda*vectorF);
                func = yHat.'*matrixD*yHat +2*vectorF.'*yHat;
                if func==0 || UBlambda-LBlambda<=precision
                    lambda = testLambda;
                elseif func<0
                    UBlambda = testLambda;
                else
                    LBlambda = testLambda;
                end
            end
            ite = ite+1;
        end

        mobileUserLocEstimate = yHat(1:2).';
    end
    
end