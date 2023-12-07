function [ mobileUserLocEstimate ] = IV( baseStationsLocs, AOAs )
%DIRECT Summary of this function goes here
%   Detailed explanation goes here

% Coded by Nil Garcia

%   Inputs:
%   - baseStationsLocs [m]
%       Array of size (numBS)x2 matrix.
%   - AOAs [radians]
%       The estimated AOA's at all BS's.
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
    
    % Check 'AOAs'
    if ~isvector(AOAs)
        error('Localization error: The AOA''s must be represented by a vector.')
    elseif any(length(AOAs)~=numBS)
        error('Localization error: The AOA''s must be specified by a vector of size 1x2.')
    elseif ~isreal(AOAs)
        error('Localization error: The AOA''s must be real numbers.')
    elseif any(~isfinite(AOAs))
        error('Localization error: The AOA''s cannot be infinite.')
    end
    
    
% -----------------------------------------------------------
%   Localization
% -----------------------------------------------------------

    select = ~isnan(AOAs);

    AOAs = AOAs(select);
    baseStationsLocs = baseStationsLocs(select,:);
    numBS = sum(select);
    
    mobileUserLocEstimate = [];
    if numBS>=2
        A = [sin(AOAs(:)), -cos(AOAs(:))];
        b = sum(baseStationsLocs.*A,2);

        mobileUserLocEstimate = (A.'*A)\A.'*b;
        mobileUserLocEstimate = mobileUserLocEstimate.';

        % Estimate angles
        temp = ones(numBS,1)*mobileUserLocEstimate -baseStationsLocs;
        newAOAs = squeeze(atan(temp(:,2)./temp(:,1)) +(temp(:,1)<0)*pi);

        G = [sin(newAOAs(:)), -cos(newAOAs(:))];

        mobileUserLocEstimate = (G.'*A)\G.'*b;
        mobileUserLocEstimate = mobileUserLocEstimate.';
    end
end