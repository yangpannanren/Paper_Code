function [ mobileUserLocEstimate ] = Stansfield( baseStationsLocs, TOAs, AOAs )
%DIRECT Summary of this function goes here
%   Detailed explanation goes here

% Coded by Nil Garcia

%   Inputs:
%   - baseStationsLocs [m]
%       Array of size (numBS)x2 matrix.
%   - TOAs [s]
%       The estimated TOA's at all BS's.
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
    
    % Check 'AOAs'
    if ~isvector(AOAs)
        error('Localization error: The TOA''s must be represented by a vector.')
    elseif any(length(AOAs)~=numBS)
        error('Localization error: The TOA''s must be specified by a vector of size 1x2.')
    elseif ~isreal(AOAs)
        error('Localization error: The TOA''s must be real numbers.')
    elseif any(~isfinite(AOAs))
        error('Localization error: The TOA''s cannot be infinite.')
    end
    
    
% -----------------------------------------------------------
%   Localization
% -----------------------------------------------------------

    select = ~isnan(TOAs) & ~isnan(AOAs);

    TOAs = TOAs(select);
    AOAs = AOAs(select);
    baseStationsLocs = baseStationsLocs(select,:);
    numBS = sum(select);
    
    mobileUserLocEstimate = [];
    if numBS>=2
        A = [sin(AOAs(:)), -cos(AOAs(:))];
        b = sum(baseStationsLocs.*A,2);
        R = diag((TOAs*3e8).^2);

        mobileUserLocEstimate = (A.'/R*A)\A.'/R*b;
        mobileUserLocEstimate = mobileUserLocEstimate.';
    end
end