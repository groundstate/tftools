function [ tn, xn ] = markgaps( t, x, dt)
%markgaps finds and marks gaps in time series data
%   x  measurement time
%   y  measurements
%   dt nominal data spacing
%   Returns expanded input sequences, with missing data flagged by NaN

% Measurement times may be inexact, so round them to the nominal spacing
trnd = round(t/dt)*dt;
% Calculate the expected measurement times
texp = trnd(1):dt:trnd(length(trnd));

[tf, loc] = ismember(trnd, texp);

tn = nan(size(texp));
xn = tn;

% Replace NaNs where there is a measured value
tn(loc)= t;
xn(loc)= x;

end

