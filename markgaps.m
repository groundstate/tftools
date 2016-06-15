function [ tn, xn ] = markgaps( t, x, dt)
%MARKGAPS finds and marks gaps in time series data
% Usage: [ tn, xn ] = MARKGAPS( t, x, dt)
% Inputs:
%   t     measurement times
%   x     measurements
%   dt    nominal data interval (s)
% Outputs:
%   Returns expanded input sequences, with missing data flagged by NaN
%
%   See also adev, freq2phase, hdev, mdev, tau2m, tdev, totdev.
%

% The MIT License (MIT)
% 
% Copyright (c) 2016 Michael J. Wouters
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
%

% Measurement times may be inexact, so round them to the nominal spacing
trnd = round(t/dt)*dt;
% Calculate the expected measurement times
texp = trnd(1):dt:trnd(length(trnd)); % a row vector
if (iscolumn(t))
    texp=texp';
end
[~, loc] = ismember(trnd, texp);

tn = nan(size(texp));
xn = tn;

% Replace NaNs where there is a measured value
tn(loc)= t;
xn(loc)= x;

end

