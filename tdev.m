function [ dev,deverr,ndev,new_tau ] = tdev( x,rate,tau,phase,gaps)
%TDEV Calculate time deviation of phase/
% or fractional frequency data
%   Usage: [ dev,deverr,ndev,new_tau ] = TDEV(x,rate,tau,phase,gaps) 
%   Inputs:
%     x       input time series
%     rate    sampling rate, in Hz
%     tau     averaging intervals 
%     phase   data is phase (=1), optional argument, default=1
%     gaps    data contains gaps, tagged with NaN, optional argument, default=0
%   Outputs:
%     dev     time deviations 
%     deverr  uncertainties
%     ndev    number of samples used
%     new_tau tau values that were used
%
%   See also adev, freq2phase, markgaps, mdev, tau2m, totdev.
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
% Credits:

if (nargin > 5)
    error('tftools:tdev:TooManyInputs', 'requires at most 2 optional arguments');
end;

defaults = {1 0};

switch nargin
    case 3
        [phase gaps] = defaults{:};
    case 4
        gaps=defaults{2};
end;

if (phase == 0)
    x=freq2phase(x,rate);
end;

[dev, deverr, ndev, new_tau] = mdev(x,rate,tau,phase,gaps);
dev = dev .* new_tau/sqrt(3.0);
deverr = deverr .* new_tau/sqrt(3.0);

end

