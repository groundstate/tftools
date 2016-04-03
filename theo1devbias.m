function [ bias, theo1corr ] = theo1devbias(m,noisefn,theo1 )
%THEO1DEVBIAS Calculates the bias correction for the Theo1
% deviation, given the noise type, and optionally corrects a vector
% of Theo1 deviation values
%
%   Usage: [ bias ] = THEO1DEVBIAS(m,noisefn) 
%   Inputs:
%     
%     m:  vector of averaging intervals (tau*rate)
%     noisefn: "white pm","flicker pm","white fm","flicker fm",
%             "random walk fm"
%     theo1: optional: data to be corrected
%        
%   Outputs:
%     bias:  vector of bias corrections, defined as
%            bias(m,tau0)= ADEV(tau)/THEO1DEV(m,tau0) - 1
% theor1corr: corrected theo1
%
%   See also adev,freq2phase, hdev, mdev, tau2m, tdev, theo1dev, totdev.

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

% Estimated degrees of freedom
% Ref: D A Howe "TheoH: a hybrid, high-confidence statistic
% that improves on the Allan deviation" Metrologia v43 S322-S331 (2006)

bias = zeros(size(m));
theo1corr = ones(size(m));

if (strcmpi(noisefn,'white pm'))
  bias = sqrt(0.09 + 0.74./ (m .^0.4)) - 1;
elseif (strcmpi(noisefn,'flicker pm'))
  bias = sqrt(0.14 + 0.82 ./ (m .^0.3)) - 1; 
elseif (strcmpi(noisefn,'white fm'))
   % nothing to do
elseif (strcmpi(noisefn,'flicker fm'))
   bias = sqrt(1.87  - 1.05 ./ (m .^0.79)) - 1;
elseif (strcmpi(noisefn,'random walk fm'))
   bias = sqrt(2.70 - 1.53 ./ (m .^0.85)) - 1;
else
   error('tftools:theo1devbias:BadInput', ['unknown noise type: ' noisefn]);
end

if (nargin == 3)
    theo1corr = theo1 .* (bias+1);
end

end

