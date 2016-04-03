function [ lowerlim, upperlim ] = theo1devci( N,m, confidence,noisefn )
%THEO1DEVCI Calculates the two-sided confidence interval for the Theo1
% deviation, given the noise type.
%
%   Usage: [ minerr,maxerr ] = THEO1DEVCI(N,m,confidence,noisefn) 
%   Inputs:
%     
%     N:  number of samples (assuming phase data)
%     m:  vector of averaging intervals (tau*rate)
%     confidence: (in %)
%     noisefn: "white pm","flicker pm","white fm","flicker fm",
%             "random walk fm"
%   Outputs:
%     lowerlim:  vector of lower limits of the error estimate
%     upperlim:  vector of upper limits of the error estimate
%   The limits are expressed as a fractional error so that 
%   given the nominal THEO1DEV, the confidence interval is
%   [1.0 - lowerlim, 1.0 + upperlim] * THEO1DEV
%   
%   The limits are only valid for a sampling period tau0 <= T/10, where T
%   is the sampled interval.
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

if (confidence <=0 || confidence >=100)
     error('tftools:theo1devci:BadInput','confidence must be between 0 and 100');
end;

% Estimated degrees of freedom
% Ref: D A Howe "TheoH: a hybrid, high-confidence statistic
% that improves on the Allan deviation" Metrologia v43 S322-S331 (2006)
if (strcmpi(noisefn,'white pm'))
   edf =  m./(m + 1.52) .* (0.86 * (N+1)*(N-m) ./ (N-0.75*m));
elseif (strcmpi(noisefn,'flicker pm'))
   edf = m ./ (m+0.4) .* ...
   (5.54*N^2 - 5.52*N*m + 10.727*m) ./ ((m+48.8) .^ 0.5 .* (N - 0.75*m));
elseif (strcmpi(noisefn,'white fm'))
   edf = m .^1.5 ./ (m .^1.5 + 8) .* ((5.5*N + 1.07) ./ m - (3.1*N + 6.5)/N); 
elseif (strcmpi(noisefn,'flicker fm'))
   edf = m .^ 3 ./ (m .^ 3 + 5.45) .* ( (2.7*N^2 -1.3*N*m -3.5*m) ./ (m*N));
elseif (strcmpi(noisefn,'random walk fm'))
    edf = ((4.4*N -2) ./ (2.175*m)) .* ...
        ((4.4*N-1)^2 - 6.45*m*(4.4*N-1) + 6.413*m.^2)/((4.4*N -3)^2);
else
   error('tftools:theo1devci:BadInput', ['unknown noise type: ' noisefn]);
end

lowerlim = zeros(size(m));
upperlim = zeros(size(m));
for i=1:length(m)
   clim = (1- confidence/100.0)/2.0;
   kind = exist('chi2inv');
   if (kind ~= 2)
      upperlim(i)=sqrt(edf(i)/approxchi2inv(clim,edf(i)))-1.0;
      lowerlim(i)=1.0-sqrt(edf(i)/approxchi2inv(1.0-clim,edf(i)));
   else
      upperlim(i)=sqrt(edf(i)/chi2inv(clim,edf(i)))-1.0;
      lowerlim(i)=1.0-sqrt(edf(i)/chi2inv(1.0-clim,edf(i)));
   end
end

end

