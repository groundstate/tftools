function [ lowerlim, upperlim ] = totdevci( N,m, confidence,noisefn )
%TOTDEVCI Calculates the two-sided confidence interval for the total
% deviation, given the noise type.
%
%   Usage: [ lowerlim,upperlim ] = TOTDEVCI(n,m,confidence,noisefn) 
%   Inputs:
%     
%     N:  number of samples (assuming phase data)
%     m:  vector of averaging intervals (tau*rate)
%     confidence: (in %)
%     noisefn: "white fm","flicker fm","random walk fm"
%   Outputs:
%     lowerlim:  vector of lower limits of the error estimate
%     upperlim:  vector of upper limits of the error estimate
%   The limits are expressed as a fractional error so that 
%   given the nominal TOTDEV, the confidence interval is
%   [1.0-lowerlim, 1.0+upperlim] * TOTDEV
%
%   See also freq2phase, hdev, mdev, tau2m, tdev, totdev.

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
     error('tftools:totdevci:BadInput','confidence must be between 0 and 100');
end;

% Estimated degrees of freedom
% Ref: NIST SP1065, Table 7 p41
if (strcmpi(noisefn,'white fm'))
   b=1.5;
   c=0.0;
elseif (strcmpi(noisefn,'flicker fm'))
   b=1.17;
   c=0.22;
elseif (strcmpi(noisefn,'random walk fm'))
   b=0.93;
   c=0.36;
else
   error('tftools:totdevci:BadInput', ['unknown noise type: ' noisefn]);
end

edf = b*N ./ m - c;

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