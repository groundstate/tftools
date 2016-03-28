function [ lowerlim, upperlim ] = adevconf( N,m, confidence,noisefn )
%ADEVCONF Calculates the two-sided confidence interval for the overlapping
% Allan deviation, given the noise type.
%
%   Usage: [ minerr,maxerr ] = ADEVCONF(n,m,confidence,noisefn) 
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
%   given the nominal ADEV, the confidence interval is
%   [1.0-lowerlim, 1.0+upperlim] * ADEV
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
     error('tftools:adevconf:BadInput','confidence must be between 0 and 100');
end;

% Estimated degrees of freedom
% Ref: NIST SP1065, Table 5 p40
if (strcmpi(noisefn,'white pm'))
   edf = (N+1)*(N-2*m) ./ (2.0*(N-m));
elseif (strcmpi(noisefn,'flicker pm'))
   edf=exp((log((N -1) ./ (2*m)) .* log((2*m+1)*(N-1)/4.0)) .^ 0.5);
elseif (strcmpi(noisefn,'white fm'))
   m2 = 4.0 * m .^ 2;
   edf = ( 3.0*(N-1) ./ (2.0*m) - 2.0*(N-2)/N ) .* m2 ./ (m2 + 5 ); 
elseif (strcmpi(noisefn,'flicker fm'))
   if (m == 1)
       edf = 2.0*((N - 2)^2)/(2.3*N - 4.9);
   else
       edf = 5.0*N^2 ./ (4.0*m .* (N + 3*m));
   end
elseif (strcmpi(noisefn,'random walk fm'))
    edf = ((N-2) ./ m) .* ( (N-1)^2  - 3*m*(N-1) + 4.0*m .^2) / (N-3)^2;
else
   error('tftools:adevconf:BadInput', ['unknown noise type: ' noisefn]);
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

% Validation data from IEEE 1139-1999 p25 Table D.2
% N=1025;
% m=[2 8 32 128];
% wpm_min=[2.9 2.9 3.0 3.1]; % white phase modulation
% wpm_max=[3.2 3.2 3.4 3.6];
% fpm_min=[2.9 3.6 5.2 8.4]; % flicker phase modulation
% fpm_max=[3.1 4.0 6.1 11.0];
% wfm_min=[2.8 4.8 8.8 16]; % white frequency modulation
% wfm_max=[3.0 5.6 12 32];
% 