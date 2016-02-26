function [ dev,deverr, ndev ] = adev( x,rate,tau,overlap,phase, gaps )
%ADEV Calculate non-overlapping/overlapping Allan deviation of phase/
% or fractional frequency data
%   Usage: ADEV(x,rate,tau,overlap,phase,gaps) 
%   x is the input time series
%   rate is the sampling rate, in Hz
%   tau is the averaging interval 
%   overlap = 1 means compute overlapping ADEV
%   phase = 1 means data is phase
%   gaps = 1 means data contains gaps, tagged with NaN
%
%   Frequency data is converted internally to phase
%   
%   See also mdev, totdev.

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
% This code is based on allantools.py, written by Anders Wallin
%

if (nargin > 6)
    error('tftools:adev:TooManyInputs', 'requires at most 3 optional arguments');
end;

defaults = {1 1 0};

switch nargin
    case 3
        [overlap phase gaps] = defaults{:};
    case 4
        [phase gaps]=defaults{2:3};
    case 5
        gaps=defaults{3};
end;

if (phase == 0)
    x=freq2phase(x,rate);
end;

ntau=length(tau);
% FIXME sanity check on tau

dev=zeros(1,ntau);
deverr = zeros(1,ntau);
ndev=zeros(1,ntau);

dt = 1; % for overlapping Allan deviation

for i=1:ntau 
   taui=tau(i);
   if (overlap==0)
       dt=taui;
   end;
   
   x2=x(2*taui+1:dt:length(x));
   x1=x(  taui+1:dt:length(x));
   x0=x(       1:dt:length(x));
   
   n = length(x0);
   if (length(x1) < n)
       n=length(x1);
   end;
   if (length(x2) < n)
       n=length(x2);
   end;
   
   if (n==0)
       display(['Not enough data for tau = ' num2str(taui) ': breaking']);
       break;
   end;
   
   varr = x2(1:n) - 2*x1(1:n) + x0(1:n);
   if (gaps == 1)
       varr = varr(~isnan(varr)); % remove NaNs
       n=length(varr); % new size
   end;
   avar = sum(varr .* varr);
   dev(i) = sqrt(avar/(2.0*n))/taui*rate; 
   deverr(i) = dev(i)/sqrt(n);
   ndev(i)=n;
   
end

end

