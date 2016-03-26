function [ dev, deverr, ndev, new_tau ] = hdev( x,rate,tau,overlap,phase, gaps )
%ADEV Calculate non-overlapping/overlapping Hadamard deviation of phase/
% or fractional frequency data
%
%   Usage: [ dev, deverr, ndev, new_tau ] = HDEV(x,rate,tau,overlap,phase,gaps) 
%   Inputs:
%     x       input time series
%     rate    sampling rate, in Hz
%     tau     averaging intervals 
%     overlap compute overlapping hdev (=1), optional argument, default=1 
%     phase   data is phase (=1), optional argument, default=1
%     gaps    data contains gaps, tagged with NaN, optional argument, default=0
%             Only works with phase data.
%   Outputs:
%     dev     Hadamard deviations 
%     deverr  uncertainties
%     ndev    number of samples used
%     new_tau tau values that were used
%   
%   See also adev, freq2phase, mdev, tau2m, tdev, totdev.

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
% 
%

if (nargin > 6)
    error('tftools:hdev:TooManyInputs', 'requires at most 3 optional arguments');
end;

if (nargin == 6)
    if (phase == 0 && gaps == 1)
        error('tftools:hdev:BadInput', 'gaps in frequency data not supported');
    end;
end;

defaults = {1 1 0};

switch nargin
    case 3
        [overlap, phase, gaps] = defaults{:};
    case 4
        [phase, gaps]=defaults{2:3};
    case 5
        gaps=defaults{3};
end;

if (phase == 0)
    x=freq2phase(x,rate);
end;

% Validate tau etc
[ new_tau,mtau ] = tau2m( tau,rate,x );

ntau=length(mtau);

dev=zeros(1,ntau);
deverr = zeros(1,ntau);
ndev=zeros(1,ntau);

dt = 1; % for overlapping Hadamard deviation

for i=1:ntau 
   taui=mtau(i);
   if (overlap==0)
       dt=taui;
   end;
   
   x3=x(3*taui+1:dt:length(x));
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
   if (length(x3) < n)
       n=length(x3);
   end;
   if (n==0)
       display(['hdev: not enough data for tau = ' num2str(new_tau(i)) ': breaking']);
       break;
   end;
   
   varr = x3(1:n) - 3*x2(1:n) + 3*x1(1:n)- x0(1:n);
   if (gaps == 1)
       varr = varr(~isnan(varr)); % remove NaNs
       n=length(varr); % new size
   end;
   
   hvar = sum(varr .* varr);
   dev(i) = sqrt(hvar/(6.0*n))/taui*rate; 
   deverr(i) = dev(i)/sqrt(n);
   ndev(i)=n;
   
end

end


