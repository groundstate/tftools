function [ dev,deverr,ndev,new_tau ] = mdev( x,rate,tau,phase,gaps)
%MDEV Calculate modified Allan deviation of phase/
% or fractional frequency data
%   Usage: [ dev,deverr,ndev,new_tau ] = MDEV(x,rate,tau,phase,gaps) 
%   Inputs:
%     x       input time series
%     rate    sampling rate, in Hz
%     tau     averaging intervals 
%     phase   data is phase (=1), optional argument, default=1
%     gaps    data contains gaps, tagged with NaN, optional argument, default=0
%             Only works with phase data.
%   Outputs:
%     dev     modified Allan deviations 
%     deverr  uncertainties
%     ndev    number of samples used
%     new_tau tau values that were used
%
%   See also adev, freq2phase, markgaps, tau2m, tdev, totdev.
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
% This code is based on allantools.py, written by Anders Wallin
%

if (nargin > 5)
    error('tftools:mdev:TooManyInputs', 'requires at most 2 optional arguments');
end;

if (nargin==5) % FIXME
    if (gaps==1)
        error('tftools:mdev:BadInput', 'Data with gaps not supported yet');
    end;
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

% Validate tau etc
[ new_tau,mtau ] = tau2m( tau,rate,x );

N=length(x);
ntau=length(mtau);

dev=zeros(1,ntau);
deverr=zeros(1,ntau);
ndev=zeros(1,ntau);

% Someone is clever
% There's a way of doing this  by unrolling the loop
for i=1:ntau 
   taui=mtau(i);
  
   if (2*taui >= N || 3*taui >= N) 
       display(['Not enough data for tau = ' num2str(taui) ': breaking']);
       break
   end;
   
   % First summation
   x2 = x(2*taui+1:3*taui); 
   x1 = x(taui+1:2*taui);
   x0 = x(1:taui);
   
   n=taui;
   
   varr = sum(x2(1:n) - 2*x1(1:n) + x0(1:n));
   mvar=  varr * varr;
   
   % Second summation 
   x3 = x(3*taui+1:N); 
   x2 = x(2*taui+1:N-taui);
   x1 = x(  taui+1:N-2*taui);
   x0 = x(       1:N-3*taui);
   
   n = N-3*taui;
   
   varr = varr + cumsum(x3(1:n) - 3*x2(1:n) + 3*x1(1:n) - x0(1:n));
   mvar = mvar+ sum(varr .* varr);
   n = n+1; % count the first one too
  
   dev(i) = sqrt(mvar/(2.0*n))/(taui*taui)*rate;
   deverr(i) = dev(i)/sqrt(n);
   ndev(i)=n;
   
end

end

