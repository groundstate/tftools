function [ dev,deverr, ndev, new_tau ] = totdev( x,rate,tau,phase,gaps )
%TOTDEV Calculate total deviation of phase or fractional frequency data
% 
%   Usage: [ dev,deverr,ndev,new_tau ] = TOTDEV(x,rate,tau,phase,gaps) 
%   Inputs:
%     x       input time series
%     rate    sampling rate, in Hz
%     tau     averaging intervals 
%     phase   data is phase (=1), optional argument, default=1
%     gaps    data contains gaps, tagged with NaN, optional argument, default=0
%             Only works with phase data.
%   Outputs:
%     dev     total deviations 
%     deverr  uncertainties
%     ndev    number of samples used
%     new_tau tau values that were used
% 
%  See also adev, freq2phase, markgaps, tau2m, tdev, mdev
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
    error('tftools:totdev:TooManyInputs', 'requires at most 2 optional arguments');
end;

if (nargin == 5)
    if (phase == 0 && gaps == 1)
        error('tftools:adev:BadInput', 'gaps in frequency data not supported');
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

if (iscolumn(x)) % change to row vector
    x=x';
end;

% Validate tau etc
[ new_tau,mtau ] = tau2m( tau,rate,x );
ntau=length(mtau);
N  = length(x);

% Create the reflected data to be inserted before the original data
%  x*(1-j) = 2x(1) - x(1+j)  for j=1..N-2
x1 = 2.0*x(1)*ones(1,N-2); % first term
x1 = x1 - x(2:1:N-1);
x1 = x1(length(x1):-1:1); % reverse the vector

% Create the reflected data to be appended after the original data
%  x*(N+j) = 2x(N) - x(N-j)  for j=1..N-2
x2= 2.0*x(N)*ones(1,N-2);
x2= x2 - x(N-1:-1:2);

xe = [x1 x x2]; % the extended data set
% Sequence length is N-2 + N + N-2 = 3N-4
dev=zeros(1,ntau);
deverr = zeros(1,ntau);
ndev=zeros(1,ntau);

for i=1:ntau 
   taui=mtau(i);
      
   xkminusm=xe(N-1-(taui-1):1:N-1-(taui-1)+(N-3));
   xk=xe(N-2+2:N-2+N-1);
   xkplusm=xe(N-2+2+taui:1:N-2+N-1+taui);
   
   n = length(xk);
   if (length(xkminusm) < n)
       n=length(xkminusm);
   end;
   if (length(xkplusm) < n)
       n=length(xkplusm);
   end;
   
   if (n==0)
       display(['Not enough data for tau = ' num2str(new_tau(i)) ': breaking']);
       break;
   end;
   
   varr = xkminusm(1:n) - 2*xk(1:n) + xkplusm(1:n);
   if (gaps == 1)
       varr = varr(~isnan(varr)); % remove NaNs
       n=length(varr); % new size
   end;
   
   tvar = sum(varr .* varr);
   dev(i) = sqrt(tvar/(2.0*n))/taui*rate; 
   deverr(i) = dev(i)/sqrt(n);
   ndev(i)=n;
end

end