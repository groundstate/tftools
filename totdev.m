function [ totdv,totdverr, ntotdv ] = totdev( x,rate,tau,phase,gaps )
%TOTDEV Calculate total deviation of phase or fractional frequency data
% 
%   d is the input time series
%   rate is the sampling rate in Hz
%   tau is the averaging interval 
%   phase = 1 means data is phase
%   gaps = 1 means data may contain gaps
%
%  Frequency data is converted to phase
%  Maximum tau allowed is N - 1
% 
%  See also adev, freq2phase, markgaps, mdev
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

if (phase == 0)
end;

ntau=length(tau);
N  = length(x);
% FIXME sanity check on tau

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
totdv=zeros(1,ntau);
totdverr = zeros(1,ntau);
ntotdv=zeros(1,ntau);

for i=1:ntau 
   taui=tau(i);
      
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
       display(['Not enough data for tau = ' num2str(taui) ': breaking']);
       break;
   end;
   
   varr = xkminusm(1:n) - 2*xk(1:n) + xkplusm(1:n);
   if (gaps == 1)
       varr = varr(~isnan(varr)); % remove NaNs
       n=length(varr); % new size
   end;
   
   tvar = sum(varr .* varr);
   totdv(i) = sqrt(tvar/(2.0*n))/taui*rate; 
   totdverr(i) = totdv(i)/sqrt(n);
   ntotdv(i)=n;
end

end