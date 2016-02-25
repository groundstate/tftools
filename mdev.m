function [ mdv,mdverr, nmdv ] = adev( x,rate,tau,phase, gaps )
%ADEV Calculate non-overlapping/overlapping Allan deviation of phase/
% or fractional frequency data
%   Usage: adev(x,rate,tau,overlap,phase,gaps) 
%   x is the input time series
%   rate is the sampling rate, in Hz
%   tau is the averaging interval 
%   overlap = 1 means compute overlapping ADEV
%   phase = 1 means data is phase
%   gaps = 1 means data contains gaps, tagged with NaN
%
%   Frequency data is converted internally to phase
%
%   See also adev,freq2phase,markgaps,totdev.
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

N=length(x);
ntau=length(tau);
% FIXME sanity check on tau 

mdv=zeros(1,ntau);
mdverr=zeros(1,ntau);
nmdv=zeros(1,ntau);

% Someone is clever
% There's a way of doing this  by unrolling the loop
for i=1:ntau 
   taui=tau(i);
  
   % First summation
   x2 = x(2*taui+1:3*taui); 
   x1 = x(taui+1:2*taui);
   x0 = x(1:taui);
   
   n=taui;
   
   varr = sum(x2(1:n) - 2*x1(1:n) + x0(1:n));
   mvar=  varr * varr;
   
   x3 = x(3*taui+1:N); 
   x2 = x(2*taui+1:N-taui);
   x1 = x(  taui+1:N-2*taui);
   x0 = x(       1:N-3*taui);
   
   n = N-3*taui;
  
   varr = varr + cumsum(x3(1:n) - 3*x2(1:n) + 3*x1(1:n) - x0(1:n));
   mvar = mvar+ sum(varr .* varr);
   n = n+1; % count the first one too
  
   mdv(i) = sqrt(mvar/(2.0*n))/(taui*taui)*rate*rate; % FIXME
   mdverr(i) = mdv(i)/sqrt(n);
   nmdv(i)=n;
   
end

end

