function [ dev,deverr, ndev, new_tau ] = theo1dev( x,rate,tau,phase,showprogress)
% THEO1DEV Calculates the Theo1 deviation
% THEO1DEV is able to 
%   Usage: [ dev, deverr, ndev, new_tau ] = THEO1DEV(x,rate,tau,overlap,phase,gaps) 
%   Inputs:
%     x       input time series
%     rate    sampling rate, in Hz
%     tau     averaging intervals (should be 0.75*(m-1)/rate, m even)
%     phase   data is phase (=1), optional argument, default=1
%     showprogress progress indicator, since theo1dev can be slow
%   Outputs:
%     dev     Theo1 deviations 
%     deverr  uncertainties (simple estimate)
%     ndev    number of samples used
%     new_tau tau values that were used
%   
%   See also adev,freq2phase, hdev, mdev, tau2m, tdev, totdev.

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

if (nargin > 5)
    error('tftools:theo1dev:TooManyInputs', 'requires 5 arguments');
end;

defaults = {1 0};

switch nargin
    case 3
        [phase,showprogress] = defaults{:};
    case 4
        showprogress = defaults{2};
end;

if (phase == 0)
    x=freq2phase(x,rate);
end;

% Theo1 requires that m be even so input sanitization has to be done here
% Assume that the input is correct
Nx=length(x);

new_tau=tau;
% mtau=round(1+new_tau*rate/0.75);
mtau=round(new_tau*rate/0.75);
mtau(mod(mtau,2) > 0) = []; % m must be even
mtau(mtau < 10) = []; % m must be >= 10 (removes negative m too)
mtau(mtau > Nx -1) = []; 
mtau=unique(mtau); % this sorts as well

if (isempty(mtau))
    error('tftools:theo1dev:BadInput', 'sanity check on tau failed');
end;

%new_tau=0.75*(mtau-1)/rate;
new_tau=0.75*mtau/rate;

ntau=length(new_tau);
dev=zeros(1,ntau);
deverr = zeros(1,ntau);
ndev=zeros(1,ntau);
if (showprogress ~= 0)
    fprintf('Calculating theo1 ...\n');
    fprintf('m=');
end
for t=1:ntau
    m=mtau(t);
    si=0;
    if (showprogress ~= 0)
        fprintf(' %i',m);
    end
    for i=1:Nx-m
        d=0:m/2 - 1;
        x1=x(i+m);
        x2=x(i+m/2:1:i+m-1); % length = (i+m-1) - (i+m/2) + 1 = m/2;
        x3=x(i+m/2:-1:i+1);  % length = (i+m/2) - (i+1) + 1 = m/2
        x4=x(i); 
        xi = x1 - x2 - x3 + x4;
        xi2 = xi .* xi;
        si = si + sum( xi2 ./ (m/2- d));
    end
    var= si/(0.75*(Nx-m)*(m*m/(rate*rate)));
    dev(t)=sqrt(var);
    ndev(t)=(Nx -m );
    deverr(t)=dev(t)/sqrt(ndev(t));
end

if (showprogress ~= 0)
        fprintf('\n');
end

end


% Test data from D A Howe "TheoH: a hybrid, high-confidence statistic
% that improves on the Allan deviation" Metrologia v43 S322-S331 (2006)
% x = [ -2.15 -0.99 1 2.5 0.65 -3.71 -3.3 1.08 0.5 2.2 4.68 3.29 ];
% Nx=12, m=10
% Theo1(10,1,12) = 0.6623
