function [acfx  ] = acf( x )
%ACF Calculate the autocorrelation function of a time series x(i),
%    defined as per eq.47 in NIST SP1065
%
%   Usage: [acfx  ] = ACF(x) 
%   Inputs:
%     x       input time series
%   Outputs:
%     acfx    the auto correlation function %   
%   See also freq2phase, hdev, mdev, tau2m, tdev, totdev.

% The MIT License (MIT)
% 
% Copyright (c) 2017 Michael J. Wouters
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
%

N = length(x);
xm = mean(x);
vm = sum((x-xm) .* (x-xm));
acfx = zeros(1,N);

for j=0:N-1
    acfx(j+1)=sum( ( x(1:N-j) - xm ) .* ( x(1+j:N) - xm ));
end

acfx = acfx/vm;

end

