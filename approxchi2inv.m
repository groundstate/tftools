function [ x ] = approxchi2inv( p, edf )
%APPROXCHI2INV Solves for chi squared, given the probability and 
% degrees of freedom
%   p    probability 
%   edf  estimated degrees of freedom

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
% Credits
% Based on the power series approximation described by Mark Weiss
% in "A simple algorithm for approximating confidence on the modified
% Allan variance and the time variance"
%

if (~isscalar(p) || ~isscalar(edf))
     error('tftools:approxchi2inv:BadInput','inputs must be scalars');
end

if (edf <=1)
     error('tftools:approxchi2inv:BadInput', 'edf must be greater than 1');
end

if (p < 0.005 || p > 0.995)
     error('tftools:approxchi2inv:BadInput', 'require 0.005 <= p <= 0.995');
end

if (p <= 0.5 && edf <= 10)
    a = edf/2.0;
    G = gamma(1.0+a);
    A=p*G;
    u=0;
    for j=1:7
        g= 1 + u/(a+1)*(1.0+u/(a+2)*(1+u/(a+3)));
        u=(A*exp(u)/g)^(1/a);
    end
    x=2*u;
else
    p1=min(p,1-p);
    t = sqrt(-2*log(p1));
    X = t - (2.30753 + 0.27601*t)/(1.0+t*(0.99229 + 0.04881*t));
    s=sign(p-0.5);
    b=2.0/(9*edf);
    x= edf*(1.0-b+s*X*sqrt(b))^3;
end

end

