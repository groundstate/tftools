function [ newtau,mtau ] = tau2m( tau,rate,x )
%TAU2M Converts tau to index, given the rate & sanitises input
%   Usage: [newtau,mtau] = taus2m(tau,rate,x)
%   tau is the vector of tau values
%   rate is data rate in Hz
%   x is the data
%   newtau is vector of sanitised tau values
%   mtau  is newtau in units of data points
%
%   See also adev, freq2phase, hdev, markgaps, mdev, tdev, totdev.
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

% A few inanity checks
if (rate == 0)
    error('tftools:tau2m:BadInput', 'the rate must be non-zero');
end;

% tau should be positive and not exceed the time spanned by the data set
newtau=tau;
newtau(~(tau>0))=[];
newtau(~(newtau < length(x)/rate))=[];

if (isempty(newtau))
    error('tftools:tau2m:BadInput', 'sanity check on tau failed');
end;

% 
mtau=round(newtau*rate);
mtau=unique(mtau); 

newtau=mtau/rate;

end

