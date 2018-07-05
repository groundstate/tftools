function [ bias, totcorr ] = totdevbias(tau,T,noisefn,tot)
%TOTDEVBIAS Calculates the bias correction for the total
% deviation, given the noise type, and optionally corrects a vector
% of TOTDEV deviation values
%
%   Usage: [ bias, totcorr ] = TOTDEVBIAS(tau,T,noisefn,tot) 
%   Inputs:
%     
%     tau:  vector of averaging intervals
%				T:  time interval spanned by the data set
% noisefn: "white pm","flicker pm","white fm","flicker fm",
%             "random walk fm" (nb biased only for flicker fm and 
%							random walk fm)
%   tot: optional: data to be corrected
%        
%   Outputs:
%     bias:  vector of bias corrections, defined as
%            bias(tau)= ADEV(tau)/TOTDEV(tau) - 1
% totcorr: corrected TOTDEV
%
%   See also totdev.

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

% 
% Ref: W.J. Riley NIST SP 1065 "Handbook of Freequency Stability Analysis" p.49
% 

bias = zeros(size(tau));
totcorr = ones(size(tau));

if (strcmpi(noisefn,'white pm'))
  % no bias 
elseif (strcmpi(noisefn,'flicker pm'))
  % no bias
elseif (strcmpi(noisefn,'white fm'))
  % no bias
elseif (strcmpi(noisefn,'flicker fm'))
   bias = 1.0 ./ sqrt(1 - 0.481 * tau/T) - 1;
elseif (strcmpi(noisefn,'random walk fm'))
   bias = 1.0 ./ sqrt(1.0 - 0.75 * tau/T) - 1;
else
   error('tftools:totdevbias:BadInput', ['unknown noise type: ' noisefn]);
end

if (nargin == 3)
    totcorr = tot .* (bias+1);
end

end

