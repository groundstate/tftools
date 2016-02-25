function [ x ] = freq2phase( y, rate )
%FREQ2PHASE Converts frequency data to phase data.
%   
%
%   See also adev, mdev, totdev.

x = cumsum(y)/rate;
x = [0 x];

end

