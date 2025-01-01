function [norm]=twonorm(x)

% It calculates 2-norm for each column.
% It returns [1,2*n].
% This funtion is for writing the other functions in a nicer way.
% x must be [2,2*n] and 2*n is the number of interpolation points.


norm=sqrt(sum(abs(x).^2,1));




end