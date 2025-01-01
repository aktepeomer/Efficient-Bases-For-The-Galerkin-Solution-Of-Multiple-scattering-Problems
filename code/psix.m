function [ y ] = psix( x )
% This is the function "psi" from Ecevit÷zen page 13.
% It returns [1,2*n].

% For x<=0 :
y=x<=0;
% when x is less than 0 gives 1 other parts are 0
% It's lenght is same as x

% For 0<x<1 :
xzerone=x.*(x>0&x<1);
% when x is between 0 and 1 it gives the value of the entry, otherwise 0
y=y+ exp( 2*exp(-1./xzerone) ./ (xzerone-1) ).*(x>0&x<1);
% Last term gets rid off outcomes coming from cases where x<=0 and x>=0

% For 1<x :
y(x>=1)=0;
%this one gets rid off NaN (or inf) at x=1 (if it exists as grid point)


end

