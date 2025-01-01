function [ rhs ] = FirstRHS( input,phase,configuration_info,obj )
%input = {1=point_number=n, 2=points=t, 3=wavenumber=eta, 4=basis_dimension=m, 5=configuration_type, 6=no_of_intervals, 7=chan_of_var, 8=kerneltype, 9=no_of_reflection}
%configuration_info = {1=X, 2=DX, 3=DDX, 4=shape_type, 5=shape_centers}
%obj=is object no

n = input{1};
eta = input{3};
kerneltype = input{8};


Phase = phase{1}; % 1 is incoming ray
psi = repmat([1;0],1,2*n);


dx = configuration_info{2}{obj} ;

% x is not used because we have phase function already.
% dx is the derivative of the shape wrt t.
% Dimensions are [2,2*n].

nu = [dx(2,:);-dx(1,:)]./[twonorm(dx);twonorm(dx)];
% nu is the outward directional normal vector of our shape.

dummy=(sum(psi.*nu,1)).';
% This is psi•nu. [2*n,1]

e_term=exp(complex(0,1)*eta*Phase);
e_term=e_term.';

% These are not from the theory but from Numerical Implementations last page.
switch kerneltype
    case 1 % 1st kind
        rhs=2*e_term;
    case 2 % 2nd kind
        rhs=2*complex(0,1)*eta*dummy.*e_term;
    case 3 % combined field
        rhs=2*complex(0,1)*eta*(dummy-1).*e_term;
end
% [2*n,1]



end

