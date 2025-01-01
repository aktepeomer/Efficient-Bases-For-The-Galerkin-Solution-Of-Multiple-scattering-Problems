function [ exact_soln ] = Exact_solution_for_one_reflection( input,kernels,rhs,reflection_no,path )
%input = {1=point_number=n, 2=points=t, 3=wavenumber=eta, 4=basis_dimension=m, 5=configuration_type, 6=no_of_intervals, 7=chan_of_var, 8=kerneltype, 9=no_of_reflection}
%path=1 first path, path=2 second path

% Nyström soln of each reflection is find here

n=input{1};

obstacle = mod(path+reflection_no, 2) +1;
% finds on which obstacle we are taking integral.

A = eye(2*n)-kernels{obstacle,obstacle};
b = rhs;
exact_soln = linsolve( A , b );


exact_soln=exact_soln.'; 
% [1,2*n]

end

