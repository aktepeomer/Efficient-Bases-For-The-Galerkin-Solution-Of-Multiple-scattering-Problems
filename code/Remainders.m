function [ remainder_1, remainder_2 ] = Remainders( input, kernels, first_cur_1, first_cur_2, no_of_reflection )
%{1=point_number=n, 2=points=t, 3=wavenumber=eta, 4=basis_dimension=m, 5=configuration_type, 6=no_of_intervals, 7=chan_of_var, 8=kerneltype, 9=no_of_reflection}


n   = input{1};
eta = input{3};
m   = no_of_reflection;

T12 = ( eye(2*n)-kernels{1,1} ) \ kernels{1,2};
T21 = ( eye(2*n)-kernels{2,2} ) \ kernels{2,1};
g1  = ( eye(2*n)-kernels{1,1} ) \ (first_cur_1.');
g2  = ( eye(2*n)-kernels{2,2} ) \ (first_cur_2.');
B = T12 * T21;
A = B^m;

b1 = A   * B   * T12     * g2;
a1 = A   * B^2           * g1;
a2 = T21 * A   * B       * g1;
b2 = T21 * A   * B * T21 * g2;


% değişebilir
constant = 1 / ( 1 - 0.5*exp(2*complex(0,1)*eta) );


remainder_1 = (a1 +b1)*constant;
remainder_2 = (a2 +b2)*constant;

remainder_1 = remainder_1.';
remainder_2 = remainder_2.';

end

