function [ sb ] = order_shadow_boundaries( input, SB , phase)
%input={1=point_number=n, 2=points=t, 3=wavenumber=eta, 4=basis_dimension=m, 5=configuration_type, 6=no_of_intervals, 7=chan_of_var, 8=kerneltype, 9=no_of_reflection}

% For SB:
% First row is for shape one.
% First column is angle. Second is x-axis, third is y-axis os shadow bndries.

no_of_reflection = input{9};
n = input{1};
sb = cell(1,no_of_reflection); % shadw boundaries for all reflections


for i=1:no_of_reflection
    
    % find which one has small angle. and set them in 0--2pi
    % a1 is first in counterclockwise parametrization and
    % a2 is second. But we dont know where the illuminated region lays.
    a1 = min( mod(SB{i}(1,1),2*pi) , mod(SB{i}(2,1),2*pi) );
    a2 = max( mod(SB{i}(1,1),2*pi) , mod(SB{i}(2,1),2*pi) );
    
    % find the physical location of the shadow boundaries.
    if a1 == SB{i}(1,1)
        a1place = [SB{i}(1,2),SB{i}(1,3)];
        a2place = [SB{i}(2,2),SB{i}(2,3)];
    else
        a1place = [SB{i}(2,2),SB{i}(2,3)];
        a2place = [SB{i}(1,2),SB{i}(1,3)];
    end
    
    % find interpolation point numbers of 2 mid points of shadow ponts
    % it is between 0 and 2*n-1
    mid1_no = floor( (a1+a2)/2 *n/pi );
    mid2_no = mod(mid1_no + n,2*n);
    % if mid points are equal to 0 the code gives error. so we just move
    % mid points little bit. These points will be used for index numbers of
    % some vectors and 0 is not considered to be as an index number in
    % matrices (in matlab).
    if mid1_no==0
        mid1_no=1;
    end
    if mid2_no==0
        mid2_no=1;
    end
    
    % The mid point with bigger phase fnc value should be in the shadow
    % region. using this inf. we list the shadow boundary points wrt
    % counterclockwise parametrization.
    
    % below "mid_no" is in illuminated region
    if phase{i}(mid1_no) < phase{i}(mid2_no)
        mid_no = mid1_no;
    else
        mid_no = mid2_no;
    end
    
    % below conditions on "mid_no" helps us to find which shadow boundary
    % point comes first in counterclockwise parametrization.
    if mid_no == mid1_no
        sb {i} = [ a1, a1place(1), a1place(2); a2, a2place(1), a2place(2)];
    elseif mid_no == mid2_no
        sb {i} = [ a2, a2place(1), a2place(2); a1, a1place(1), a1place(2)];
    end 
end
% now sb{i} is 3*3 matrix such that first row is first shadow bndry point
% in the parametrization. First is angle, and second and third one is x and
% y coordinates.


end

