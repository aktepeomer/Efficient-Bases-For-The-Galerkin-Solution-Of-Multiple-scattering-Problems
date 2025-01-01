function [ z,ksi,validity ] = intervals_with_optimization( t1,t2,noofintervals,k,shapetype,chan_of_var,ksi_changes )
% This function returns a [noofintervals,4] matrix of intervals whose jth
% column means: (j,1)=a, (j,2)=a' (j,3)=b', (j,4)=b for jth interval.

% This fnc determines intervals for galerkin method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HERE WE ALWAYS USE CIRCLE CASE.
shapetype = 1;
% CHANGES FOR CHANGE OF VARIABLES CASE IS NOT DONE YET. TO CREATE A GOOD
% METHOD WHICH WORKS FOR BOTH CIRCLE AND ELLIPSE I NEED TO TRY A COUPLE OF
% THINGS. AFTER SMALL MODIFICATIONS WE CAN HAVE A SHORTER AND BETTER WAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% t1 is before t2 in counterclockwise parametrization. in [0,2pi]
% k is wave number


L=2*pi;

if chan_of_var==0 % If we do not have change of variables on Transition regions.
    
    switch shapetype
        case 1 % For circle
            
            switch noofintervals
                % case j means j many intervals are used
                % epsilon(j) and epsilon(k) are same as Ecevit&Çaðan section 5.1
                case 3
                    epsilon=1/9;
                    ksi=[L/10,7*L/10,7*L/30];
                case 4
                    epsilon=[1/9,1/9,1/9];
                    ksi=[L/10,L/3,L/2,7*L/30];
                case 5
                    epsilon=[1/9,1/9,1/9];
                    ksi=[L/10,L/3,L/2,7*L/30];
                case 6
                    epsilon=[1/5,1/15,1/15];
                    ksi=[L/4,3*L/10,7*L/15,7*L/30,L/15,7*L/20];
                case 7
                    epsilon=[1/5,1/15,1/15];
                    ksi=[L/4,3*L/10,7*L/15,7*L/30,L/15,7*L/20];
                case 8
                    epsilon=[1/5,1/15,1/5];
                    ksi=[L/4,3*L/10,7*L/15,7*L/30,L/15,7*L/20,2*L/5,3*L/10;L/4,3*L/10,7*L/15,7*L/30,L/15,7*L/20,2*L/5,3*L/10];
                    %ksi=[L/4,3*L/10,7*L/15,7*L/30] yazýyor hh'de fark ne bilmiyorum
            end
            
        case 2 % For ellipse
            switch noofintervals
                % case j means j many intervals are used
                % epsilon(j) and epsilon(k) are same as Ecevit&Çaðan section 5.1
                case 3
                    epsilon=1/9;
                    ksi=[1/10,1/3,1/4]*L;
                case 4
                    epsilon=[1/9,1/9,1/9];
                    ksi=[1/10,3/20,3/10,1/4]*L;
                case 5
                    epsilon=[1/9,1/9,1/9];
                    ksi=[1/10,3/20,3/10,1/4]*L;
                case 6
                    epsilon=[2/9,1/9,1/9];
                    ksi=[3/20,1/5,1/4,3/20,1/20,1/5]*L;
                case 7
                    epsilon=[2/9,1/9,1/9];
                    ksi=[3/20,1/5,1/4,3/20,1/20,1/5]*L;
                case 8
                    epsilon=[2/9,1/9,2/9];
                    ksi=[3/20,1/5,1/4,3/20,1/20,1/5,7/30,1/5]*L;
                    % HH inki farklý olabilir.
            end
            
    end
elseif chan_of_var==1 % if there is change of variables
    
    switch shapetype
        case 1 % For circle
            % If we do have change of variables on Transition regions.
            % So this works when there are 6 or 5 regions.
            switch noofintervals
                % case j means j many intervals are used
                case 5
                    %epsilon=[1/3,0,1/3];
                    epsilon=[1/3,0,1/3];
                    %ksi=[0, 0.3, 5.25, 3.70, 1.70, 2.0, 0, 0];
                    ksi=[0, 0.3, 5.25, 3.70, 1.70, 0, 0, 0];
                case 6
                    %epsilon=[1/3,0,1/3];
                    epsilon=[1/3,0,1/3];
                    %ksi=[0.7,1.0,5.15,3.65,2.10,2.0,1.8,1.8];
                    %ksi=[0, 1.0, 5.15, 3.65, 2.10, 2.0, 0, 0];
                    ksi=[0, 1.0, 5.05, 3.75, 2, 0, 0, 0];
            end
        case 2 % For ellipse
            switch noofintervals
                % case j means j many intervals are used
                case 6
                    epsilon=[1/3,0,1/3];
                    %ksi=[0.7,1.0,5.15,3.65,2.10,2.0,1.8,1.8];
                    %ksi=[0, 1.05, 5.75, 3.45, 2.30, 0, 0, 0;...
                    %    0, 0.8, 5.25, 3.35, 2.30, 0, 0, 0];
                    ksi=[0, 1.05, 5.75, 3.45, 2.0, 0, 0, 0; 0, 0.8, 5.25, 3.35, 2.0, 0, 0, 0];
            end
    end
    
end

% Optimization part. we change ksi values
ksi=ksi+ksi_changes*L;

y=zeros(noofintervals,4);
% Intervals below are listed counterclockwise. So, one region , SB for
% exmp, can be in a different place in two different region cases.
% y(j,1) is "a" of jth interval
% y(j,2) is "a'" of jth interval
% y(j,3) is "b'" of jth interval
% y(j,4) is "b" of jth interval

% they are now in [0,2pi]
t1=mod(t1,2*pi);
t2=mod(t2,2*pi);

% We use these inf. on interval construction and scaling (in reflection)
if t1<t2
    relativezero=mod((t1+t2)/2+pi,2*pi);
    relative_length=mod(t2-t1,2*pi)/pi;
elseif t1>t2
    relativezero=mod((t1+t2)/2,2*pi);
    relative_length=mod(t2-t1,2*pi)/pi;
end

% So, y is the list of all "a,a',b',b"s (from top to bottom)
if chan_of_var==0
    % If we do not have change of variables on Transition regions.
    switch noofintervals
        case 3
            %shadow boundary 1
            y(1,1)=t1-ksi(2)*k^(-1/3);
            y(1,4)=t1+ksi(3)*k^(epsilon-1/3);
            
            %illuminated region
            y(2,1)=t1+ksi(1)*k^(epsilon-1/3);
            y(2,4)=t2-ksi(1)*k^(epsilon-1/3);
            
            %shadow boundary 2
            y(3,1)=t2-ksi(3)*k^(epsilon-1/3);
            y(3,4)=t2+ksi(2)*k^(-1/3);
            
        case 4
            %deep shadow region
            y(1,1)=t2+ (ksi(2)*k^(epsilon(3)-1/3))*(2-relative_length);
            y(1,4)=t1- (ksi(2)*k^(epsilon(3)-1/3))*(2-relative_length);
            
            %shadow boundary 1
            y(2,1)=t1- (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
            y(2,4)=t1+ (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
            
            %illuminated region
            y(3,1)=t1+ (ksi(1)*k^(epsilon(1)-1/3))*relative_length;
            y(3,4)=t2- (ksi(1)*k^(epsilon(1)-1/3))*relative_length;
            
            %shadow boundary 2
            y(4,1)=t2- (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
            y(4,4)=t2+ (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
            
        case 5
            switch shapetype
                case 1 % For circle
                    %deep shadow region 1
                    y(1,1)=relativezero-L/32;
                    y(1,4)=t1- (ksi(2)*k^(epsilon(3)-1/3))*(2-relative_length);
                    
                    %shadow boundary 1
                    y(2,1)=t1- (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
                    y(2,4)=t1+ (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
                    
                    %illuminated region
                    y(3,1)=t1+ (ksi(1)*k^(epsilon(1)-1/3))*relative_length;
                    y(3,4)=t2- (ksi(1)*k^(epsilon(1)-1/3))*relative_length;
                    
                    %shadow boundary 2
                    y(4,1)=t2- (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
                    y(4,4)=t2+ (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
                    
                    %deep shadow region 2
                    y(5,1)=t2+ (ksi(2)*k^(epsilon(3)-1/3))*(2-relative_length);
                    y(5,4)=relativezero+L/32;
                    
                case 2 % For ellipse
                    %deep shadow region 1
                    y(1,1)=relativezero-L/32;
                    y(1,4)=t1- (ksi(2)*k^(epsilon(3)-1/3))*(2-relative_length);
                    
                    %shadow boundary 1
                    y(2,1)=t1- (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
                    y(2,4)=t1+ (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
                    
                    %illuminated region
                    y(3,1)=t1+ (ksi(1)*k^(epsilon(1)-1/3))*relative_length;
                    y(3,4)=t2- (ksi(1)*k^(epsilon(1)-1/3))*relative_length;
                    
                    %shadow boundary 2
                    y(4,1)=t2- (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
                    y(4,4)=t2+ (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
                    
                    %deep shadow region 2
                    y(5,1)=t2+ (ksi(2)*k^(epsilon(3)-1/3))*(2-relative_length);
                    y(5,4)=relativezero+L/32;
            end
            
        case 6
            %deep shadow region
            y(1,1)=t2+ (ksi(2)*k^(epsilon(3)-1/3))*(2-relative_length);
            y(1,4)=t1- (ksi(2)*k^(epsilon(3)-1/3))*(2-relative_length);
            
            %shadow boundary 1
            y(2,1)=t1- (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
            y(2,4)=t1+ (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
            
            %illuminated transition 1
            y(3,1)=t1+ (ksi(5)*k^(epsilon(2)-1/3))*relative_length;
            y(3,4)=t1+ (ksi(6)*k^(epsilon(1)-1/3))*relative_length;
            
            %illuminated region
            y(4,1)=t1+ (ksi(1)*k^(epsilon(1)-1/3))*relative_length;
            y(4,4)=t2- (ksi(1)*k^(epsilon(1)-1/3))*relative_length;
            
            %illuminated transition 2
            y(5,1)=t2- (ksi(6)*k^(epsilon(1)-1/3))*relative_length;
            y(5,4)=t2- (ksi(5)*k^(epsilon(2)-1/3))*relative_length;
            
            %shadow boundary region 2
            y(6,1)=t2- (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
            y(6,4)=t2+ (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
            
        case 7
            switch shapetype
                case 1 % For circle
                    %deep shadow region 1
                    y(1,1)=relativezero-L/32;
                    y(1,4)=t1- (ksi(2)*k^(epsilon(2)-1/3))*(2-relative_length);
                    
                    %shadow boundary 1
                    y(2,1)=t1- (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
                    y(2,4)=t1+ (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
                    
                    %illuminated transition 1
                    y(3,1)=t1+ (ksi(5)*k^(epsilon(2)-1/3))*relative_length;
                    y(3,4)=t1+ (ksi(6)*k^(epsilon(1)-1/3))*relative_length;
                    
                    %illuminated region
                    y(4,1)=t1+ (ksi(1)*k^(epsilon(1)-1/3))*relative_length;
                    y(4,4)=t2- (ksi(1)*k^(epsilon(1)-1/3))*relative_length;
                    
                    %illuminated transition 2
                    y(5,1)=t2- (ksi(6)*k^(epsilon(1)-1/3))*relative_length;
                    y(5,4)=t2- (ksi(5)*k^(epsilon(2)-1/3))*relative_length;
                    
                    %shadow boundary region 2
                    y(6,1)=t2- (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
                    y(6,4)=t2+ (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
                    
                    %deep shadow region 2
                    y(7,1)=t2+ (ksi(2)*k^(epsilon(2)-1/3))*(2-relative_length);
                    y(7,4)=relativezero+L/32;
                    
                case 2 % For ellipse
                    %deep shadow region 1
                    y(1,1)=relativezero-L/32;
                    y(1,4)=t1- (ksi(2)*k^(epsilon(2)-1/3))*(2-relative_length);
                    
                    %shadow boundary 1
                    y(2,1)=t1- (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
                    y(2,4)=t1+ (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
                    
                    %illuminated transition 1
                    y(3,1)=t1+ (ksi(5)*k^(epsilon(2)-1/3))*relative_length;
                    y(3,4)=t1+ (ksi(6)*k^(epsilon(1)-1/3))*relative_length;
                    
                    %illuminated region
                    y(4,1)=t1+ (ksi(1)*k^(epsilon(1)-1/3))*relative_length;
                    y(4,4)=t2- (ksi(1)*k^(epsilon(1)-1/3))*relative_length;
                    
                    %illuminated transition 2
                    y(5,1)=t2- (ksi(6)*k^(epsilon(1)-1/3))*relative_length;
                    y(5,4)=t2- (ksi(5)*k^(epsilon(2)-1/3))*relative_length;
                    
                    %shadow boundary region 2
                    y(6,1)=t2- (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
                    y(6,4)=t2+ (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
                    
                    %deep shadow region 2
                    y(7,1)=t2+ (ksi(2)*k^(epsilon(2)-1/3))*(2-relative_length);
                    y(7,4)=relativezero+L/32;
            end
            
        case 8
            %deep shadow region
            y(1,1)=t2+ (ksi(2,2)*k^(epsilon(3)-1/3))*(2-relative_length);
            y(1,4)=t1- (ksi(1,2)*k^(epsilon(3)-1/3))*(2-relative_length);
            
            %shadow transition region 1
            y(2,1)=t1- (ksi(1,7)*k^(epsilon(1)-1/3))*(2-relative_length);
            % epsilon(3) olmasý lazým olabilir!!!
            y(2,4)=t1- (ksi(1,8)*k^(epsilon(2)-1/3))*(2-relative_length);
            
            %shadow boundary 1
            y(3,1)=t1- (ksi(1,3)*k^(epsilon(2)-1/3))*(2-relative_length);
            y(3,4)=t1+ (ksi(1,4)*k^(epsilon(2)-1/3))*relative_length;
            
            %illuminated transition 1
            y(4,1)=t1+ (ksi(1,5)*k^(epsilon(2)-1/3))*relative_length;
            y(4,4)=t1+ (ksi(1,6)*k^(epsilon(1)-1/3))*relative_length;
            
            %illuminated region
            y(5,1)=t1+ (ksi(1,1)*k^(epsilon(1)-1/3))*relative_length;
            y(5,4)=t2- (ksi(2,1)*k^(epsilon(1)-1/3))*relative_length;
            
            %illuminated transition 2
            y(6,1)=t2- (ksi(2,6)*k^(epsilon(1)-1/3))*relative_length;
            y(6,4)=t2- (ksi(2,5)*k^(epsilon(2)-1/3))*relative_length;
            
            %shadow boundary 2
            y(7,1)=t2- (ksi(2,4)*k^(epsilon(2)-1/3))*relative_length;
            y(7,4)=t2+ (ksi(2,3)*k^(epsilon(2)-1/3))*(2-relative_length);
            
            %shadow transition region 2
            y(8,1)=t2+ (ksi(2,8)*k^(epsilon(2)-1/3))*(2-relative_length);
            y(8,4)=t2+ (ksi(2,7)*k^(epsilon(1)-1/3))*(2-relative_length);
    end
    
elseif chan_of_var==1
    % If we do have change of variables on Transition regions.
    % So this works when there are more than 5 regions.
    switch noofintervals
        case 5
            switch shapetype
                case 1 % For circle
                    %deep shadow region
                    %y(1,1)=t2+ksi(2)*k^(epsilon(3)-1/3);
                    y(1,1)=t2+ (ksi(2)*k^(epsilon(2)-1/3))*(2-relative_length);
                    %y(1,4)=t1-ksi(2)*k^(epsilon(3)-1/3);
                    y(1,4)=t1- (ksi(2)*k^(epsilon(2)-1/3))*(2-relative_length);
                    
                    %shadow boundary 1
                    y(2,1)=t1- (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
                    y(2,4)=t1+ (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
                    
                    %illuminated transition 1
                    y(3,1)=t1+ (ksi(5)*k^(epsilon(2)-1/3))*relative_length;
                    %y(3,4)=t1+ksi(6)*k^(epsilon(1)-1/3);
                    y(3,4)=t1+ (ksi(5)*k^(epsilon(1)-1/3))*relative_length;
                    
                    %illuminated transition 2
                    %y(4,1)=t2-ksi(6)*k^(epsilon(1)-1/3);
                    y(4,1)=t2- (ksi(5)*k^(epsilon(1)-1/3))*relative_length;
                    y(4,4)=t2- (ksi(5)*k^(epsilon(2)-1/3))*relative_length;
                    
                    %shadow boundary 2
                    y(5,1)=t2- (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
                    y(5,4)=t2+ (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
                    
                case 2 % For ellipse
                    
                    
            end
            
        case 6
            switch shapetype
                case 1 % For circle
                    %deep shadow region 1
                    y(1,1)=relativezero-L/32;
                    %y(1,4)=t1-ksi(2)*k^(epsilon(3)-1/3);
                    y(1,4)=t1- (ksi(2)*k^(epsilon(2)-1/3))*(2-relative_length);
                    
                    %shadow boundary 1
                    y(2,1)=t1- (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
                    y(2,4)=t1+ (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
                    
                    %illuminated transition 1
                    y(3,1)=t1+ (ksi(5)*k^(epsilon(2)-1/3))*relative_length;
                    %y(3,4)=t1+ksi(6)*k^(epsilon(1)-1/3);
                    y(3,4)=t1+ (ksi(5)*k^(epsilon(1)-1/3))*relative_length;
                    
                    %illuminated transition 2
                    %y(4,1)=t2-ksi(6)*k^(epsilon(1)-1/3);
                    y(4,1)=t2- (ksi(5)*k^(epsilon(1)-1/3))*relative_length;
                    y(4,4)=t2- (ksi(5)*k^(epsilon(2)-1/3))*relative_length;
                    
                    %shadow boundary region 2
                    y(5,1)=t2- (ksi(4)*k^(epsilon(2)-1/3))*relative_length;
                    y(5,4)=t2+ (ksi(3)*k^(epsilon(2)-1/3))*(2-relative_length);
                    
                    %deep shadow region 2
                    %y(6,1)=t2+ksi(2)*k^(epsilon(3)-1/3);
                    y(6,1)=t2+ (ksi(2)*k^(epsilon(2)-1/3))*(2-relative_length);
                    y(6,4)=relativezero+L/32;
                    
                case 2 % For ellipse
                    %deep shadow region 1
                    y(1,1)=relativezero-L/32;
                    %y(1,4)=t1-ksi(2)*k^(epsilon(3)-1/3);
                    y(1,4)=t1- (ksi(1,2)*k^(epsilon(2)-1/3))*(2-relative_length);
                    
                    %shadow boundary 1
                    y(2,1)=t1- (ksi(1,3)*k^(epsilon(2)-1/3))*(2-relative_length);
                    y(2,4)=t1+ (ksi(1,4)*k^(epsilon(2)-1/3))*relative_length;
                    
                    %illuminated transition 1
                    y(3,1)=t1+ (ksi(1,5)*k^(epsilon(2)-1/3))*relative_length;
                    %y(3,4)=t1+ksi(6)*k^(epsilon(1)-1/3);
                    y(3,4)=t1+ (ksi(1,5)*k^(epsilon(1)-1/3))*relative_length;
                    
                    %illuminated transition 2
                    %y(4,1)=t2-ksi(6)*k^(epsilon(1)-1/3);
                    y(4,1)=t2- (ksi(2,5)*k^(epsilon(1)-1/3))*relative_length;
                    y(4,4)=t2- (ksi(2,5)*k^(epsilon(2)-1/3))*relative_length;
                    
                    %shadow boundary region 2
                    y(5,1)=t2- (ksi(2,4)*k^(epsilon(2)-1/3))*relative_length;
                    y(5,4)=t2+ (ksi(2,3)*k^(epsilon(2)-1/3))*(2-relative_length);
                    
                    %deep shadow region 2
                    %y(6,1)=t2+ksi(2)*k^(epsilon(3)-1/3);
                    y(6,1)=t2+ (ksi(2,2)*k^(epsilon(2)-1/3))*(2-relative_length);
                    y(6,4)=relativezero+L/32;
            end
    end
    
end


% This part sets a' and b' values. [a,a',b',b]
for i=1:noofintervals
    y(i,2)=y(mod(i-1-1,noofintervals)+1,4);
    y(i,3)=y(mod(i+1-1,noofintervals)+1,1);
end

z=mod(y,2*pi); % Whole intervals' points are in between 0 and 2*pi.

% If a<a'<b'<b is not satisfied, this part fixes the problem. We assume in
% this part that 0 or 2pi is between a and b, and there are no other
% problems such as : a in [-2*pi,0], a' in [0,2*pi] and b2 in[2*pi,4*pi]".
for i=1:noofintervals
    for j=1:3
        if z(i,4-j)>z(i,4-j+1)
            z(i,1:4-j)=z(i,1:4-j)-L;
        end
    end
end
% Now [a,a',b',b] is in incresing order given that ksi changes are not
% problematic. If they are we deal with this below.

% If there was an unwanted situation about intervals this parts eliminates
% that case and dont continue calculations in the "Reflection solver" fnc.
for i=1:noofintervals
    if z(i,4)<z(i,3) || z(i,3)<z(i,2) || z(i,2)<z(i,1) || abs( z(i,4)-z(i,1) ) > 2*pi
        validity=0;
    else
        validity=1;
    end
        
end
% colums are regions and  rows of z are region information. 

end

