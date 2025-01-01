function [ Y ] = Basisfncs( points,n,m,chan_of_var,ksi,eta,no_of_intervals,i,t1,t2,shapetype )

% From left to right it runs over interpolation points.
% From top to buttom we have functions that look like x^0, x^1,...,x^m.
% This function returns a [m+1,2*n] matrix.

%points=[a,a',b',b]
%m:maximum polynomial degree
%n:there are 2n many points
%a<b<=c<d is known

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%My assumption here is "b-a<2*pi" and this it is trivial.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=points(1);
b=points(2);
c=points(3);
d=points(4);


x_n=-4*pi:pi/n:4*pi-pi/n;
dummy=x_n;
%x_n is [1,8*n], extended interpolation points
% they are defined on -4pi--4pi

% This part is for finding which ksi should be used in change of variables.
if chan_of_var==1 && (i==3||i==4) % if "Basisfncs" is called for a transition 
                                  % region find ksi values we are working on
                                  % so that we can apply change of
                                  % variables
    if shapetype==1 % For circle
        switch no_of_intervals
            case 5
                if i==3||i==4
                    ksim=ksi(5);
                end
            case 6
                if i==3||i==4
                    ksim=ksi(5);
                end
        end
    elseif shapetype==2 % For ellipse
        switch no_of_intervals
            case 5
                if i==3
                    ksim=ksi(1,5);
                elseif i==4
                    ksim=ksi(2,5);
                end
            case 6
                if i==3
                    ksim=ksi(1,5);
                elseif i==4
                    ksim=ksi(2,5);
                end
                
        end
    end
end
if chan_of_var==1 % Change of variable (only for Ýllnum. Trans. regions)
    if (no_of_intervals==6&&i==3)||(no_of_intervals==5&&i==3)
        % In this case our region is IT1
        x_n= log((x_n-t1)/ksim)*(3*(d-a)/log(eta)) +d;
    elseif (no_of_intervals==6&&i==4)||(no_of_intervals==5&&i==4)
        % In this case our region is IT2
        x_n= log((t2-x_n)/ksim)*(3*(a-d)/log(eta)) +a;
%     elseif (no_of_intervals==8&&i==2)
%         % In this case our region is ST1
%         x_n= (( log(x_n-t1) - log(-ksim) )/log(eta))*(-3*(d-a)) +a;
%     elseif (no_of_intervals==8&&i==8)
%         % In this case our region is ST2
%         x_n= (( log(x_n-t2) - log(ksim) )/log(eta))*(3*(d-a)) +d;
    end
    x_n(dummy<a)=dummy(dummy<a);
    x_n(dummy>d)=dummy(dummy>d);
    % Bu son iki satýr gereksiz olabilir bilmiyorum!!!
end


x=repmat(x_n,m+1,1);
%x is [m+1,8*n]
%from left to right it is extended interpolation points from -4pi to 4pi

powervec=(0:m)';
%powervec is [m+1,1] it is for polynomial degree
power=repmat(powervec,1,8*n);
%power is [m+1,8*n]
%from left to right it is the same power


%This part is for bump function
phi_r=(dummy-b)/(a-b);
phi_p=(dummy-c)/(d-c);
%phi_r,[a,b] & phi_p,[c,d] are from Ecevit&Çaðan page 13
bumpfnc=(1/4)*(psix(phi_r)+1-psix(1-phi_r)).*(psix(phi_p)+1-psix(1-phi_p));
% bumpfnc is [1,6*n]
bump=repmat(bumpfnc,m+1,1);
%This part was for bump function

y=(( (2*x-a-d) / (d-a) ).^power).*bump;
%y is [m+1,8*n]

soln=y(:,1:2*n)+y(:,2*n+1:4*n)+y(:,4*n+1:6*n)+y(:,6*n+1:8*n);

% Normalization==>
Y=( ( sqrt( (pi/n)*(soln*soln') ) .* eye(min(size(soln))) )^(-1) )*soln;


end