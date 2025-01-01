function [  ] = Picture( input, configuration_info, sb, shape_no, ksi_changes, no_of_ref, obs )
%{1=point_number=n, 2=points=t, 3=wavenumber=eta, 4=basis_dimension=m, 5=configuration_type, 6=no_of_intervals, 7=chan_of_var, 8=kerneltype, 9=no_of_reflection}
%configuration_info = {1=X, 2=DX, 3=DDX, 4=shape_type, 5=shape_centers}

% This fnc gives pictures about intervals.
% It was disabled in the code.

t1= sb(1,1)  ;
t2= sb(2,1)  ;

n=input{1};
k=input{3};
noofintervals=input{6};
chan_of_var=input{7};


shapetype=configuration_info{4}(shape_no);
x=configuration_info{1}{shape_no};
shape_centers=configuration_info{5};

center_x=shape_centers{shape_no}(1);
center_y=shape_centers{shape_no}(2);

[y,~]=intervals_with_optimization( t1,t2,noofintervals,k,shapetype,chan_of_var,ksi_changes );

% a=zeros(noofintervals,2);
% for i=1:noofintervals
%     a(i,1)=y(i,1);
%     a(i,2)=y(i,4);
% end
% a/(pi);
% mod(a,2)

dummy=-4*pi:pi/n:4*pi-pi/n;
bump=zeros(noofintervals,2*n);

%This part is for bump function
for i=1:noofintervals
    a=y(i,1);
    b=y(i,2);
    c=y(i,3);
    d=y(i,4);
    phi_r=(dummy-b)/(a-b);
    phi_p=(dummy-c)/(d-c);
    %phi_r,[a,b] & phi_p,[c,d] are from Ecevit&Çaðan page 13
    bumpfnc=(1/4)*(psix(phi_r)+1-psix(1-phi_r)).*(psix(phi_p)+1-psix(1-phi_p));
    bump(i,:)=bumpfnc(1:2*n)+bumpfnc(2*n+1:4*n)+bumpfnc(4*n+1:6*n)+bumpfnc(6*n+1:8*n);
end


% This figure name is in order to determine which figure belongs to which
% case. read orde : (shapeno reflection number eta)
figure (k + 1000*input{4} + 1000000*no_of_ref + 100000000*shape_no)
hold on
for i=1:noofintervals
    plot( (bump(i,:)*0.1+1).*(x(1,:)-center_x)+center_x , center_y+(x(2,:)-center_y).*(bump(i,:)*0.1+1) )
end
plot(x(1,:),x(2,:),'k')
text(sb(1,2),sb(1,3),'t1')
text(sb(2,2),sb(2,3),'t2')
axis('equal');
hold off

end

