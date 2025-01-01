function [ configuration_info ] = Configurations( input )
%{1=point_number=n, 2=points=t, 3=wavenumber=eta, 4=basis_dimension=m, 5=configuration_type, 6=no_of_intervals, 7=chan_of_var, 8=kerneltype, 9=no_of_reflection}

% configuration inf. is entered here by hand



configuration_type=input{5};
t=input{2};

if configuration_type==0
    X11 = cos(t);
    X12 = sin(t);
    X21 = cos(t);
    X22 = sin(t)-3;
    
    DX11 = -sin(t);
    DX12 = cos(t);
    DX21 = -sin(t);
    DX22 = cos(t);
    
    DDX11= -cos(t);
    DDX12= -sin(t);
    DDX21= -cos(t);
    DDX22= -sin(t);
    
    X{1} = [X11;X12];
    X{2} = [X21;X22];
    DX{1} = [DX11;DX12];
    DX{2} = [DX21;DX22];
    DDX{1} = [DDX11;DDX12];
    DDX{2} = [DDX21;DDX22];
    
    shape_type = [1,1];
    % 1 = cicle, 2 = ellipse
    % This part is not necessary now but we can not delate it yet.Code uses
    % this in a lot places and would certainly give error.
    
    shape_centers = { [0,0] , [0,-3] };
    % this is used for drawing shapes. It has no use in problem solving
    
elseif configuration_type==3
    
    EX11 = 3/2*cos(t);
    EX12 = 1/2*sin(t);
    EX21 = 3/2*cos(t);
    EX22 = 1/3*sin(t);
    
    EDX11 = -3/2*sin(t);
    EDX12 = 1/2*cos(t);
    EDX21 = -3/2*sin(t);
    EDX22 = 1/3*cos(t);
    
    theta_1 = pi/6;
    theta_2 = pi/12;
    
    X11 = cos(theta_1)*EX11-sin(theta_1)*EX12;
    X12 = sin(theta_1)*EX11+cos(theta_1)*EX12;
    X21 = 0.5+cos(theta_2)*EX21-sin(theta_2)*EX22;
    X22 = -1.6+sin(theta_2)*EX21+cos(theta_2)*EX22;
    
    DX11  = cos(theta_1)*EDX11-sin(theta_1)*EDX12;
    DX12  = sin(theta_1)*EDX11+cos(theta_1)*EDX12;
    DX21  = cos(theta_2)*EDX21-sin(theta_2)*EDX22;
    DX22  = sin(theta_2)*EDX21+cos(theta_2)*EDX22;
    
    DDX11= cos(theta_1)*EX11+sin(theta_1)*EX12;
    DDX12= sin(theta_1)*EX11-cos(theta_1)*EX12;
    DDX21= cos(theta_2)*EX21+sin(theta_2)*EX22;
    DDX22= sin(theta_2)*EX21-cos(theta_2)*EX22;    
    
    X{1} = [X11;X12];
    X{2} = [X21;X22];
    DX{1} = [DX11;DX12];
    DX{2} = [DX21;DX22];
    DDX{1} = [DDX11;DDX12];
    DDX{2} = [DDX21;DDX22];
    
    
    shape_type = [1,1];
    
    shape_centers = { [0,0] , [0.5,-1.6] };
    
end

configuration_info={X,DX,DDX,shape_type,shape_centers};



end

