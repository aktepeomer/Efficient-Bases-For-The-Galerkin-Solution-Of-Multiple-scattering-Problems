function [ kernels ] = Kernels( input , configuration_info)
%input = {1=point_number=n, 2=points=t, 3=wavenumber=eta, 4=basis_dimension=m, 5=configuration_type, 6=no_of_intervals, 7=chan_of_var, 8=kerneltype, 9=no_of_reflection}
%configuration_info = {1=X, 2=DX, 3=DDX, 4=shape_type, 5=shape_centers}

% kernels are without "I" part

kernels=cell(2,2);% kernel{i}{j} means the usual definition (without identity.)


n = input{1};
t = input{2};
eta = input{3};
kerneltype = input{8};

C=0.5772156649015328606065120900824024310421;

% x1 = configuration_info{1}{1}(1,:);
% x2 = configuration_info{1}{1}(2,:);
% dx1 = configuration_info{2}{1}(1,:);
% dx2 = configuration_info{2}{1}(2,:);
% ddx1 = configuration_info{3}{1}(1,:);
% ddx2 = configuration_info{3}{1}(2,:);
%
% y1 = configuration_info{1}{2}(1,:);
% y2 = configuration_info{1}{2}(2,:);
% dy1 = configuration_info{2}{2}(1,:);
% dy2 = configuration_info{2}{2}(2,:);
% ddy1 = configuration_info{3}{2}(1,:);
% ddy2 = configuration_info{3}{2}(2,:);

for i=1:2 % At this loop we define kernels 11 and 22
    
    % i-th object information
    x1 = configuration_info{1}{i}(1,:);
    x2 = configuration_info{1}{i}(2,:);
    dx1 = configuration_info{2}{i}(1,:);
    dx2 = configuration_info{2}{i}(2,:);
    ddx1 = configuration_info{3}{i}(1,:);
    ddx2 = configuration_info{3}{i}(2,:);
    
    
    % Dimensions are [2*n,2*n]. "t" is from Colton-Kress and Numerical
    % Implementation. Definitions of these terms are obvious.
    t_t=repmat(t.',1,2*n);
    x1t=repmat(x1.',1,2*n);
    x2t=repmat(x2.',1,2*n);
    dx1t=repmat(dx1.',1,2*n);
    dx2t=repmat(dx2.',1,2*n);
    ddx1t=repmat(ddx1.',1,2*n);
    ddx2t=repmat(ddx2.',1,2*n);
    norm_dxt=repmat(twonorm([dx1;dx2]).',1,2*n);
    %                                                               .
    % this specific one is just for following the steps more easly.|x(t)|.
    % Top to buttom "t" increses.
    % Columns are same.
    
    % Dimensions are [2*n,2*n]. "tao" is from Colton-Kress and Numerical
    % Implementation. Definitions of these terms are obvious.
    t_tao=repmat(t,2*n,1);
    x1tao=repmat(x1,2*n,1);
    x2tao=repmat(x2,2*n,1);
    norm_dxtao=repmat(twonorm([dx1;dx2]),2*n,1);
    %                                                               .
    % this specific one is just for following the steps more easly.|x(tao)|.
    % Left to right "t" increses.
    % Rows are same.
    
    % Dimensions are [2*n,2*n].
    dummy1=x1t-x1tao;
    dummy2=x2t-x2tao;
    r=sqrt(dummy1.^2+dummy2.^2);
    % r is |x(t)-x(tao)| and dimension is [2*n,2*n].
    % Name r comes from Colton-Kress page 68.
    % From left to right "tao" changes and from top to buttom "t" changes.
    
    
    % This part is for finding R^n_|i-j| in Colton-Kress page 70 which is used
    % in approximating a "log" term in Colton-Kress page 69.
    m=1:1:n-1;
    % m is in the definition of R^n_|i-j|.
    R_first_row=(-2*pi/n)*(1./m)*cos(m'*t)-(pi/n^2)*((-1).^(0:2*n-1));
    % This is the first row of R. The rest of the rows are just combinations of
    % the entries of the first row. So we can find them by this one.
    comb=abs(repmat((0:2*n-1)',1,2*n)-repmat((0:2*n-1),2*n,1))+ones(2*n,2*n);
    % This is |i-j|+1. +1 is needed to find the corresponding entry in 1st row.
    % R_first_row( (comb)ij )=(R)ij where "R" is our R^n_|i-j| in Colton-Kress.
    R_n=R_first_row(comb.');
    % This is the R^n_|i-j| matrix[2*n,2*n] we want for i,j=0,...,2*n-1.
    % The reason we write Comb' is a matlab technicality (but we do not need to
    % do this because Comb is a symmetric matrix).
    
    if kerneltype==1||kerneltype==3
        % For finding M(t,tao): [2*n,2*n], from Numerical Implementation.
        M=(0.5*complex(0,1))*(besselh(0,eta*r));
        % We have a different definition for (t,t) so we shall change the diogonal.
        M(logical(eye(size(M)))) = 0;
        % This part is for setting diogonal part 0.
        M1=(-0.5/pi)*(besselj(0,eta*r));
        M2=(M - M1.*log(4*sin(0.5*(t_t-t_tao)).^2));
        M2(logical(eye(size(M2)))) = 0;
        M2d=0.5*complex(0,1)-(1/pi)*(C+log(0.5*eta*norm_dxt));
        % This is the diogonal part of M2. All columns are equal to each other. The
        % reason is just a technicality. It is .* by identity matrix below so it
        % turns into the form we want.
        M2=M2 + M2d.*eye(size(M2));
    end
    if kerneltype==2||kerneltype==3
        % For finding L(t,tao): [2*n,2*n], from Numerical Implementation.
        L=(0.5*complex(0,1)*eta)*besselh(1,eta*r).*((dx2t.*(x1t-x1tao)-dx1t.*(x2t-x2tao))./r)./norm_dxt;
        % "nu" in the definition is outward normal vector of the shape which is
        % [dx_2,-dx_1]/|[dx_1,dx_2]|
        % We have a different definition for (t,t) so we shall change the diogonal.
        L(logical(eye(size(L)))) = 0;
        % This part is for setting diogonal part 0.
        % There is no need to define diogonal terms here since we will define
        L1=(-0.5*eta/pi)*(besselj(1,eta*r)).*((dx2t.*(x1t-x1tao)-dx1t.*(x2t-x2tao))./r)./norm_dxt;
        L1(logical(eye(size(L1)))) = 0;
        L2=L - L1.*log(4*sin(0.5*(t_t-t_tao)).^2);
        L2(logical(eye(size(L2)))) = 0;
        L2d=(-0.5/pi)*(dx2t.*ddx1t-dx1t.*ddx2t)./(norm_dxt.^3);
        % This is the diogonal part of L2. All rows are equal to each other. The
        % reason is just a technicality. It is .* by identity matrix below so it
        % turns into the form we want.
        L2=L2 + L2d.*eye(size(L2));
    end
    
    switch kerneltype
        case 1 % 1st kind
            K1=M1;
            K2=M2;
        case 2 % 2nd kind
            K1=L1;
            K2=L2;
        case 3 % combined field
            K1=L1 + complex(0,1)*eta*M1;
            K2=L2 + complex(0,1)*eta*M2;
    end
    
    kernels{i,i}=(R_n.*K1 + (pi/n)*K2).*norm_dxtao;
    % The last part comes from change of variable (Numerical Implementation).
    % It is[2*n,2*n].
end

% Burdan aþaðýsý hatalý

for i=1:2 % At this loop we define kernel 12 and 21
    % y üzerinden integral alýp x üzerine tanýmlý
    x1 = configuration_info{1}{i}(1,:);
    x2 = configuration_info{1}{i}(2,:);
    dx1 = configuration_info{2}{i}(1,:);
    dx2 = configuration_info{2}{i}(2,:);
    ddx1 = configuration_info{3}{i}(1,:);
    ddx2 = configuration_info{3}{i}(2,:);
    
    j=mod(i,2)+1; % i=1 then j=2   or   i=2 then j=1
    
    y1 = configuration_info{1}{j}(1,:);
    y2 = configuration_info{1}{j}(2,:);
    dy1 = configuration_info{2}{j}(1,:);
    dy2 = configuration_info{2}{j}(2,:);
    ddy1 = configuration_info{3}{j}(1,:);
    ddy2 = configuration_info{3}{j}(2,:);
    
    
    % Definition of the following variables are different from prev. fncs.
    X1=repmat(x1.',1,2*n);
    X2=repmat(x2.',1,2*n);
    Y1=repmat(y1,2*n,1);
    Y2=repmat(y2,2*n,1);
    DX1=repmat(dx1.',1,2*n);
    DX2=repmat(dx2.',1,2*n);
    
    normDX=repmat(twonorm([dx1;dx2]).',1,2*n);
    normDY=repmat(twonorm([dy1;dy2])  ,2*n,1);
    
    dummy1=X1-Y1;
    dummy2=X2-Y2;
    r=sqrt(dummy1.^2+dummy2.^2);
    
    if kerneltype==2 || kerneltype==3
        L=complex(0,1)*eta*(1/2)*besselh(1,eta*r).* ( (DX2.*(X1-Y1)-DX1.*(X2-Y2))./r) ./normDX;
    end
    if kerneltype==1 || kerneltype==3
        M=complex(0,1)*(1/2)*besselh(0,eta*r);
    end
    
    switch kerneltype
        case 1 % 1st kind
            K=M;
        case 2 % 2nd kind
            K=L;
        case 3 % combined filed
            K=L + complex(0,1)*eta*M;
    end
    
    kernels{i,j}=K*(pi/n).*normDY;
end


end

