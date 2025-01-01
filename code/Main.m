function [  ] = Main(  )

eta_values       = [50,100,200,400,800];
n_values         = [1024,1024,1024,1024,1024*2]; % hslf the number of discretization points
poly_degrees     = [4,8,12,16,20];
no_of_reflection = 40;

configuration_type = 0;    % The ceometrical configuration
no_of_intervals    = 8;
% if change of variable = 0
    %8 means 4 transitions, 1 Illuminated 1 deep shadow
    %7 means 2 transitions only in the illuniminated region, 1 Illuminated 2 deep shadow
    %6 means 2 transitions, 1 Illuminated 1 deep shadow
    %5 means 0 transitions, 1 Illuminated 2 deep shadow
    %4 means 0 transitions, 1 Illuminated 1 deep shadow
    %3 means 0 transitions, 1 Illuminated 0 deep shadow
% if change of variable = 1
    %6 means 2 transitions only in the illuniminated region, 0 Illuminated 2 deep shadow
    %5 means 2 transitions only in the illuniminated region, 0 Illuminated 1 deep shadow

chan_of_var        = 0;    % 0 means no change of variables. 1 otherwise.
kerneltype         = 3;    % 3=combined field, 2= second kind, 1= first kind


all_errors         = cell(2,1); % in all reflection relative errors
gal_ref_sum_errors = cell(1,1); % sum of reflections

% remainder ~~~~~~

for eta_index = 1:max(size(eta_values))
    
    n   = n_values(eta_index);           % number of interpolation points between 0 and pi
    eta = eta_values(eta_index);            % Wavenumber
    t   = 0:pi/n:2*pi-pi/n;
    m=0; % This particular line has no meaning
    
    input = {n,t,eta,m,configuration_type,no_of_intervals,chan_of_var,kerneltype,no_of_reflection};
    %{1=point_number=n, 2=points=t, 3=wavenumber=eta, 4=basis_dimension=m, 5=configuration_type, 6=no_of_intervals, 7=chan_of_var, 8=kerneltype, 9=no_of_reflection}
    
    
    configuration_info = Configurations(input);  %{1=X, 2=DX, 3=DDX, 4=shape_type}
    kernels = Kernels(input,configuration_info); % Usual matrices for kernels. They are without idedtity part.
    
    
    phase = cell(2,1);  % phase function of each reflection.
    sb    = cell(2,1);  % shadow boundary points for each reflection.
    if configuration_type==0
        
        load('Omer-Config0-Path1.mat');
        for i=1:no_of_reflection
            phase{1}{i} = interpft(PHASE{i},n*2);
        end
        % phase{1} now have enough points. for 1st path
        
        sb{1} = order_shadow_boundaries(input,SB,phase{1});
        % shadow boundary pnts for 1st path are ordered according to counterclockwise parametrization
        
        
        load('Omer-Config0-Path2.mat');
        for i=1:no_of_reflection
            phase{2}{i} = interpft(PHASE{i},n*2);
        end
        % phase{2} now have enough points. for 2nd path
        
        sb{2} = order_shadow_boundaries(input,SB,phase{2});
        % shadow boundary pnts for 2nd path are ordered according to counterclockwise parametrization
        
    elseif configuration_type==3
        
        load('Omer-Config3-Path1.mat');
        for i=1:no_of_reflection
            phase{1}{i} = interpft(PHASE{i},n*2);
        end
        sb{1} = order_shadow_boundaries(input,SB,phase{1});
        
        load('Omer-Config3-Path2.mat');
        for i=1:no_of_reflection
            phase{2}{i} = interpft(PHASE{i},n*2);
        end
        sb{2} = order_shadow_boundaries(input,SB,phase{2});
        
    end
    
    
    first_rhs = cell(2,1);
    first_rhs{1} = FirstRHS(input,phase{1},configuration_info,1); 
    % the right hand side in the main integral equation for first  obj.
    first_rhs{2} = FirstRHS(input,phase{2},configuration_info,2); 
    % the right hand side in the main integral equation for second obj.
    
    
    % Nysterm solution of the whole problem. denoted as "exact_soln"
    big_kernel    = eye(4*n) - [ kernels{1,1},kernels{1,2} ; kernels{2,1},kernels{2,2} ];
    big_first_rhs = [first_rhs{1} ; first_rhs{2}];
    exact_soln = linsolve(big_kernel,big_first_rhs);
    exact_soln = exact_soln.';
    
    
    % Nysterm solutions part
    exact_reflctn_rhs  = cell (1,2);
    exact_reflctn      = cell (1,2);
    % Nysterm solutions of the reflection problems
    
    for i=1:no_of_reflection
        j=mod(i,2)+1;   % for the path starting from object 2. Gives the object number we are working on i-th step
        k=mod(i+1,2)+1; % for the path starting from object 1. Gives the object number we are working on i-th step
        if i==1
            exact_reflctn_rhs{1}{i} = first_rhs{1}; % for first  path
            exact_reflctn_rhs{2}{i} = first_rhs{2}; % for second path
        else
            exact_reflctn_rhs{1}{i} = kernels{k,j} * exact_reflctn{1}{i-1}.'; % for first  path ith reflection
            exact_reflctn_rhs{2}{i} = kernels{j,k} * exact_reflctn{2}{i-1}.'; % for second path ith reflection
        end
        
        exact_reflctn{1}{i} = Exact_solution_for_one_reflection( input,kernels,exact_reflctn_rhs{1}{i},i,1 ); % 1 here represent it is from first  path.
        exact_reflctn{2}{i} = Exact_solution_for_one_reflection( input,kernels,exact_reflctn_rhs{2}{i},i,2 ); % 2 here represent it is from second path.
    end
    
    % Sum of the Nyström solutions
    exact_reflection_sum{1}   = zeros(1,2*n); % total sum on 1st object.
    exact_reflection_sum{2}   = zeros(1,2*n); % total sum on 2nd object.
    % exact_reflection_sum_hist = cell(1,1); % {i}{2} sum of i many reflections for the second object.
    
%     %Summation of the Nysterm solutions. We dont really need this part.
%     for i=1:no_of_reflection
%         k=mod(i+1,2)+1 ;
%         j=mod( i ,2)+1 ;
%         exact_reflection_sum{1} = exact_reflection_sum{1} + exact_reflctn{k}{i}; % 1st path. Nyström solns
%         exact_reflection_sum{2} = exact_reflection_sum{2} + exact_reflctn{j}{i}; % 2nd path. Nyström solns
%         
%         % exact_reflection_sum_hist{i}    = {exact_reflection_sum{1},exact_reflection_sum{2}};    
%         % % This last one records each sum(sum of 1 reflection,...,sum of all reflections).
%         % % It is not used. It was for research.
%     end
    
    
%     [gal_rem_1,gal_rem_2] = Remainders(input,kernels,exact_reflctn{1}{1},exact_reflctn{1}{2},no_of_reflection);
% 
%     compare = exact_soln - [exact_reflection_sum{1},exact_reflection_sum{2}];
%     log10( sqrt( (compare*compare') / (exact_soln*exact_soln') ) )
%         
%         compare = exact_soln - [exact_reflection_sum{1},exact_reflection_sum{2}] - [gal_rem_1,gal_rem_2];
%         log10( sqrt( (compare*compare') / (exact_soln*exact_soln') ) )
    
%     compare=exact_soln(1:2*n) - exact_reflection_sum{1};
%     log10( sqrt( (compare*compare') / (exact_soln(1:2*n)*exact_soln(1:2*n)') ) )
%     compare=exact_soln(2*n+1:4*n) - exact_reflection_sum{2};
%     log10( sqrt( (compare*compare') / (exact_soln(1+2*n:4*n)*exact_soln(1+2*n:4*n)') ) )
    
%     figure
%     subplot(2,1,1)
%     plot(input{2},exact_soln(1:2*n))
%     subplot(2,1,2)
%     plot(input{2},exact_soln(2*n+1:4*n))
%     figure
%     subplot(2,2,3)
%     plot(input{2},exact_reflection_sum{1})
%     subplot(2,2,4)
%     plot(input{2},exact_reflection_sum{2})
    
    
    % Galerkin solutions
    % runs over all poly degrees for a specific wave number
    for poly_index = 1:max(size(poly_degrees))
        
        t   = 0:pi/n:2*pi-pi/n;             % "load file" disturbs this variable
        m   = poly_degrees(poly_index);     % 0,x,...x^m are our basis functions
                                            % input now includes poly. degrees It was zero before for technical reasons 
                                            
        input = {n,t,eta,m,configuration_type,no_of_intervals,chan_of_var,kerneltype,no_of_reflection};
        %{1=point_number=n, 2=points=t, 3=wavenumber=eta, 4=basis_dimension=m, 5=configuration_type, 6=no_of_intervals, 7=chan_of_var, 8=kerneltype, 9=no_of_reflection}
        
        
        galerkin_soln = cell(2,1);  % Galerkin soln for each reflection
        ksi_changes   = cell(2,1);  % We use this after regions are almost fixed so that we do not optimize ksi values again again. They are same after some reflection
        error         = cell(2,1);  % Just temporary variable. Holds relative errors of each reflection.
        galerkin_rhs  = cell(2,1);  % We use this to get i-th Galerkin solution from the i-1-th galerkin solution
        % Here none of these are not recorded. As wave number changes old informations are lost.
        % We store only the errors under the name "all_errors"
        
        % Galerkin solution part
        for i=1:no_of_reflection
            j=mod(i,2)+1;   % for reflection starting from object 2. Gives the object no we are working on i-th step
            k=mod(i+1,2)+1; % for reflection starting from object 1. Gives the object no we are working on i-th step
            
            % This part computes RHS
            if i==1         % In first reflection we use right hand side of the main integral eqn
                galerkin_rhs{1}{1} = first_rhs{1};  % path 1
                galerkin_rhs{2}{1} = first_rhs{2};  % path 2
            else            % In the other reflection we use right hand siden coming from Neumann expansion
                % HERE WE USE GALERKIN SOLUTIONS AS A RIGHT HAND SIDE
                % IF WE WANT TO SOLVE THE REFLECTION BY NYSTRÖM RIGHT
                % HAND SIDE "galerkin_soln{1}{i-1}" SHOULD BE CHANGED BY 
                % "exact_reflctn{1}{i-1}". I KNOW IT LOOKS UGLY :(
                galerkin_rhs{1}{i} = kernels{k,j} * galerkin_soln{1}{i-1}.'; % 1st path
                galerkin_rhs{2}{i} = kernels{j,k} * galerkin_soln{2}{i-1}.'; % 2nd path
            end
            
            % In this part we find ith {galerkin solution, relative error
            % and ksi changes}. Ksi changes are used after reflection is
            % stabilized.
            if i < 11 % This depends configuration and must be obtained for each conf. by looking at shadow boundaries.
                
                % 0's at the very end means we dont fix ksi values because shadow boundaries are not stable yet
                [galerkin_soln{1}{i}, error{1}{i}, ksi_changes{1}{i}] ... 
                    = Reflection_solver (input,configuration_info,kernels,...
                    galerkin_rhs{1}{i},phase{1}{i},sb{1}{i},k,exact_reflctn{1}{i},i,0); % 1st path
                
                [galerkin_soln{2}{i}, error{2}{i}, ksi_changes{2}{i}] ...
                    = Reflection_solver (input,configuration_info,kernels,...
                    galerkin_rhs{2}{i},phase{2}{i},sb{2}{i},j,exact_reflctn{2}{i},i,0); % 2nd path
            else
                
                % ksi_changes's at the very end means shadow boundaries are stabilized so we can use previous intervals.
                [galerkin_soln{1}{i}, error{1}{i}, ksi_changes{1}{i}] ... 
                    = Reflection_solver (input,configuration_info,kernels,...
                    galerkin_rhs{1}{i},phase{1}{i},sb{1}{i},k,exact_reflctn{1}{i},i,ksi_changes{1}{i-2}); % 1st path
                
                [galerkin_soln{2}{i}, error{2}{i}, ksi_changes{2}{i}] ...
                    = Reflection_solver (input,configuration_info,kernels,...
                    galerkin_rhs{2}{i},phase{2}{i},sb{2}{i},j,exact_reflctn{2}{i},i,ksi_changes{2}{i-2}); % 2nd path
            end
            
            % In this part we store relative errors for each reflection
            all_errors{1}{eta_index}(poly_index,i)=error{1}{i}; % 1st path
            all_errors{2}{eta_index}(poly_index,i)=error{2}{i}; % 2nd path
        end
        
        
        galerkin_sum{1} = zeros(1,2*n); % total sum on 1st object.
        galerkin_sum{2} = zeros(1,2*n); % total sum on 2nd object.
        gal_ref_sum_errors{eta_index,poly_index} = zeros(1,no_of_reflection); % {i}{2} = ith relative on the second object.
        % galerkin_sum_hist = cell(1,1);  % {i}{2} = ith summation on the
        % second object. We dont use this anymore
        
        
        % Galerkin solutions are summed here
        for i=1:no_of_reflection
            k=mod(i+1,2)+1 ; % for reflection starting from object 1. Gives the object no we are working on i-th step
            j=mod( i ,2)+1 ; % for reflection starting from object 2. Gives the object no we are working on i-th step
            
            galerkin_sum{1} = galerkin_sum{1} + galerkin_soln{k}{i}; % For 1st object
            galerkin_sum{2} = galerkin_sum{2} + galerkin_soln{j}{i}; % For 2nd object
            
            % galerkin_sum_hist{i}{1} = galerkin_sum{1}; % Stores whole sumation history. path 1. We dont use it any more.
            % galerkin_sum_hist{i}{2} = galerkin_sum{2}; % Stores whole sumation history. path 2. We dont use it any more.
            
            % Relative error of the big picture
            compare = exact_soln - [galerkin_sum{1},galerkin_sum{2}] ;
            gal_ref_sum_errors{eta_index,poly_index}(i) = log10( sqrt( (compare*compare') / (exact_soln*exact_soln') ) );
            % Stores whole relative error history. Here we compare the
            % galerkin sums with Nyström solution which is found WITHOUT
            % Neumann expansion. As "i" gets bigger we get better solutions
            % until some reflection number. After that this relative error
            % stabilizes.
        end
        
%         [gal_rem_1,gal_rem_2] = Remainders(input,kernels,galerkin_soln{1}{1},galerkin_soln{2}{1},no_of_reflection);
%         
%         compare = exact_soln - [galerkin_sum{1},galerkin_sum{2}] - [gal_rem_1,gal_rem_2];
%         log10( sqrt( (compare*compare') / (exact_soln*exact_soln') ) )
%         
%         figure
%         subplot(2,1,1)
%         plot(t,gal_rem_1)
%         subplot(2,1,2)
%         plot(t,gal_rem_2)
        
%         compare=exact_soln(1:2*n) -  galerkin_sum{1};
%         log10( sqrt( (compare*compare') / (exact_soln*exact_soln') ) )
%         compare=exact_soln(2*n+1:4*n) -  galerkin_sum{2};
%         log10( sqrt( (compare*compare') / (exact_soln*exact_soln') ) )
%         
%         figure
%         subplot(2,1,1)
%         plot(input{2},galerkin_sum{1})
%         subplot(2,1,2)
%         plot(input{2},galerkin_sum{2})
    end
end

% saving is good and safe :)
save('all_errors.mat','all_errors') % relative errors using Nyström reflection solution to compare.
save('gal_ref_sum_errors.mat','gal_ref_sum_errors') % Reflection sum errors using "Nyström solution without using Neumann expansion"
% save('galerkin_sum_hist.mat','galerkin_sum_hist')   % I really dont use this. ama hatrý var iyi iþe yaramýþtý


% plotting styles...
mrkr = {'g','r','m','b','k'};
titles = {'Reflection : 0','Reflection : 1','Reflection : 10','Reflection : 11','Reflection : 20','Reflection : 21'};
titles2 = {'Local Polynomial degree: 4','Local Polynomial degree: 8','Local Polynomial degree: 12','Local Polynomial degree: 16','Local Polynomial degree: 20'};
line = {'-+g','-sr','-dm','-ob','->k'};
line2 = {'+g','sr','dm','ob','>k'};
line3 = {'-g','-r','-m','-b','-k'};


figure % summation figure
subplot(2,3,1)
plot(configuration_info{1}{1}(1,:),configuration_info{1}{1}(2,:),...
    configuration_info{1}{2}(1,:),configuration_info{1}{2}(2,:),...
    'LineWidth',1.5)
text(-2,-0.5,'\rightarrow','FontSize',30)
text(-2,-1.5,'\rightarrow','FontSize',30)
text(-2,-2.5,'\rightarrow','FontSize',30)

title('Configuration', 'FontSize', 16,'FontWeight','normal')
axis equal
set(gca,'FontSize',13)
box on


for poly_index = 1:max(size(poly_degrees))
    
    subplot(2,3,poly_index+1)
    hold on
    for eta_index = 1:max(size(eta_values))
        plot(0:0,gal_ref_sum_errors{eta_index,poly_index}(1:1),...
            line{eta_index},'MarkerFaceColor',mrkr{eta_index}, 'MarkerSize',6, 'LineWidth',1.3)
    end
    title(titles2{ poly_index }, 'FontSize', 16,'FontWeight','normal')
    if poly_index==1
        legend('k=50','k=100','k=200','k=400','k=800','Location','SouthEast')
    else
        legend('k=50','k=100','k=200','k=400','k=800','Location','NorthEast')
    end
    xlabel('Number of reflections','Fontsize',16)
    set(gca,'FontSize',13)
    
    for eta_index = 1:max(size(eta_values))
        plot(0:(no_of_reflection-1),gal_ref_sum_errors{eta_index,poly_index}(1:no_of_reflection),...
            line3{eta_index},'MarkerFaceColor',mrkr{eta_index}, 'MarkerSize',6, 'LineWidth',1.3)
    end
    ac=gca;
    ac.XTick = 0:10:100;
    axis([0 50  -7 -0.3])
    set(gca,'FontSize',13)
    for eta_index = 1:max(size(eta_values))
        plot(0:5:(no_of_reflection-1),gal_ref_sum_errors{eta_index,poly_index}(1:5:no_of_reflection),...
            line2{eta_index}, 'MarkerFaceColor',mrkr{eta_index}, 'MarkerSize',6)
    end    
    set(gca,'FontSize',13)
    box on
    hold off
    
end

% reflrction figures of paths
for obstacle = 1:2
    figure
    for i=1:3 % for 0th, 10th, and 20th reflection
        subplot(2,3,i)
        hold on
        for eta_index = 1:max(size(eta_values))
            plot( poly_degrees, all_errors{obstacle}{eta_index}( :, (i-1)*10 +1 ),line{eta_index}, 'MarkerFaceColor',mrkr{eta_index}, 'MarkerSize',6, 'LineWidth',1.3  )
        end
        title(titles{ 2*i-1 }, 'FontSize', 16,'FontWeight','normal')
        legend('k=50','k=100','k=200','k=400','k=800','Location','SouthWest')
        axis([3 21 -8 -2])
        xlabel('Local polynomial degree','Fontsize',16)
        ac=gca;
        ac.XTick = 4:4:20;
        ac.YTick = -8:-2;
        set(gca,'FontSize',13)
        box on
        hold off
        
        subplot(2,3,i+3)
        hold on
        for eta_index = 1:max(size(eta_values))
            plot( poly_degrees, all_errors{obstacle}{eta_index}( :, (i-1)*10 +2 ),line{eta_index}, 'MarkerFaceColor',mrkr{eta_index}, 'MarkerSize',6, 'LineWidth',1.3 )
        end
        title( titles{ 2*i }, 'FontSize', 16,'FontWeight','normal')
        legend('k=50','k=100','k=200','k=400','k=800','Location','SouthWest')
        axis([3 21 -8 -2])
        xlabel('Local polynomial degree','Fontsize',16)
        ac=gca;
        ac.XTick = 4:4:20;
        ac.YTick = -8:-2;
        set(gca,'FontSize',13)
        box on
        hold off
    end
end


end

