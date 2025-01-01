function [ galerkin_soln,errormin,ksi_changes ] = Reflection_solver( input, configuration_info, kernels, rhs, phase, sb, shape_no, exact_reflctn, no_of_ref, stabilized_ksi )
%{1=point_number=n, 2=points=t, 3=wavenumber=eta, 4=basis_dimension=m, 5=configuration_type, 6=no_of_intervals, 7=chan_of_var, 8=kerneltype, 9=no_of_reflection}
%{1=X, 2=DX, 3=DDX, 4=shape_type}

% Here phase is just for obstacle = shape_no and at reflection = reflection_no
% Here sb    is just for obstacle = shape_no and at reflection = reflection_no

n                  = input{1};
eta                = input{3};
m                  = input{4};
no_of_intervals    = input{6};
chan_of_var        = input{7};
kerneltype         = input{8};

shape_type = configuration_info{4}(shape_no);
% circle=1 or ellipse=2
% it has no use now. all set to 1 

t1 = sb(1,1); % sb points as angles
t2 = sb(2,1); % sb points as angles

% choose krnel to work
kernel = kernels{shape_no,shape_no};
exact = exact_reflctn; % just to use short notation

% standard e^... part
e_term = exp(complex(0,1)*eta*phase);

% ýncreses degrees for deep shadow. may be unnecessary
max_poly_degree=repmat(m,no_of_intervals,1);
if (chan_of_var==1 && no_of_intervals==5) || (chan_of_var==0 && (no_of_intervals==4 || no_of_intervals==6 )); % Maybe we can add no_of_intervals==8 too?
    max_poly_degree(1)=max_poly_degree(1)*2+1;
end


basisfncs=zeros(sum(max_poly_degree+1),2*n);


% width=[0.1 0.05]; % The limit of changes we are going to try.
% inc  =[ 0.05 0.025];
% width=[ 0.1];
% inc  =[ 0.1 0.025];% changes increases by this amount
if no_of_ref==1 && eta==50 && chan_of_var==1 % for some reason first reflection is computed good but is not ploted at all
    width=0;% 0.05]; % The limit of changes we are going to try.
    inc  =0.025;% 0.025];
else
    width=[0.05 0.05]; % The limit of changes we are going to try.
    inc  =[0.05 0.025];
end
% this part may chance by rearranging optimization process
% This is what we are doing, consider it as intervals:
% w..............originalksi..............w
% c....c....c....originalksi....c....c....c
% we are trying each "c" value as new ksi below and choose the best one

% this one finds dimentions of the ksi matrix
[~,ksi_orj,~]=intervals_with_optimization(t1,t2,no_of_intervals,eta,shape_type,chan_of_var,0);
ksi_changes=zeros( size(ksi_orj) );
% If ksi_changes is one row it means circle. If 2, then ellipse.
% The difference comes from loosing the symmetry.
% ksi_orj is original ksi values we take from previous works.
% "ksi_changes" RESTROES ALL OPTIMAL KSI CHANGES. INITIALLY IT IS ZERO
% OFCOURSE


if no_of_ref < 11 % untill this many reflection shadow boundaries are not stabilized
    for w=1:length(width)
        
        % For circle and ellipse it varies (technical issue discussed above)
        for k=1:min(size(ksi_changes))
            for ksi_index=1:length(ksi_changes(k,:))
                
                % If ksi_orj(i)=0, this means that ksi is not used in the
                % computation. So no need to consider that case.
                % technical issue it will be fixed when intervalswithoptm.
                % fnc is changed
                if ksi_orj(k,ksi_index)~=0
                    
                    change=-width(w):inc(w):width(w);
                    % All the ksi changes we try this step. we just
                    % increase or decrease ksi values this much
                    
                    errors=zeros( size(change) );
                    % All realtive errors we get by corresponding ksi change.
                    
                    % Now we try some changes for one specific ksi.
                    for c=1:length(change)
                        
                        ksi_changes(k,ksi_index)=ksi_changes(k,ksi_index)+change(c);
                        % This is not ksi, it is the change of the ksi.
                        % below we will undo this change so that we do not
                        % have a technical problems.
                        
                        [intervalMatrix,ksi_on_use,validity]=intervals_with_optimization(t1,t2,no_of_intervals,eta,shape_type,chan_of_var,ksi_changes);
                        % ksi_on_use is all our current ksi values.
                        % matrix made by intervals. jth column means:
                        % (j,1)=a, (j,2)=a' (j,3)=b', (j,4)=b for jth interval.
                        % [noofintervals,4]
                        
                        % if the ksi changes are appropriate we continue
                        if validity == 1
                            
                            j=1;
                            for i=1:no_of_intervals
                                basisfncs(j:j+max_poly_degree(i),:)...
                                    =Basisfncs(intervalMatrix(i,:),n,max_poly_degree(i),chan_of_var,ksi_on_use,eta,no_of_intervals,i,t1,t2,shape_type);
                                j=j+max_poly_degree(i)+1;
                            end
                            % On each loop basis functions on one interval is find.
                            % each row is a func.
                            % "j" is just a technicality.
                            % From Numerical Implementation page 14. [1,2*n]
                            longbasis=basisfncs.*repmat(e_term,sum(max_poly_degree+1),1);
                            % basis fncs * usual "e^..." term for each row
                            
                            % we will solve the problem A*x=b coming from
                            % galerkin matrix
                            switch kerneltype
                                case 1 % for 1st kind galerkin matrix
                                    A=conj(longbasis)*kernel*longbasis.';
                                case 2 % for 2nd kind galerkin matrix
                                    A=(conj(longbasis)*longbasis.')-(conj(longbasis)*kernel*longbasis.');
                                case 3 % for combined field galerkin matrix
                                    A=(conj(longbasis)*longbasis.')-(conj(longbasis)*kernel*longbasis.');
                            end
                            % From Numerical Implementation page 14, the whole matrix in the end. Dimention is [basisnumber,basisnumber]
                            
                            b=conj(longbasis)*rhs;
                            % Dimention is [basisnumber,1]
                            
                            % Preconditioning
                            S=eye(size(A)).*A;
                            A=S^(-1)*A;
                            b=S^(-1)*b;
                            
                            alfa=linsolve(A,b);
                            % Coefficients of basis functions. Dimention is [basisnumber,1]
                            
                            soln=sum(longbasis.*repmat(alfa,1,2*n),1);
                            % This is our approximation to the exact soln. [1,2*n]
                            
                            
                            compare=(soln-exact);
                            y=log10( sqrt((compare*compare') / (exact*exact')) );
                            % relative error
                            
                            % store the relative error
                            errors(c)=y;
                            
                            % undo the ksi changes done in this loop for this spesific ksi so
                            % that prev ksi changes are not effected and
                            % the next trial for this specisic ksi change
                            % is like never done. so that "ksi_changes(k,ksi_index)=ksi_changes(k,ksi_index)+change(c);"
                            % part is not effected.
                            ksi_changes(k,ksi_index)=ksi_changes(k,ksi_index)-change(c);
                            
                        else % if ksi changes are not appropriate we reject it like this
                            
                            % ksi change we applied was bad so we assign
                            % error really big
                            errors(c)=200;
                            
                            % thia part is explained above
                            ksi_changes(k,ksi_index)=ksi_changes(k,ksi_index)-change(c);
                        end
                    end
                    
                    % This part gets rid of +-Nan and +-inf cases.
                    errors(isinf(errors))=100;
                    errors(isnan(errors))=100;
                    
                    errormin=min(errors);
                    % We find the minimum relative eror after trying different values for one specific ksi.
                    
                    % this part was really problematic i dont know how to
                    % write it in an elegent way
                    ksi_changes(k,ksi_index)=ksi_changes(k,ksi_index)+ change*( errors==min(errors) )'/sum(errors==min(errors));%change*( errors==min(errors) ).';
                    %change*( errors==min(errors) )'/sum(errors==min(errors))
                    %change( errors==min(errors) )
                    %these commends are for writing this line in a
                    %different way. As it mentioned i could not find an
                    %elegant way
                end
            end
        end
    end
    
else % when reflection is more than 10 shadow boundaries are stabilized 
     % so no need for optimization. use ksi values that were used 2 
     %reflections before. we get the ksi changes data from "Main" fnc.
    
    % optimized ksi values of the 2 reflection before is applied to this
    % problem
    [intervalMatrix,ksi_on_use,~]=intervals_with_optimization(t1,t2,no_of_intervals,eta,shape_type,chan_of_var,stabilized_ksi);
    
    % explanization of the below part is done above
    
    j=1;
    for i=1:no_of_intervals
        basisfncs(j:j+max_poly_degree(i),:)...
            =Basisfncs(intervalMatrix(i,:),n,max_poly_degree(i),chan_of_var,ksi_on_use,eta,no_of_intervals,i,t1,t2,shape_type);
        j=j+max_poly_degree(i)+1;
    end
    longbasis=basisfncs.*repmat(e_term,sum(max_poly_degree+1),1);
    
    switch kerneltype
        case 1
            A=conj(longbasis)*kernel*longbasis.';
        case 2
            A=(conj(longbasis)*longbasis.')-(conj(longbasis)*kernel*longbasis.');
        case 3
            A=(conj(longbasis)*longbasis.')-(conj(longbasis)*kernel*longbasis.');
    end
    
    b=conj(longbasis)*rhs;
    
    S=eye(size(A)).*A;
    A=S^(-1)*A;
    b=S^(-1)*b;
    
    alfa=linsolve(A,b);
    
    soln=sum(longbasis.*repmat(alfa,1,2*n),1);
    
    compare=(soln-exact);
    y=log10( sqrt((compare*compare') / (exact*exact')) );
    
    errormin=y;
    ksi_changes=stabilized_ksi;
end

%display(ksi_changes)
%display(errormin)
%Picture( input, configuration_info, sb, shape_no, ksi_changes, no_of_ref )
% these parts are not necessary. they were for understanding the checking
% some things. "Picture" gives a nice picture about intervals.

galerkin_soln = soln;


end

