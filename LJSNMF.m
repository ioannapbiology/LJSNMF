% set the working working directory
%cd('/home/ioanna/Documents/first_epoch/LJSNMF/scripts/lungCancerTCGA')

function [result_labels] = LJSNMF(A1, A2, A3, A4, W, D, L, beta, gamma, delta, ground_truth_labels)

 %% set some variables
n = size(A1,1);
k = 2;
max_iter = 1000;
TolFun = 1e-5;
obj_vector_values = nan(max_iter,1);

% initialize H, Hs
H1 = rand(n,k);
H2 = rand(n,k);
H3 = rand(n,k);
H4 = rand(n,k);
Hs = rand(n,k); %Hstar

% the Q
Q1 = diag(sum(H1,1));
Q2 = diag(sum(H2,1));
Q3 = diag(sum(H3,1));
Q4 = diag(sum(H4,1));

% ============= OPTIMIZATION PART ==========================

for i=1:max_iter
    
    % update rules part -------------
    
    % update rule for H1
    Q1 = diag(sum(H1,1));
    numer1 = A1*H1 + gamma{1}*(Hs*Q1');
    H1 = H1 .*  (1 - delta + delta*(numer1./...
        ((H1*(H1'*H1)) + gamma{1}*((H1*Q1)*Q1') + eps(numer1))));


    % update rule for H2
    Q2 = diag(sum(H2,1));
    numer2 = A2*H2 + gamma{2}*(Hs*Q2');
    H2 = H2 .*  (1 - delta + delta*(numer2 ./...
        ((H2*(H2'*H2)) + gamma{2}*((H2*Q2)*Q2') + eps(numer2))));

    
    % update rule for H3
    Q3 = diag(sum(H3,1));
    numer3 = A3*H3 + gamma{3}*(Hs*Q3');
    H3 = H3 .*  (1 - delta + delta*(numer3 ./...
        ((H3*(H3'*H3)) + gamma{3}*((H3*Q3)*Q3') + eps(numer3))));
    
    % update rule for H4
    Q4 = diag(sum(H4,1));
    numer4 = A4*H4 + gamma{4}*(Hs*Q4');
    H4 = H4 .*  (1 - delta + delta*(numer4./...
        ((H4*(H4'*H4)) + gamma{4}*((H4*Q4)*Q4') + eps(numer4))));

    
    % update rule  Hstar (H*)
    numers = gamma{1}*H1*Q1 + gamma{2}*H2*Q2 + gamma{3}*H3*Q3 + gamma{4}*H4*Q4 + beta*W*Hs;
    Hs = Hs .*  (numers ./...
        ((gamma{1}*Hs + gamma{2}*Hs + gamma{3}*Hs + gamma{4}*Hs)...
        + beta*D*Hs  + eps(numers)));
    
    % end of update rules part ------
    
    
    %----------------------------------------------------------------------------------------
    % calculate the value of the objective function, store it in
    % 'obj_vector_values',  for every iteration 'i'
    
    p1 = (A1 - H1*H1');
    p2 = (A2 - H2*H2');
    p3 = (A3 - H3*H3');
    p4 = (A4 - H4*H4');
    
    p5 = (H1*Q1 - Hs);
    p6 = (H2*Q2 - Hs);
    p7 = (H3*Q3 - Hs);
    p8 = (H4*Q4 - Hs);
    
    p9 = (Hs'*L)*Hs;
    
    
    obj_vector_values(i) = (sum(p1(:).^2)) + (sum(p2(:).^2)) + (sum(p3(:).^2)) + (sum(p4(:).^2))+... % fix it!!!
        gamma{1}*(sum(p5(:).^2)) + gamma{2}*(sum(p6(:).^2)) + gamma{3}*(sum(p7(:).^2)) + gamma{4}*(sum(p8(:).^2)) +...
        beta*trace(p9);
    
    
    %----------------------------------------------------------------------------------------
    
          
    if mod(i,10) == 0 || mod(i+1,10) == 0
        
        obj = obj_vector_values(i);
        
        if mod(i+1,10) == 0
            obj0=obj;
            continue
        end
        
        if mod(i,100) == 0
            display(sprintf('...LJ-SNMF iteration #%d out of %d, error: %f\n',...
               i, max_iter, obj ));
        end
        
        if exist('obj0')
            assert(obj <= obj0, sprintf('Objective Function Increasing! From %f to %f.',...
                obj0, obj));
        end
        
        %check for convergence
        if exist('obj0') && obj0 - obj <= TolFun*max(1,obj0)
            display(sprintf('Stopped at iteration #%d, obj: %f, obj0: %f',  i, obj, obj0));
            break
        end
        
    end
    
    
end

%% =================== PLOTS ====================

plot(obj_vector_values, 'LineWidth',2)

title(['LJ-SNMF: Objective Value:' num2str(obj_vector_values(i))...
    ', β = ' num2str(beta)  ', γ =' num2str(gamma{1}) ', δ = ' num2str(delta) ])


%% =================  EVALUATION METRICS ============

% get the raw result labels from Hs
[~, result_labels] = max(Hs, [], 2);

result = result_labels;

% get the accuracy
AC = accuracy(ground_truth_labels, result_labels);

% get the Normalized Mutual Information
nmi_val = nmi(ground_truth_labels, result_labels);

% get the Adjusted Rand Index
Adj_RI = rand_index(ground_truth_labels, result_labels);

%print the evulation metrics results
fprintf('AC: %0.4f\tNMI:%0.4f\tAdj_RI: %0.4f\n', AC, nmi_val, Adj_RI)

end
