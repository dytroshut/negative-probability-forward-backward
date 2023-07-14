function [P,Pi,cost,Err_1,Err_2] = fb(Adj_pos,Adj_neg,p,iter)

%% The Forward-Backward algorithm  -- Anqi Dong, Tryphon T. Georgiou, Allen Tannenbaum
% Input:
% 1. Adj_pos: generated positive adjacence matrix (symmetric)
% 2. Adj_neg: generated random edges, assigned with negative $-1$
% 3. p: probability (column) vector p
% 4. iter: number of iterations.
%
% Output:
% 1. P: transition matrix
% 2. obj: vector of objective values
% 3. Err_1, Err_2: vector of marginal errors.


[num_node,~] = size(Adj_pos);

A = Adj_pos + Adj_neg;
A(A>=1) = 1;

P = A;
a = A;
a(Adj_neg~=0) = -1;

a_pos = (A + a)./2;
a_neg = (A - a)./2;


%% algorithm initialization

% intial value of variable \mu and \nu, set as all-ones vectors.
mu = ones(num_node,1);
nu = ones(num_node,1);

% element scaling by exp(-1).
A_pos = a_pos.*exp(-1);
A_neg = a_neg.*exp(-1);

% empty vector for iteration tracking.
cost = [];
Err_1 = [];
Err_2 = [];


%% forward-backward (sinkhorn-like) iteration
for t = 1:iter


%% 1. update of nu_i
% compute a,b,c for \nu
a_nu = A_pos'*(exp(-mu).*p);
b_nu = A_neg'*(exp(mu).*p);
c_nu = p;

% address the positive & negative elements, using location indicator
nu_neg_loc = (b_nu~=0);
nu_pos_loc = ~nu_neg_loc;

% update \nu_i in two cases, i.e.,
nu_var = zeros(num_node,1);
% 1) when a_{ij} < 0, compute \nu_i using function 'FB_solver'
nu_var(nu_neg_loc) = FBsolver(a_nu(nu_neg_loc),b_nu(nu_neg_loc),c_nu(nu_neg_loc));
% 2) when a_{ij} > 0, compute \nu_i using sinkhorn iteration
nu_var(nu_pos_loc) = log(a_nu(nu_pos_loc)./c_nu(nu_pos_loc));
% 3) update \nu from the two cases
nu = nu_var;



%% 2. update of mu_j
% compute a,b,c for \nu
a_mu = A_pos*(exp(-nu));
b_mu = A_neg*(exp(nu));
c_mu = ones(num_node,1);

% address the positive & negative elements, using location indicator
mu_neg_loc = (b_mu~=0);
mu_pos_loc = ~mu_neg_loc;

% update \mu_j in two cases, i.e.,
mu_var = zeros(num_node,1);
% 1) when a_{ij} < 0, compute \mu_j using function 'FB_solver'
mu_var(mu_neg_loc) = FBsolver(a_mu(mu_neg_loc),b_mu(mu_neg_loc),c_mu(mu_neg_loc));
% 2) when a_{ij} > 0, compute \mu_j using sinkhorn iteration
mu_var(mu_pos_loc) = log(a_mu(mu_pos_loc)./c_mu(mu_pos_loc));
% 3) update \nu from the two cases
mu = mu_var;


%% 3. convergence 
% update \P according to the element-wise closed-form expression.
P = diag(exp(-mu))*A_pos*diag(exp(-nu)) + diag(exp(mu))*A_neg*diag(exp(nu));

% compute and save the value of the objective function.
obj = p.*P.*log(P);
obj(isnan(obj)) = 0;
obj_value = sum(real(obj),'all');
cost(end+1) = obj_value;

% compuate and save the violation/error of the two marginals (marginal constraints).
Pi = a.*P;
Err_1(end+1) = norm(Pi'*p-p)/norm(p);
Err_2(end+1)=norm(Pi*ones(num_node,1)-ones(num_node,1))/norm(ones(num_node,1));



end







end