function [P,Pi,cost,Err_1,Err_2] = gd(Adj_pos,Adj_neg,p,gamma,iter)


%% The gradient descent method  
% Input:
% 1. Adj_pos: generated positive adjacence matrix (symmetric)
% 2. Adj_neg: generated random edges, assigned with negative $-1$
% 3. p: probability (column) vector p
% 4. iter: number of iterations.
% 5. gamma: fixed step size
%
% Output:
% 1. P: transition matrix
% 2. obj: vector of objective values
% 3. Err_1, Err_2: vector of marginal errors.
%
% @Copyright: Anqi Dong, Tryphon T. Georgiou, Allen Tannenbaum

%% Adjancecy matrix construction

[num_node,~] = size(Adj_pos);

A = Adj_pos + Adj_neg;
A(A>=1) = 1;

P = A;
a = A;

a(Adj_neg~=0) = -1;

a_pos = (A + a)./2;
a_neg = (A - a)./2;


%% Initializations

mu = zeros(num_node,1);
lambda = zeros(num_node,1);

A_pos = a_pos.*exp(-1);
A_neg = a_neg.*exp(-1);

cost = [];
Err_1 = [];
Err_2 = [];



%% 
for t = 1:iter

delta_mu = (a.*P)'*p - p;
mu = mu + gamma.*delta_mu;

P = A_pos.*(exp(-lambda./p)*exp(-mu)') + A_neg.*(exp(lambda./p)*exp(mu)');

delta_lambda = sum(a.*P,2) - 1;
lambda = lambda + gamma.*delta_lambda;

P = A_pos.*(exp(-lambda./p)*exp(-mu)') + A_neg.*(exp(lambda./p)*exp(mu)');


obj = p.*P.*log(P);
obj(isnan(obj)) = 0;
obj_value = sum(real(obj),'all');
cost(end+1) = obj_value;

Pi = a.*P;
Err_1(end+1) = norm(Pi'*p-p)/norm(p);
Err_2(end+1) = norm(Pi*ones(num_node,1)-ones(num_node,1))/norm(ones(num_node,1));
end

end 


