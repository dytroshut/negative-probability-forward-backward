%% Gradient descent for negative probability

clear all 
clc

% p = [0.2;0.3;0.1;0.4];
% 
% 
% a = [1,-1,1,1; 1,1,1,1; 1,1,1,1; 1,1,1,1];
% a = [0,-1,1,1; 1,0,1,1; 1,1,0,1; 1,1,1,0];

% p = [0.3;0.3;0.4];
% a = [1,1,1;
%      -1,0,1;
%      1,0,1];
% 
% A = abs(a);

p = [0.1;0.05;0.05;0.15;0.2;0.05;0.03;0.07;0.25;0.05];

a = [1  0 0 0 -1 0 1 0  1 0;
     0  0 1 0  0 0 0 0  0 1;
     0  1 0 0  1 0 1 0  0 0;
     0  0 0 0  1 0 0 1  0 0;
     -1 0 1 1  1 -1 0 0  0 0;
     0  0 0 0  -1 0 1 0  1 0;
     1  0 1 0  0 1 1 0  0 0;
     0  0 0 1  0 0 0 0  1 -1;
     1  0 0 0  0 1 0 1  0 1;
     0  1 0 0  0 0 0 -1 1 0];

% p = rand(10,1);
% p = p./sum(p);


% gradient descent
n = length(p); 

% sign-definite adjacency matrix
A = abs(a);

% algorithm intitialization
gamma = 0.05; % step size
T = 5000;   % number of iteration

mu = zeros(n,1);
lambda = zeros(n,1);

neg_A = (A-a)./2;
pos_A = (A+a)./2;
P = A;

% convergence check
cost = [];

A_pos = pos_A.*exp(-1);
A_neg = neg_A.*exp(-1);

%%
% gradient descent method iteration
for t = 1:T

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
end 

%%
Pi = a.*P;

%%
disp('objective value:')
sum(obj_value,'all')


disp('stochastic matrix:')
disp(Pi'*p-p)
sum(Pi,2)

disp('Pi matrix:')
disp(Pi)
disp(P)
disp('Ajacency matrix:')
disp(a)

figure()
plot(1:T,cost);