%% Forward-Backward algorithm 
clear all 
clc
clf

format short

% p = [0.2;0.3;0.1;0.4];
% a = [0,-1,1,1; 1,0,1,1; 1,1,0,1; 1,1,1,0];

% p = [0.3;0.3;0.4];
% a = [1,1,1;
%      -1,0,1;
%      1,0,0];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = abs(a);
a_pos = (A + a)./2;
a_neg = (A - a)./2;


% Initialization
[n,~] = size(a);
mu = ones(n,1);
nu = ones(n,1);
T = 200;
cost = [];

A_pos = a_pos.*exp(-1);
A_neg = a_neg.*exp(-1);


%% Iteration

for t = 1:T

%%%%%%%%%%%%%%%%%%
% update of nu
a_nu = A_pos'*(exp(-mu).*p);
b_nu = A_neg'*(exp(mu).*p);
c_nu = p;
% location indicator
nu_neg_loc = (b_nu~=0);
nu_pos_loc = ~nu_neg_loc;
% two cases 
nu_var = zeros(n,1);
nu_var(nu_neg_loc) = FBsolver(a_nu(nu_neg_loc),b_nu(nu_neg_loc),c_nu(nu_neg_loc));
nu_var(nu_pos_loc) = log(a_nu(nu_pos_loc)./c_nu(nu_pos_loc));
nu = nu_var;


%%%%%%%%%%%%%%%%%%
% update of mu
a_mu = A_pos*(exp(-nu));
b_mu = A_neg*(exp(nu));
c_mu = ones(n,1);
% location indicator
mu_neg_loc = (b_mu~=0);
mu_pos_loc = ~mu_neg_loc;
% two cases 
mu_var = zeros(n,1);
mu_var(mu_neg_loc) = FBsolver(a_mu(mu_neg_loc),b_mu(mu_neg_loc),c_mu(mu_neg_loc));
mu_var(mu_pos_loc) = log(a_mu(mu_pos_loc)./c_mu(mu_pos_loc));
mu = mu_var;


%%%%%%%%%%%%%%%%%%
% convergence 
P = diag(exp(-mu))*A_pos*diag(exp(-nu)) + diag(exp(mu))*A_neg*diag(exp(nu));

obj = p.*P.*log(P);
obj(isnan(obj)) = 0;
obj_value = sum(real(obj),'all');
cost = [cost, obj_value];

end



%%
Pi = a.*P;

%% Result checking -- prelimaries
disp('objective value:')
sum(obj_value,'all')


disp('stochastic matrix:')
disp(Pi'*p)
disp(p)
sum(Pi,2)

disp('Pi matrix:')
disp(Pi)
% % disp(P)
% disp('Ajacency matrix:')
% disp(a)

figure(1)
plot(1:T,cost);
% axis tight;
xlabel('Iteration') 
ylabel('Objective') 
title('Convergence')