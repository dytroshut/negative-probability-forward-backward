%% Gradient Descent method:

clear all 
clc
clf
format short

addpath('function')
addpath('example_data')

%% Four examples are included and solved them by Gradient Descent method

%  Through out, the variables are defined as follows:
%  num_node: number of nodes of the network/graph  
%  p: the probability (column) vector 
%  A: sign-definite {0,1}-adjacency matrix
%  a: sign-indefinite {0,-1,1}-adjacency matrix
%  Adj_pos: location indicator of positive elements {1} of $a$
%  Adj_neg: location indicator of negative elements {-1} of $a$
%  gamma: step size used in the gradient descent

%  Moreover, we define the outputs:
%  P: the sign-definite transition matrix, i.e., P_{ij}>0
%  Pi=a.*P: sign-indefinite transition matrix


%% Case 1: 3-node network %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% num_node = 3;
% 
% p = [0.3; 0.3; 0.4];
% a = [ 1, 1, 1;
%      -1, 0, 1;
%       1, 0, 1];
% 
% Adj_pos = (a>0);
% Adj_neg = (a<0);


%% Case 2: 4-node network %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% num_node = 4;
% 
% p = [0.2; 0.3; 0.1; 0.4];
% a = [1,-1,1,1; 
%     -1, 1,1,1; 
%      1, 1,1,1; 
%      1, 1,1,1];
% 
% Adj_pos = (a>0);
% Adj_neg = (a<0);


%% Case 3: 10-node network %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_node = 10;

p = [0.1;0.05;0.05;0.15;0.2;0.05;0.03;0.07;0.25;0.05];

a = [1  0 0 0 -1  0 1 0  1  0;
     0  0 1 0  0  0 0 0  0  1;
     0  1 0 0  1  0 1 0  0  0;
     0  0 0 0  1  0 0 1  0  0;
    -1  0 1 1  1 -1 0 0  0  0;
     0  0 0 0 -1  0 1 0  1  0;
     1  0 1 0  0  1 1 0  0  0;
     0  0 0 1  0  0 0 0  1 -1;
     1  0 0 0  0  1 0 1  0  1;
     0  1 0 0  0  0 0 -1 1  0];

Adj_pos = (a>0);
Adj_neg = (a<0);


%% Case 4 (option 1): 100-node random network %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Random generator using functions: graph_generator,negative_generator

% num_node = 100;
% [Adj_pos,Gp] = graph_generator(num_node, 0.17);
% Adj_neg = negative_generator(num_node,0.07);
% 
% disp('number of positive edges')
% sum(Adj_pos,'all')
% disp('number of positive edges')
% sum(Adj_neg,'all')
% 
% p = rand(1,num_node);
% p = p'./sum(p);


% %% Data saving

% writematrix(p,'pvector');
% writematrix(Adj_pos,'Ap');
% writematrix(Adj_neg,'Am');


%% Case 4 (option 2): Reproduce the example in paper %%%%%%%%%%%%%%%%%%%%%

% num_node = 100;
% p = readmatrix('example_data/pvector');
% Adj_pos = readmatrix('example_data/Ap');
% Adj_neg = readmatrix('example_data/Am');
% 
% disp('number of positive edges')
% sum(Adj_pos,'all')
% disp('number of positive edges')
% sum(Adj_neg,'all')




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gradient descent
%% Attention: 
% number of iteration and the step size needs a case-by-case adjust.  

T = 3000;   % number of iteration
gamma = 0.05; % step size

[P,Pi,cost,Err_1,Err_2] = gd(Adj_pos,Adj_neg,p,gamma,T);


%% Result & Constraints check

disp('..............Gradient-Descent.............')

disp('..............objective value')
disp(cost(end))

disp('.....Marginal constraint and totol sum')
% disp(Pi'*p-p)
% sum( Pi,2 )
sum(Pi,'all')


figure(1);
plot(1:T,cost);
axis tight;
xlabel('Iteration') 
ylabel('Objective') 
title('Convergence')

figure(2);
subplot(2,1,1)
plot(1:T,log10(Err_1));axis tight; title('log||\Pi^T p - p||');
subplot(2,1,2)
plot(1:T,log(Err_2));axis tight; title('log||\Pi 1 - 1||');


%% Network Visualization (symmetric adjancenct matrix A is required)

A_graph = (P~=0);
g = graph(A_graph);

[neg_row,neg_col] = find(Pi < 0 );
redEdge = [neg_row, neg_col];
colormap = lines(num_node);
% colormap = linspecer(num_node,'sequential');
figure(3)
gplot = plot(g,'Layout','circle');
gplot.NodeColor = colormap;
gplot.MarkerSize = 10;
gplot.NodeFontSize = 14;
gplot.Marker = 'o';
gplot.LineWidth = 0.8;
highlight(gplot,neg_row,neg_col,'EdgeColor','r','LineWidth',2.5)


