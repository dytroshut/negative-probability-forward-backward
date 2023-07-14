%% https://www.mathworks.com/matlabcentral/fileexchange/48576-circulargraph
%% graph plot and convergence comparison

clear all 
clc
clf
format short

%% Example 3 

% num_node = 10;
% 
% p = [0.1;0.05;0.05;0.15;0.2;0.05;0.03;0.07;0.25;0.05];
% 
% a = [1  0 0 0 -1 0 1 0  1 0;
%      0  0 1 0  0 0 0 0  0 1;
%      0  1 0 0  1 0 1 0  0 0;
%      0  0 0 0  1 0 0 1  0 0;
%      -1 0 1 1  1 -1 0 0  0 0;
%      0  0 0 0  -1 0 1 0  1 0;
%      1  0 1 0  0 1 1 0  0 0;
%      0  0 0 1  0 0 0 0  1 -1;
%      1  0 0 0  0 1 0 1  0 1;
%      0  1 0 0  0 0 0 -1 1 0];
% 
% Adj_pos = (a>0);
% Adj_neg = (a<0);

%% Random generator using functions: graph_generator,negative_generator


% num_node = 100;
% [Adj_pos,Gp] = graph_generator(num_node, 0.17);
% Adj_neg = negative_generator(num_node,0.07);
% 
% sum(Adj_pos,'all')
% sum(Adj_neg,'all')
% 
% p = rand(1,num_node);
% p = p'./sum(p);
% 
% %% Data saving
% 
% writematrix(p,'pvector');
% writematrix(Adj_pos,'Ap');
% writematrix(Adj_neg,'Am');

%%

num_node = 100;
p = readmatrix('pvector');
Adj_pos = readmatrix('Ap');
Adj_neg = readmatrix('Am');

sum(Adj_pos,'all')
sum(Adj_neg,'all')


%% 1. Forward-Backward algorithm 

T_fb = 50;
[P_fb,Pi_fb,cost_fb,Err_1_fb,Err_2_fb] = fb(Adj_pos,Adj_neg,p,T_fb);


%% Result checking -- prelimaries
disp('..............Forward-Backward.............')
disp('objective value:')
disp(cost_fb(end))

% disp('stochastic matrix:')
disp(Pi_fb'*p-p);
sum(Pi_fb,2);
sum(Pi_fb,'all')

figure(1);
plot(1:T_fb,cost_fb);
axis tight;
xlabel('Iteration') 
ylabel('Objective') 
title('Convergence')

figure(2);
subplot(2,1,1)
con1 = plot(1:T_fb,cost_fb); 
% ylabel('Objective value','FontSize',10);
% xlabel('iterations','FontSize',15);
ylim([-2.94 -2.91]);
con1.Color = "#0072BD";
con1.LineWidth = 2;
set(gca,'fontsize',12);
subplot(2,1,2)
con2 = plot(1:T_fb,log(Err_1_fb));
axis tight; 
% ylabel('Marginal','FontSize',15);
xlabel('iterations','FontSize',15)
con2.Color = "#A2142F";
con2.LineWidth = 2;
set(gca,'fontsize',12);
% figure(2);
% subplot(2,1,1)
% plot(1:T_fb,log10(Err_1_fb));axis tight; %title('log|\pi^T p - p|');
% subplot(2,1,2)
% plot(1:T_fb,log(Err_2_fb));axis tight; %title('log|\pi 1 - 1|');


%% Gene network plot

% A_gene = (P_fb~=0);
% 
% g = graph(A_gene);
% 
% [neg_row,neg_col] = find(Pi_fb < 0 );
% redEdge = [neg_row, neg_col];
% colormap = lines(num_node);
% % colormap = linspecer(num_node,'sequential');
% figure(3)
% gplot = plot(g,'Layout','circle');
% gplot.NodeColor = colormap;
% gplot.MarkerSize = 10;
% gplot.NodeFontSize = 14;
% gplot.Marker = 'o';
% gplot.LineWidth = 0.8;
% % highlight(gplot,neg_row,neg_col,'MarkerSize',15,'NodeColor','red','Marker','diamond','EdgeColor','r','LineWidth',3)
% highlight(gplot,neg_row,neg_col,'EdgeColor','r','LineWidth',2.5)


%% Vertification by gradient descent

% T_gd = 100000;
% gamma = 0.00001;

% T_gd = 3000;   % number of iteration
% gamma = 0.05; % step size
% 
% [P_gd,Pi_gd,cost_gd,Err_1_gd,Err_2_gd] = gd(Adj_pos,Adj_neg,p,gamma,T_gd);


%% gradient descent plots

% disp('..............Gradient Descent.............')
% disp('objective value:')
% disp(cost_gd(end))
% disp('stochastic matrix:')
% disp(Pi_gd'*p-p);
% sum(Pi_gd,2);
% sum(Pi_gd,'all')
% 
% figure(4);
% plot(1:T_gd,cost_gd);
% axis tight;
% xlabel('Iteration') 
% ylabel('Objective') 
% title('Convergence')

% figure(5);
% subplot(2,1,1)
% plot(1:T_gd,log10(Err_1_gd));axis tight; %title('log|\pi p - p|');
% subplot(2,1,2)
% plot(1:T_gd,log(Err_2_gd));axis tight; %title('log|\pi 1 - 1|');









