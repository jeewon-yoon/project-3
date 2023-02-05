function [MST,mtx_MST] = plot_minimum_spanning_tree(D,color_list) 
% Plot minimum spanning tree 
% This function calls 'arch_layout.m' for estimating the location of nodes
% in a tree. 
% D : p-by-p distance matrix 

% number of nodes 
p = size(D,1); 


% Estimate the minimum spanning tree 
[Tree,pred] = graphminspantree(sparse(D)); 
[row,col] = find(Tree); 
val = Tree(find(Tree)); 
MST = [row col val]; 
[tval, tind] = sort(val,'ascend'); 
MST = full(MST(tind,:)); 
mtx_MST = full(Tree+Tree'); 

% Find the position of nodes using TreeVis 
loc_node = arch_layout(full(Tree+Tree')); 
% Find the position of nodes using Kamada-Kawai algorithm 
% cd('.\matlab_bgl'); 
% [X,data] = kamada_kawai_spring_layout((Tree+Tree'));
% loc_node = X'; 
% cd .. 

% Color of nodes 
if nargin < 2, 
    color_list = colormap(jet(p)); 
end 

% Draw edgs
for i = 1:p-1, 
    line(loc_node(1,MST(i,1:2)),loc_node(2,MST(i,1:2)),'Color', ... 
        mean(color_list(MST(i,1:2),:),1),'LineWidth',2); 
    hold on; 
end 

% Draw nodes
scatter(loc_node(1,:),loc_node(2,:),30,color_list,'fill'); 

