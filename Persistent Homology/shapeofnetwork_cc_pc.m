function [MST, A, Dx, Cx] = shapeofnetwork_cc_pc(X,VAR)

% X : p*n matrix (p: # of ROIs, n: # of samples)
% MST : minimum spanning tree (p-1)*3 matrix (each row contains [i j r]
% where i and j are connected with the distance r. Total (p-1) edges are
% exist in MST. 
% A : p*p adjacency matrix of MST
% Dx : p*p single linkage matrix %(sqrt)HSG
% Cx : p*p distance matrix %(sqrt)HSG 

[n p] = size(X); %[p n]= size(X); 

%C = partialcorr(X,VAR)  %C = corrcoef(X);  %C=partialcorr(X',VAR);
C = corr(X,'Type','Pearson');  
%C = partialcorr(X',VAR); 
%C = corr(X,'Type','Spearman');
%C = corrcoef(X');

C = C.*(C>0); % Find only the positive matrix and compute the distance
Cx = sqrt(1 - C); % added sqrt command

[Tree, pred] = graphminspantree(sparse(Cx),'Method','Kruskal');
[row col] = find(Tree);
val = nonzeros(Tree);
MST = [row col val];
[tval tind] = sort(val,'ascend');
MST = MST(tind,:);
A = full(Tree);
% Dx = single_linkage_distance(MST,size(Cx,1));
Dx = single_linkage_distance(Cx); 


