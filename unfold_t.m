function [B]=unfold_t(X)
B=[];
[a,b,c]=size(X);
% A=cell(c,1);
for i=1:c
%     B{i}=X(:,:,i);
    B=[B;X(:,:,i)];
end
% B=cell2mat(A);