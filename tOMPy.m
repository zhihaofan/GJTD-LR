function [X]=tOMPy(D,Y,L)

[m1,n1,k1]=size(D);
[m2,n2,k2]=size(Y);

X=zeros(n1,n2,k1);
for i=1:k1
    X(:,:,i)=OMP(D(:,:,i),Y(:,:,i),L);
end
