function [MT]=Trans(M)
[m1,n1,k1]=size(M);
result=zeros(m1,m1,k1);
MT=zeros(n1,m1,k1);
for i=1:k1
    MT(:,:,i)=M(:,:,i)';
end