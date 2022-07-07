function [z]=fold(X,f1,f2,f3)
% [f1*f3,f2]=size(X);
z=zeros(f1,f2,f3);
for i=1:f3
    z(:,:,i)=X((i-1)*f1+1:i*f1,:);
end


