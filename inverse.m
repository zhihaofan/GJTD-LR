function [result]=inverse(MMT)
result=zeros(size(MMT));
fMMT=fft(MMT,[],3);
k1=size(MMT,3);
for j=1:k1
    result(:,:,j)=pinv(fMMT(:,:,j));
end
result=ifft(result,[],3);

% hyper2=hyperConvert2D(MMT);
% 
% ihyper=pinv(hyper2);
% result=hyperConvert3D(ihyper',size(MMT,1),size(MMT,2));