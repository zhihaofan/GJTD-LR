clear
a = rand(5,5,3);
b = rand(5,5,3);
c = tprod(a,b);
f_a = fft(a,[],3);
f_b = fft(b,[],3);
f_c = fft(c,[],3);

%% -----------------------------1----------------------
a_res1=zeros(size(a));
tic;
halfn3 = round(size(a,3)/2); % round：四舍五入，若结果是.5，则向无穷
a_res1(:,:,1)=f_c(:,:,1)/f_b(:,:,1);
for i = 2 : halfn3
    a_res1(:,:,i) = f_c(:,:,i)/f_b(:,:,i);
    a_res1(:,:,size(a,3)+2-i) = conj(a_res1(:,:,i));
end
% for i = 1:size(a,3)
%     a_res1(:,:,i)=f_c(:,:,i)/f_b(:,:,i);
% end
a_res1 = ifft(a_res1,[],3);
toc;

%% -----------------------------2------------------------
tic;
a_res2=tprod(c,inverse(b));
toc;


