function [im]=is_nan(InImg)

InImg=real(InImg);

aa=0;
[m1,m2,m3]=size(InImg);
for i=1:m1
    for j=1:m2
        for k=1:m3
            if isnan(InImg(i,j,k))
%                 aa=aa+1
                InImg(i,j,k)=InImg(i,j-1,k);
%                 i
%                 j
%                 k
            end
        end
    end
end
% aa
im=InImg;