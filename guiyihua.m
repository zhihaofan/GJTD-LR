function OutImg = guiyihua(InImg,maxim,minim)
InImg=abs(InImg);
[m1,m2,m3]=size(InImg);
% for i=1:m1
%     for j=1:m2
%         for k=1:m3
%             if InImg(i,j,k)<=0.1
%                 InImg(i,j,k)=2*rand(1);
%             end
%             if isnan(InImg(i,j,k))
%                 InImg(i,j,k)=0;
%             end
%         end
%     end
% end


ymax=maxim;ymin=minim;
xmax = max(max(max(InImg)));
xmin = min(min(min(InImg)));
OutImg = (ymax-ymin)*(InImg-xmin)/(xmax-xmin) + ymin; 

end
