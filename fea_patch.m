function [X_m]=fea_patch(b,lImfea_Z,r,c)

k1=size(lImfea_Z,4);
k2=size(lImfea_Z,3);

X_m=zeros(b*b*k2,size(r,2)*size(c,2),k1);

for kk=1:k1
    lImG11=lImfea_Z(:,:,1,kk);
    lImG12=lImfea_Z(:,:,2,kk);
    lImG21=lImfea_Z(:,:,3,kk);
    lImG22=lImfea_Z(:,:,4,kk);
    k    =  0;
    for i  = 1:b
        for j  = 1:b
            k       =  k+1;
            blk11     = lImG11(r-1+i,c-1+j);
            Px11(k,:) =  blk11(:)';
        end
    end
    
    k    =  0;
    for i  = 1:b
        for j  = 1:b
            k       =  k+1;
            blk12     = lImG12(r-1+i,c-1+j);
            Px12(k,:) =  blk12(:)';
        end
    end
    
    k    =  0;
    for i  = 1:b
        for j  = 1:b
            k       =  k+1;
            blk21     = lImG21(r-1+i,c-1+j);
            Px21(k,:) =  blk21(:)';        
        end
    end

    k    =  0;
    for i  = 1:b
        for j  = 1:b
            k       =  k+1;
            blk22    = lImG22(r-1+i,c-1+j);
            Px22(k,:) =  blk22(:)';        
        end
    end
    temp=[Px11;Px12;Px21;Px22];
    X_m(:,:,kk)=temp;

end



