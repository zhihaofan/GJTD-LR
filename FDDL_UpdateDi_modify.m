function [Di_new,delet] = FDDL_UpdateDi_modify (A1,X1,index,classn,Fish_ipts,Fish_par)
% ========================================================================
% Dictionary updating of FDDL, Version 1.0
% Copyright(c) 2011  Meng YANG, Lei Zhang, Xiangchu Feng and David Zhang
% All Rights Reserved.
% ----------------------------------------------------------------------- 
% This is an implementation of the algorithm for updating the
% dictionary of FDDL (fix the coefficient)
% Please refer to the following paper
% Meng Yang, Lei Zhang, Xiangchu Feng, and David Zhang,"Fisher Discrimination 
% Dictionary Learning for Sparse Representation", In IEEE Int. Conf. on
% Computer Vision, 2011.
% Meng Yang, Lei Zhang, Jian yang and David Zhang, "Metaface learning for 
% sparse representation based face recognition", In ICIP, 2010.
%----------------------------------------------------------------------
%  Inputs :   (1) A :    the training data
%             (2) X :    the coefficient matrix of the training data
%             (3) index:   the label of the class being processed
%             (4) classn:  the total number of classes
%             (5) Fish_ipts
%                        . D  the dictionary in the last interation
%                        . trls  the labels of training data
%             (6) Fish_par
%                        .dls   the labels of the dicitonary atoms
% Outputs:    (1) Di_new :  the updated dictionary of the index-th class
%             (2) delet:    the indication of deleted atoms
%--------------------------------------------------------------------
k1=size(A1,3);
Dit=zeros(size(Fish_ipts.D,1),size(Fish_ipts.D,2)/classn,size(Fish_ipts.D,3));
for ii=1:k1
    D    =  Fish_ipts.D(:,:,ii);
    A=A1(:,:,ii);
    X=X1(:,:,ii);
    tau3 =  1;
    tau2 =  1;
    trls =  Fish_ipts.trls;
    dls =  Fish_par.dls;
    delet = [];

    Do = D(:,dls~=index);
    Di = D(:,dls==index);
    Ai = A(:,trls==index);
    Ao = A(:,trls~=index);
    Xi = X(:,trls==index);
    Xo = X(:,trls~=index);

    Xi_i = Xi(dls==index,:);
    Xi_o = Xi(dls~=index,:);
    Xo_i = Xo(dls==index,:);
    Xo_o = Xo(dls~=index,:);
    % X_i  = X(dls==index,:);

    Zi = Ai-Do*Xi_o;
    Zo = Ao-Do*Xo_o;

    for i = 1: size(Di,2)
        Yi = Zi - Di*Xi_i + Di(:,i)*Xi_i(i,:);
        Ui = Ai - Di*Xi_i + Di(:,i)*Xi_i(i,:);
    %     Ua = A  - Di*X_i+Di(:,i)*X_i(i,:);
        Vo = Zo - Di*Xo_i + Di(:,i)*Xo_i(i,:);

        UaXti = zeros(size(Yi*(Xi_i(i,:))'));
        for t_i = 1:classn
            if t_i~=index
            Xt_i = X(dls==index,trls==t_i);
            UaXti = UaXti + (0  - Di*Xt_i + Di(:,i)*Xt_i(i,:))*(Xt_i(i,:))';
            end
        end

        tem1 = -Yi*(Xi_i(i,:))' - (tau2)*Ui*(Xi_i(i,:))' - Vo*(Xo_i(i,:))' - tau3*UaXti;
    %     tem1 = - (tau2)*Ui*(Xi_i(i,:))';
    %     tem1 = -Yi*(Xi_i(i,:))' - (tau2)*Ui*(Xi_i(i,:))'  - tau3*Ua*(Xo_i(i,:))';
        tem  = -tem1;        %norm(a,2) 返回A的最大奇异值，若Di 字典原址的最大奇异值小于10的-6次方，将此列原子置为0
    %     if norm(tem,2)<1e-6
    %         Di(:,i) = zeros(size(tem));
    %         delet = [delet i];
    %     else
            Di(:,i) = tem./norm(tem,2);   
    %     end
    end
    Dit(:,:,ii)=Di;
end






Di_new = Dit;