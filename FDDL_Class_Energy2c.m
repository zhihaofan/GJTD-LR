function [gap] = FDDL_Class_Energy2c(Ai,D,Xi,drls,trls,index,lambda1,classn,tau2,tau3)
%  ���ϵ�� Xi ��Ŀ�꺯��ֵ
%=======================================================================
% Class energy computation of FDDL, Version 1.0 
% Input :   (1) Ai :  the data matrix of this class
%           (2) D :   the whole dictionary
%           (3) Xi:   the coefficient matrix of this class
%           (4) Xa:   the coefficient matrix of the whole class
%           (5) drls: labels of dictionary's column
%           (6) trls: labels of training samples
%           (7) index: label of class being processed
%           (8) lambda1 : parameter of l1-norm energy of coefficient
%           (9) lambda2 : parameter of within-class scatter
%           (10) lambda3 : parameter of between-class scatter
%           (11) lambda4:  parameter of l2-norm energy of coefficient
%           (12) classn:   the number of class
%           (13) tau2: parameter of ||A_i-D_iX_i^i||_F^2 in fidelity term
%           (14) tau3: parameter of ||D_jX_i^j||_F^2 in  fidelity term
% 
% Outputs : (1) gap  :    the total energy of some class
%
%------------------------------------------------------------------------
 %Ai ��i����������  D�����ֵ�  Xi:��ǰ��i��������Ӧ��ϵ����xi 
 %drls:���ֵ��ǩ��1*700) trls:���������ı�ǩ  index����ǰ���ڴ����������
 %classn:�������    % fish_tau2:  parameter of ||D_jX_i^j||_F^2
                     % fish_tau3:  parameter of ||A_i-D_iX_i^i||_F^2
                     
GAP1  =   norm((Ai-D*Xi),'fro')^2;      %||A_i-DX_i||_F^2
GAP2  =   lambda1*sum(abs(Xi(:)));      %||Xi||1    
GAP5  =   (tau2)*norm((Ai-D(:,drls==index)*Xi(drls==index,:)),'fro')^2;    %||A_i-D_iX_i^i||_F^2
   
GAP6  =   0;                              %��||D_jX_i^j||_F^2 
    for i = 1:classn
        if i~=index
        GAP6 = GAP6+tau3*norm(D(:,drls==i)*Xi(drls==i,:),'fro')^2;
        end
    end
    gap = GAP1+GAP2+GAP5+GAP6;