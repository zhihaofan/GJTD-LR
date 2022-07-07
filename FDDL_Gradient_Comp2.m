function grad = FDDL_Gradient_Comp2(Xi,classn,index,tau2,tau3,trls,drls,newpar)
% ========================================================================
% IPM's Gradient computation of FDDL, Version 1.0
% -----------------------------------------------------------------------  
% Input :   (1) Xi:   the coefficient matrix of this class
%           (2) Xa:   the coefficient matrix of the whole class
%           (3) classn:   the number of class
%           (4) index: label of class being processed
%           (8) tau2: parameter of ||A_i-D_iX_i^i||_F^2 in fidelity term
%           (9) tau3: parameter of ||D_jX_i^j||_F^2 in  fidelity term
%           (10) trls: labels of training samples
%           (11) drls: labels of dictionary's column      
% Outputs : (1) grad  :    the gradient vector of coding model
%
%------------------------------------------------------------------------
n_d             =      newpar.n_d;            % the sample number of i-th training data
DD              =      newpar.DD;
DAi             =      newpar.DAi;
Di0Di0          =      newpar.Di0Di0;
Di0Ai           =      newpar.Di0Ai;
m               =      newpar.m;
DoiDoi          =      newpar.DoiDoi;     % %Doi'*Doi ΪDj'*Dj 

XiT      =   Xi';
tem      =   2*DD*Xi-2*DAi;      % 700*7  || Ai -DXi ||F^2 ��Xi ��  
grad1    =   tem(:);             %4900*1     || Ai -DXi ||F^2 ���ֵĵ�����

t1 = (tau2)*2*Di0Di0*Xi-(tau2)*2*Di0Ai;      %|| Ai -DiXi_i ||F^2 ���� ��Xi�ĵ�����
t2 = DoiDoi;
t3 = DoiDoi*Xi;    %Doi'*Doi ΪDj'*Dj    ����i��֮�⣬������ ||D_jX_i^j||_F^2 ���� ��Xi ����
t4 = tau3*2*DoiDoi*Xi;   %tau3: parameter of ||D_jX_i^j||_F^2 in  fidelity term
tem      =   tau3*2*DoiDoi*Xi+(tau2)*2*Di0Di0*Xi-(tau2)*2*Di0Ai;
grad7    =   tem(:);
%tau2: parameter of ||A_i-D_iX_i^i||_F^2 in fidelity term
grad = grad1+grad7;   %������fisher ���־��Ƕ�ǰ����ĵ��� ��grad1+grad7)
end