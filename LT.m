function  [Z] =  LT(X,sf,sz)
F=create_F();
% sz=[512 512];
% sz=[64 64];
s0=1;
psf        =    fspecial('gaussian',7,2);
par.fft_B      =    psf2otf(psf,sz);
par.fft_BT     =    conj(par.fft_B);
par.H          =    @(z)H_z(z, par.fft_B, sf, sz,s0 );
par.HT         =    @(y)HT_y(y, par.fft_BT, sf, sz,s0);
par.P=create_F();
% F=F(:,3:31);
for band = 1:size(F,1) 
    div = sum(F(band,:));
    for i = 1:size(F,2)
        F(band,i) = F(band,i)/div;
    end
end
S=X;
[M,N,L] = size(S);
S_bar = hyperConvert2D(S);
hyper= par.H(S_bar);
Y_h = hyperConvert3D(hyper, M/sf, N/sf);
Y = hyperConvert3D((F*S_bar), M, N);
para.K=100;
para.eta=1e-2;
HSI=Y_h;
MSI=Y;
R=F;
K=para.K;
C=para.eta;
FBm=par.fft_B;

mu=1e-3;
patchsize=8;
overlap=4;
bparams.block_sz = [patchsize, patchsize];
bparams.overlap_sz=[overlap overlap];
[nr, nc,~]=size(MSI);
L=size(HSI,3);
num1=(nr-patchsize)/(patchsize-overlap)+1;
num2=(nc-patchsize)/(patchsize-overlap)+1;
bparams.block_num=[num1 num2]; 

fkmeans_omaxpt.careful = 1;
predenoised_blocks = ExtractBlocks(MSI, bparams);
Y2=Unfold(predenoised_blocks,size(predenoised_blocks),4);   
[aa ]=fkmeans(Y2,K,fkmeans_omaxpt);


HSI_int=zeros(nr,nc,L);
HSI_int(1:sf:end,1:sf:end,:)=HSI;

FBmC  = conj(FBm);
FBs  = repmat(FBm,[1 1 L]);

FBCs1=repmat(FBmC,[1 1 L]);
HHH=ifft2((fft2(HSI_int).*FBCs1));
HHH1=hyperConvert2D(HHH);
psfY.w=nr/sf;
psfY.h=nc/sf;
psfY.W=nr;
psfY.H=nc;

Zt=HSI_int;
HR_load1=imresize(HSI, sf,'bicubic');
ZE=hyperConvert2D(HR_load1);

%% iteration

MSI3=Unfold(MSI,size(MSI),3);
n_dr=nr/sf;
n_dc=nc/sf;

V1=ZE;
V2=ZE;
V3=ZE;
G1=zeros(size(V1));
G2=zeros(size(V2));
G3=zeros(size(V3));
CCC=R'*MSI3+HHH1;
C1=R'*R+3*mu*eye(size(R,2)); 
[Q,Lambda]=eig(C1);                        
Lambda=reshape(diag(Lambda),[1 1 L]);
InvLbd=1./repmat(Lambda,[ sf*n_dr  sf*n_dc 1]);   
B2Sum=PPlus(abs(FBs).^2./( sf^2),n_dr,n_dc);
InvDI=1./(B2Sum(1:n_dr,1:n_dc,:)+repmat(Lambda,[n_dr n_dc 1]));

for i=1:3
    HR_HSI3=mu*(V1+G1/(2*mu)+V2+G2/(2*mu)+V3+G3/(2*mu));
    C3=CCC+HR_HSI3;
    C30=fft2(reshape((Q\C3)',[nr nc L   ])).*InvLbd;
    temp  = PPlus_s(C30/( sf^2).*FBs,n_dr,n_dc);        
    invQUF = C30-repmat(temp.*InvDI,[ sf  sf 1]).*FBCs1;
    VXF    = Q*reshape(invQUF,[nc*nc L])';
    ZE = reshape(real(ifft2(reshape(VXF',[nr nc L]))),[nc*nc L])'; 

    %   [ZE1] = Sylvester(C1,psfY.B, sf,n_dr,n_dc,C3);  

    Zt=hyperConvert3D(ZE,nr, nc );
    % aa1=C1*ZE-C3+ creat_HSI_T(creat_HSI(ZE ,psfY),psfY);
    % norm(aa1(:))

    %% spatial low rank
    B1=ZE-G1/(2*mu);
    B1=hyperConvert3D(B1,nr, nc );
    predenoised_blocks1 = ExtractBlocks(B1, bparams);
    Z_block1=zeros( bparams.block_sz(1), bparams.block_sz(2),L, bparams.block_num(1)* bparams.block_num(2));


    B2=ZE-G2/(2*mu);
    B2=hyperConvert3D(B2,nr, nc );
    predenoised_blocks2 = ExtractBlocks(B2, bparams);
    Z_block2=zeros( bparams.block_sz(1), bparams.block_sz(2),L, bparams.block_num(1)* bparams.block_num(2));

    B3=ZE-G3/(2*mu);
    B3=hyperConvert3D(B3,nr, nc );
    predenoised_blocks3 = ExtractBlocks(B3, bparams);
    Z_block3=zeros( bparams.block_sz(1), bparams.block_sz(2),L, bparams.block_num(1)* bparams.block_num(2));


    for mn=1:max(aa)
        gg=find(aa==mn);
        XES=predenoised_blocks1(:,:,:,gg);
        [a b c d ]=size(XES);

        a1=min(a*b,c*d);
        a2=min(a*b*c,d);
        a3=a;
        c1=sqrt(a1)/(sqrt(a1)+sqrt(a2)+sqrt(a3));
        c2=sqrt(a2)/(sqrt(a1)+sqrt(a2)+sqrt(a3));
        c3=sqrt(a3)/(sqrt(a1)+sqrt(a2)+sqrt(a3));


        D1=C*c1;
        D2=C*c2;
        D3=C*c3;
        nsig2=1;
        XES = reshape(XES,[a*b c*d]);
        V11=WNNM( XES, D1/2/mu, nsig2 );
        V11=reshape(V11,[a b c d]);
        Z_block1(:,:,:,gg)=V11;

        XES=predenoised_blocks2(:,:,:,gg);
        [a b c d ]=size(XES);
        XES = reshape(XES,[a*b*c d]);
        V22=WNNM( XES, D2/2/mu, nsig2 );
        V22=reshape(V22,[a b c d]);
        Z_block2(:,:,:,gg)=V22;


        XES=predenoised_blocks3(:,:,:,gg);
        [a b c d ]=size(XES);
        XES = reshape(XES,[a b*c*d]);
        V33=WNNM( XES, D3/2/mu, nsig2 );
        V33=reshape(V33,[a b c d]);
        Z_block3(:,:,:,gg)=V33;
    end

    V1= JointBlocks(Z_block1, bparams);
    V1=hyperConvert2D(V1);
    V2= JointBlocks(Z_block2, bparams);
    V2=hyperConvert2D(V2);
    V3= JointBlocks(Z_block3, bparams);
    V3=hyperConvert2D(V3);
    % tt1(i)=(norm(2*mu*(V1-ZE),'fro')+norm(2*mu*(V2-ZE),'fro'))/(norm(G1,'fro')+norm(G2,'fro'))
    % if tt1<0.15
    %     break
    % end
    G1=G1+2*mu*(V1-ZE);
    G2=G2+2*mu*(V2-ZE);
    G3=G3+2*mu*(V3-ZE);
    % if i>=2
    % 
    %    indictor(i)= abs(norm(Ztold(:))-norm(Zt(:)))/sqrt(numel(Zt));
    % end
    % Ztold=Zt;
end
Z=Zt;
