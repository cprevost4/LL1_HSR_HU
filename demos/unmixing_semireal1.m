clear all
close all
clc

load('init_laketahoe_unmixing')

%% Load data

load('Tahoe_preproc_t1.mat')
clear MSI; load('MSI_tahoe_better_registered.mat') % load a better registered MSI
Pm = SRF; clear SRF %Spectral response function
Z = HSI; clear HSI %Get SRI

d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(Z, q, d1, d2);
HSI = tmprod(tmprod(Z,P1,1),P2,2);

load('endmembers_tahoe2.mat')
C = M0; clear M0; clear M_VCA
Cbar_ref = vca(tens2mat(MSI,[],3)',3);
load('data_Tahoe_abundances.mat')

%% Initialization

R = 3; L = 20; 
nIter = 20; lambda = 0; innerIter = 5; rho = 1e-3;

tic;
[ZS,ZC,ZCbar,cost] = cnn_btd_regul(Z,HSI,MSI,B0,C0,Cbar0,P1,P2,L,R,nIter,innerIter,rho,lambda);
t1 = toc;

Zhat = reshape(ZS*ZC',size(Z));

[ZC, ind] = sort_endmembers_to_ref(C,ZC);
ZS = ZS(:,ind);

%% CNMF

M_em = R;

tic;
[Out,W_hyper,H_hyper,W_multi,H_multi] = CNMF_fusion(HSI,MSI,M_em);
t2 = toc;
W_hyper = W_hyper(1:end-1,:); H_multi = H_multi';

[W_hyper, ind] = sort_endmembers_to_ref(C,W_hyper);
H_multi = H_multi(:,ind);

%% Mu-ACC

U0 = rand(size(Z,3),R); V0 = rand(R,size(Z,1)*size(Z,2));
maxiter = 1e6; timelimit = 5;
tic;
[Uma,Vma] = MUacc(tens2mat(Zhat,[],3)',U0,V0,2,0.1,maxiter,timelimit);
t4 = toc;

[Uma, ind] = sort_endmembers_to_ref(C,Uma);
Vma = Vma'; Vma = Vma(:,ind);

%% BMDR-ADMM

mu = 50;
tic;
[V,U] = adaptor_bmdr(Zhat,R,mu);
t5 = toc;
[U, ind] = sort_endmembers_to_ref(C,U);
V = V'; V = V(:,ind);

%% Figures
 
figure(1)
for r=1:R
    subplot(R,5,5*(r-1)+1); imagesc(A_VCA(:,:,r))
    set(gca,'FontName','Times','FontSize',16)
    subplot(R,5,5*(r-1)+2); imagesc(reshape(ZS(:,r),size(MSI,1),size(MSI,2))) 
    set(gca,'FontName','Times','FontSize',16)
    subplot(R,5,5*(r-1)+3); imagesc(reshape(H_multi(:,r),size(MSI,1),size(MSI,2))) 
    set(gca,'FontName','Times','FontSize',16)
    subplot(R,5,5*(r-1)+4); imagesc(reshape(Vma(:,r),size(MSI,1),size(MSI,2))) 
    set(gca,'FontName','Times','FontSize',16)
    subplot(R,5,5*(r-1)+5); imagesc(reshape(V(:,r),size(MSI,1),size(MSI,2))) 
    set(gca,'FontName','Times','FontSize',16)
end
%saveas(gcf,'figures/abundance_tahoe.fig')   

for r=1:R
    ZC(:,r) = ZC(:,r)/norm(ZC(:,r));
    C(:,r) = C(:,r)/norm(C(:,r));
    W_hyper(:,r) = W_hyper(:,r)/norm(W_hyper(:,r));
    Uma(:,r) = Uma(:,r)/norm(Uma(:,r));
    U(:,r) = U(:,r)/norm(U(:,r));
end

figure(2)
for r=1:R
    subplot(R,5,5*(r-1)+1); plot(C(:,r),'k--','LineWidth',1);  
    set(gca,'FontName','Times','FontSize',16)
     subplot(R,5,5*(r-1)+2); plot(ZC(:,r),'r--','LineWidth',1);
     set(gca,'FontName','Times','FontSize',16)
    subplot(R,5,5*(r-1)+3); plot(W_hyper(:,r),'g--','LineWidth',1); 
    set(gca,'FontName','Times','FontSize',16)
    subplot(R,5,5*(r-1)+4); plot(Uma(:,r),'b--','LineWidth',1); 
    set(gca,'FontName','Times','FontSize',16)
    subplot(R,5,5*(r-1)+5); plot(U(:,r),'m--','LineWidth',1); 
    set(gca,'FontName','Times','FontSize',16)
end
%saveas(gcf,'figures/spectra_tahoe.fig') 


S  = tens2mat(Z,[],3)'\C;

sad_tahoe = [sad(C,ZC) sad(C,W_hyper) sad(C,Uma) sad(C,U)]
rmse_tahoe = [rmse(S,ZS) rmse(C,H_multi) rmse(C,Vma) rmse(C,V)]
