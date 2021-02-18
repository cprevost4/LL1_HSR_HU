clc
clear all
close all

%% Load data

load('Playa_preproc_subim3_.mat')
Pm = SRF; clear SRF %Spectral response function
SRI = HSI; clear HSI %Get SRI 
Z = SRI;

d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);

load('endmembers_playa1.mat')
Cref = M0; C = Cref; clear M0;
Sref = tens2mat(SRI,[],3)/Cref'; Sref = max(Sref,0); S = Sref;

%% Algo

R = 4; L = 18; 
load('init_ivanpah_unmixing')

nIter = 20; lambda = 0; innerIter = 5; rho = 1e-3;

%You need to uncomment the rho lines for the algorithm
tic;
[ZS,ZC,ZCbar,cost] = cnn_btd_regul(Z,HSI,MSI,B0,C0,Cbar0,P1,P2,L,R,nIter,innerIter,rho,lambda);
t1 = toc;

Zhat = reshape(ZS*ZC',size(SRI));

%% CNMF

M_em = R;

tic;
[Out,W_hyper,H_hyper,W_multi,H_multi] = CNMF_fusion(HSI,MSI,M_em);
t2 = toc;
W_hyper = W_hyper(1:end-1,:); H_multi = H_multi';

%% Mu-Acc

U0 = rand(size(Z,3),R); V0 = rand(R,size(Z,1)*size(Z,2));
maxiter = 1e6; timelimit = 5;
tic;
[Uma,Vma] = MUacc(tens2mat(Zhat,[],3)',U0,V0,2,0.1,maxiter,timelimit);
t4 = toc;

Vma = Vma';

%% BMDR-ADMM

mu = 1;
tic;
[V,U] = adaptor_bmdr(Zhat,R,mu);
t5 = toc;
 V = V';

%% Maps 

H_multi = sort_endmembers_to_ref(Sref,H_multi);
Vma = sort_endmembers_to_ref(Sref,Vma);
V = sort_endmembers_to_ref(Sref,V);
ZS = sort_endmembers_to_ref(Sref,ZS);

for r=1:R
    Sref(:,r) = Sref(:,r)/norm(Sref(:,r));
    ZS(:,r) = ZS(:,r)/norm(ZS(:,r));
    H_multi(:,r) = H_multi(:,r)/norm(H_multi(:,r));
    Vma(:,r) = Vma(:,r)/norm(Vma(:,r));
    V(:,r) = V(:,r)/norm(V(:,r));
end

figure(1)
for r=1:R
    subplot(R,5,5*(r-1)+1); imagesc(reshape(Sref(:,r),size(MSI,1),size(MSI,2)))
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
saveas(gcf,'figures/abundance_ivanpah.fig') 

for r=1:R
    Cref(:,r) = Cref(:,r)/norm(Cref(:,r));
    ZC(:,r) = ZC(:,r)/norm(ZC(:,r));
    W_hyper(:,r) = W_hyper(:,r)/norm(W_hyper(:,r));
    Uma(:,r) = Uma(:,r)/norm(Uma(:,r));
    U(:,r) = U(:,r)/norm(U(:,r));
end

figure(2)
for r=1:R
    subplot(R,5,5*(r-1)+1); plot(Cref(:,r),'k--','LineWidth',1);  
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
saveas(gcf,'figures/spectra_ivanpah.fig') 
