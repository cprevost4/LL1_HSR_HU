clear all
close all
clc

%% Fusion

load('Tahoe_preproc_t1.mat')
load('init_laketahoe_fusion.mat')
clear MSI; load('MSI_tahoe_better_registered.mat') 
Pm = SRF; clear SRF %Spectral response function

MSI(MSI<1e-3) = 1e-3; HSI(HSI<1e-3) = 1e-3; MSI(MSI>1)=1; HSI(HSI>1)=1;
MSim = MSI; HSim = HSI;
MSp = zeros(size(MSim));
for i=1:size(MSim,3)
    x = MSim(:,:,i);
    xmax  = quantile(x(:),0.999);
    % Normalize to 1
    x = x/xmax;
    MSp(:,:,i) = x;
end
MSI = MSp;
HSp = zeros(size(HSim));
for i=1:size(HSim,3)
    x = HSim(:,:,i);
    xmax  = quantile(x(:),0.999);
    %normalize to 1
    x = x/xmax;
    HSp(:,:,i) = x;
end
HSI = HSp;

SRI = HSI; clear HSI %Get SRI
SRI = denoising(SRI);

d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);

Psi = MSI - tmprod(SRI,Pm,3);

HSI = awgn(HSI,30,'measured');
MSI = awgn(MSI,30,'measured');

% BTD-Var
R = 3; L = 20; 
nIter = 20; lamda = 1;
tic;
[A_hat,B_hat,S,C_hat,C_tilde,cost,valid] = BTD_Var(SRI,HSI,MSI,P1,P2,Pm,R,B0,C0,Cbar0,nIter,lamda);
t1 = toc;
Zhat1 = ll1gen({A_hat,B_hat,C_hat},L*ones(1,R));
Psihat1 = ll1gen({A_hat,B_hat,C_tilde-Pm*C_hat},L*ones(1,R));
err1 = compute_metrics(SRI,Zhat1,d1,d2);
err1_psi = compute_metrics(Psi,Psihat1,d1,d2);

% CNN-BTD-Var
options.strategy = false;
nIter = 20; innerIter = 5; rho = 1e-3; lambda = 0; % No regularization
tic;
[ZS,ZC,ZCbar,cost] = cnn_btd_regul(SRI,HSI,MSI,B0,C0,Cbar0,P1,P2,L,R,nIter,innerIter,rho,lambda,options);
t2 = toc;
Zhat2 = reshape(ZS*ZC',size(SRI));
Psihat2 = reshape(ZS*(ZCbar-Pm*ZC)',size(MSI));
err2 = compute_metrics(SRI,Zhat2,d1,d2);
err2_psi = compute_metrics(Psi,Psihat2,d1,d2);

% Results
res_Z = ["Algorithm" "R-SNR" "CC" "SAD" "ERGAS" "Time";
        "Alg. 1" err1{:} t1; "Alg. 2" err2{:} t2;]
res_Psi = ["Algorithm" "R-SNR" "CC" "SAD" "ERGAS" "Time";
        "Alg. 1" err1_psi{:} t1; "Alg. 2" err2_psi{:} t2;]
figure(1)
subplot(2,2,1); imagesc(SRI(:,:,40)); colorbar; caxis([min(min(SRI(:,:,40))) max(max(SRI(:,:,40)))])
title('Reference'); set(gca,'FontName','Times','FontSize',16)
subplot(2,2,2); imagesc(Zhat1(:,:,40)); colorbar;  caxis([min(min(SRI(:,:,40))) max(max(SRI(:,:,40)))])
title('Alg. 1'); set(gca,'FontName','Times','FontSize',16)
subplot(2,2,3); imagesc(Zhat2(:,:,40)); colorbar;  caxis([min(min(SRI(:,:,40))) max(max(SRI(:,:,40)))])
title('Alg. 2'); set(gca,'FontName','Times','FontSize',16)
figure(2)
subplot(2,2,1); imagesc(Psi(:,:,6)); colorbar; caxis([min(min(Psi(:,:,6))) max(max(Psi(:,:,6)))])
title('Reference'); set(gca,'FontName','Times','FontSize',16)
subplot(2,2,2); imagesc(Psihat1(:,:,6)); colorbar;  caxis([min(min(Psi(:,:,6))) max(max(Psi(:,:,6)))])
title('Alg. 1'); set(gca,'FontName','Times','FontSize',16)
subplot(2,2,3); imagesc(Psihat2(:,:,6)); colorbar;  caxis([min(min(Psi(:,:,6))) max(max(Psi(:,:,6)))])
title('Alg. 2'); set(gca,'FontName','Times','FontSize',16)

%% Unmixing

clear all
load('Tahoe_preproc_t1.mat')
load('init_laketahoe_unmixing')
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

R = 3; L = 20; 
nIter = 20; lambda = 0; innerIter = 5; rho = 1e-3;
tic;
[ZS,ZC,ZCbar,cost] = cnn_btd_regul(Z,HSI,MSI,B0,C0,Cbar0,P1,P2,L,R,nIter,innerIter,rho,lambda);
t1 = toc;

Zhat = reshape(ZS*ZC',size(Z));
[ZC, ind] = sort_endmembers_to_ref(C,ZC);
ZS = ZS(:,ind);

 
figure(1)
for r=1:R
    subplot(R,2,2*(r-1)+1); imagesc(A_VCA(:,:,r))
    set(gca,'FontName','Times','FontSize',16)
    subplot(R,2,2*(r-1)+2); imagesc(reshape(ZS(:,r),size(MSI,1),size(MSI,2))) 
    set(gca,'FontName','Times','FontSize',16)
end

for r=1:R
    ZC(:,r) = ZC(:,r)/norm(ZC(:,r));
    C(:,r) = C(:,r)/norm(C(:,r));
end

figure(2)
for r=1:R
    subplot(R,2,2*(r-1)+1); plot(C(:,r),'k--','LineWidth',1);  
    set(gca,'FontName','Times','FontSize',16)
     subplot(R,2,2*(r-1)+2); plot(ZC(:,r),'r--','LineWidth',1);
     set(gca,'FontName','Times','FontSize',16)
end

