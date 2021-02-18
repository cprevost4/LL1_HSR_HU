clear all
close all
clc

load('init_ivanpah_fusion')

%% Load data

load('Playa_preproc_subim3_.mat')
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


%% BTD-Var

R = 4; L = 18; 
nIter = 20; lamda = 1;

tic;
[A_hat,B_hat,S,C_hat,C_tilde,cost,valid] = BTD_Var(SRI,HSI,MSI,P1,P2,Pm,R,B0,C0,Cbar0,nIter,lamda);
t1 = toc;
Zhat1 = ll1gen({A_hat,B_hat,C_hat},L*ones(1,R));
Psihat1 = ll1gen({A_hat,B_hat,C_tilde-Pm*C_hat},L*ones(1,R));

err1 = compute_metrics(SRI,Zhat1,d1,d2);
err1_psi = compute_metrics(Psi,Psihat1,d1,d2);

%% CNN-BTD-Var

R = 4; L = 18; 
nIter = 20; lambda = 0; innerIter = 5; rho = 1e-3;

tic;
[ZS,ZC,ZCbar,cost] = cnn_btd_regul(SRI,HSI,MSI,B0,C0,Cbar0,P1,P2,L,R,nIter,innerIter,rho,lambda);
t2 = toc;

Zhat2 = reshape(ZS*ZC',size(SRI));
Psihat2 = reshape(ZS*(ZCbar-Pm*ZC)',size(MSI));

err2 = compute_metrics(SRI,Zhat2,d1,d2);
err2_psi = compute_metrics(Psi,Psihat2,d1,d2);

%% STEREO

maxit = 25; t_rank = 10;
MAXIT = 10; lamda = 1;

tic;
[A_hat,B_hat,C_hat,A_tilde,B_tilde,C_tilde] = TenRec(MSI,tens2mat(HSI,[],3),maxit,t_rank,P1,P2);
[A,B,C,cost,valid] = STEREO(SRI,HSI,MSI,P1,P2,Pm,MAXIT,lamda,A_hat,B_hat,C_hat,C_tilde);
t3 = toc;

Zhat3 = cpdgen({A,B,C});
err3 = compute_metrics(SRI,Zhat3,d1,d2);

%% SCOTT

R = [30 30 10];
tic;
Zhat4 = scott2(HSI, MSI, P1, P2, Pm, R);
t4 = toc;
err4 = compute_metrics(SRI,Zhat4,d1,d2);

%% CNN-BTD

R = 4; L = 18; 
nInit = 20; obj = 10^50;

nIter = 20; rho = 1e-3; innerIter = 5; lamda = 1;
[ZA,ZB,ZC,cost,valid] = ll1_als_cnn1(SRI,HSI,MSI,B0,C0,P1,P2,Pm,L,R,lamda,rho,nIter,innerIter);
t5 = toc;
Zhat5 = reshape(pw_vec(ZB,ZA,R)*ZC',size(SRI));

err5 = compute_metrics(SRI,Zhat5,d1,d2);

%% CNMF

opts.P = 30;

tic;
Zhat6 = cnmf_adaptor(HSI, MSI, P1, P2, Pm, R, opts);
t6 = toc;
err6 = compute_metrics(SRI,Zhat6,d1,d2);

%% GLP-HS

opts.P = 30;

tic;
Zhat7 = glphs_adaptor(HSI, MSI, P1, P2, Pm, R, opts);
t7 = toc;
err7 = compute_metrics(SRI,Zhat7,d1,d2);

%% HySure

opts.P = 30;

tic;
Zhat8 = hysure_adaptor(HSI, MSI, P1, P2, Pm, R, opts);
t8 = toc;
err8 = compute_metrics(SRI,Zhat8,d1,d2);

%% CT-STAR

Rz = [15 15 8]; Rpsi = [3 3 2];

tic;
[Zhat10, Psihat10] = hosvdvar(MSI,HSI,Rz,Rpsi,P1,P2,Pm);
t10 = toc;

err10 = compute_metrics(SRI,Zhat10,d1,d2);
err10_psi = compute_metrics(Psi,Psihat10,d1,d2);

%% CB-STAR

Rz = [40 40 4]; Rpsi = [40 40 5];

tic;
[Zhat11, Psihat11] = hosvdvaropt(MSI,HSI,P1,P2,Pm,Rz,Rpsi,1);
t11 = toc;

err11 = compute_metrics(SRI,Zhat11,d1,d2);
err11_psi = compute_metrics(Psi,Psihat11,d1,d2);

%% Results

res_Z = ["Algorithm" "R-SNR" "CC" "SAD" "ERGAS" "Time";
        "Alg. 1" err1{:} t1; "Alg. 2" err2{:} t2;
       "STEREO" err3{:} t3; "SCOTT" err4{:} t4;
       "CNN-BTD" err5{:} t5; "CNMF" err6{:} t6;
       "GLP-HS" err7{:} t7; "HySure" err8{:} t8;
       "CT-STAR" err10{:} t10;
       "CB-STAR" err11{:} t11;
       ]
   
res_Psi = ["Algorithm" "R-SNR" "CC" "SAD" "ERGAS" "Time";
        "Alg. 1" err1_psi{:} t1; "Alg. 2" err2_psi{:} t2;
       "CT-STAR" err10_psi{:} t10;
       "CB-STAR" err11_psi{:} t11;
       ]  
   
%% Figures

figure(1)
subplot(2,2,1); imagesc(SRI(:,:,40)); colorbar; caxis([min(min(SRI(:,:,40))) max(max(SRI(:,:,40)))])
title('Reference'); set(gca,'FontName','Times','FontSize',16)
subplot(2,2,2); imagesc(Zhat1(:,:,40)); colorbar;  caxis([min(min(SRI(:,:,40))) max(max(SRI(:,:,40)))])
title('Alg. 1'); set(gca,'FontName','Times','FontSize',16)
subplot(2,2,3); imagesc(Zhat2(:,:,40)); colorbar;  caxis([min(min(SRI(:,:,40))) max(max(SRI(:,:,40)))])
title('Alg. 2'); set(gca,'FontName','Times','FontSize',16)
subplot(2,2,4); imagesc(Zhat11(:,:,40)); colorbar;  caxis([min(min(SRI(:,:,40))) max(max(SRI(:,:,40)))])
title('CB-STAR'); set(gca,'FontName','Times','FontSize',16)
saveas(gcf,'figures/sri_ivanpah.fig')   

figure(2)
subplot(2,2,1); imagesc(Psi(:,:,6)); colorbar; caxis([min(min(Psi(:,:,6))) max(max(Psi(:,:,6)))])
title('Reference'); set(gca,'FontName','Times','FontSize',16)
subplot(2,2,2); imagesc(Psihat1(:,:,6)); colorbar;  caxis([min(min(Psi(:,:,6))) max(max(Psi(:,:,6)))])
title('Alg. 1'); set(gca,'FontName','Times','FontSize',16)
subplot(2,2,3); imagesc(Psihat2(:,:,6)); colorbar;  caxis([min(min(Psi(:,:,6))) max(max(Psi(:,:,6)))])
title('Alg. 2'); set(gca,'FontName','Times','FontSize',16)
subplot(2,2,4); imagesc(Psihat11(:,:,6)); colorbar;  caxis([min(min(Psi(:,:,6))) max(max(Psi(:,:,6)))])
title('CB-STAR'); set(gca,'FontName','Times','FontSize',16)
saveas(gcf,'figures/psi_ivanpah.fig')   