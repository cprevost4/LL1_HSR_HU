clear all
close all
clc

load('init_ip_fusion.mat')

%% Load data

SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption
SRI(1,:,:) = []; SRI(:,1,:) = [];
SRI = denoising(SRI);
Pm = spectral_deg(SRI,"LANDSAT");
MSI = tmprod(SRI,Pm,3);

d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);

HSI = awgn(HSI,30,'measured');
MSI = awgn(MSI,30,'measured');

%% BTD-Var

R = 6; L = 8; 
nIter = 20; lamda = 1;

tic;
[A_hat,B_hat,S,C_hat,C_tilde,cost,valid] = BTD_Var(SRI,HSI,MSI,P1,P2,Pm,R,B0,C0,Cbar0,nIter,lamda);
t1 = toc;
Zhat1 = ll1gen({A_hat,B_hat,C_hat},L*ones(1,R));

err1 = compute_metrics(SRI,Zhat1,d1,d2)

%% Alg. 2

options.strategy = 'fusion';
R = 6; L = 8; 
nIter = 20; lambda = 0; innerIter = 5; rho = 1e-3;

tic;
[ZS,ZC,ZCbar,cost] = cnn_btd_regul(SRI,HSI,MSI,B0,C0,Cbar0,P1,P2,L,R,nIter,innerIter,rho,lambda, options);
t2 = toc;

Zhat = reshape(ZS*ZC',size(SRI));
err2 = compute_metrics(SRI,Zhat,d1,d2)

%% STEREO

maxit = 25; t_rank = 50;
MAXIT = 10; lamda = 1;

tic;
[A_hat,B_hat,C_hat,A_tilde,B_tilde,C_tilde] = TenRec(MSI,tens2mat(HSI,[],3),maxit,t_rank,P1,P2);
[A,B,C,cost,valid] = STEREO(SRI,HSI,MSI,P1,P2,Pm,MAXIT,lamda,A_hat,B_hat,C_hat,C_tilde);
t3 = toc;

Zhat = cpdgen({A,B,C});
err3 = compute_metrics(SRI,Zhat,d1,d2)

%% SCOTT

R = [30 30 16];
tic;
Zhat = scott2(HSI, MSI, P1, P2, Pm, R);
t4 = toc;
err4 = compute_metrics(SRI,Zhat,d1,d2)

%% CNN-BTD

R = 6; L = 8;
nIter = 20; rho = 1e-5; innerIter = 5; lamda = 1;

[ZA,ZB,ZC,cost,valid] = ll1_als_cnn1(SRI,HSI,MSI,B0,C0,P1,P2,Pm,L,R,lamda,rho,nIter,innerIter);
t5 = toc;
Zhat5 = reshape(pw_vec(ZB,ZA,R)*ZC',size(SRI));

err5 = compute_metrics(SRI,Zhat5,d1,d2)

%% CNMF

opts.P = 30;

tic;
Zhat = cnmf_adaptor(HSI, MSI, P1, P2, Pm, R, opts);
t6 = toc;
err6 = compute_metrics(SRI,Zhat,d1,d2)

%% GLP-HS

opts.P = 30;

tic;
Zhat = glphs_adaptor(HSI, MSI, P1, P2, Pm, R, opts);
t7 = toc;
err7 = compute_metrics(SRI,Zhat,d1,d2)

%% HySure

opts.P = 30;

tic;
Zhat = hysure_adaptor(HSI, MSI, P1, P2, Pm, R, opts);
t8 = toc;
err8 = compute_metrics(SRI,Zhat,d1,d2)

%% CT-STAR

Rz = [15 15 8]; Rpsi = [3 3 2];

tic;
Zhat10 = hosvdvar(MSI,HSI,Rz,Rpsi,P1,P2,Pm);
t10 = toc;
err10 = compute_metrics(SRI,Zhat10,d1,d2)

%% CB-STAR

Rz = [40 40 4]; Rpsi = [40 40 5];

tic;
Zhat11 = hosvdvaropt(MSI,HSI,P1,P2,Pm,Rz,Rpsi,1);
t11 = toc;
err11 = compute_metrics(SRI,Zhat11,d1,d2)

%% Results

res = ["Algorithm" "R-SNR" "CC" "SAD" "ERGAS" "Time";
        "Alg. 1" err1{:} t1; "Alg. 2" err2{:} t2;
       "STEREO" err3{:} t3; "SCOTT" err4{:} t4;
       "CNN-BTD" err5{:} t5; "CNMF" err6{:} t6;
       "GLP-HS" err7{:} t7; "HySure" err8{:} t8;
       "CT-STAR" err10{:} t10;
       "CB-STAR" err11{:} t11;
       ]
   
   %% Figures
   
   SRI(:,:,40) = SRI(:,:,40)./max(max(SRI(:,:,40)));
   Zhat1(:,:,40) = Zhat1(:,:,40)./max(max(Zhat1(:,:,40)));
   Zhat(:,:,40) = Zhat(:,:,40)./max(max(Zhat(:,:,40)));
   Zhat5(:,:,40) = Zhat5(:,:,40)./max(max(Zhat5(:,:,40)));
   Zhat10(:,:,40) = Zhat10(:,:,40)./max(max(Zhat10(:,:,40)));
   Zhat11(:,:,40) = Zhat11(:,:,40)./max(max(Zhat11(:,:,40)));

figure
subplot(1,6,1); imagesc(SRI(:,:,40));
colorbar; caxis([min(min(SRI(:,:,40))) max(max(SRI(:,:,40)))])
title('Reference'); set(gca,'FontName','Times','FontSize',16)
subplot(1,6,2); imagesc(Zhat1(:,:,40))
colorbar; caxis([min(min(SRI(:,:,40))) max(max(SRI(:,:,40)))])
title('Alg. 1'); set(gca,'FontName','Times','FontSize',16)
subplot(1,6,3); imagesc(Zhat(:,:,40))
colorbar; caxis([min(min(SRI(:,:,40))) max(max(SRI(:,:,40)))])
title('Alg. 2'); set(gca,'FontName','Times','FontSize',16)
subplot(1,6,4); imagesc(Zhat5(:,:,40))
colorbar; caxis([min(min(SRI(:,:,40))) max(max(SRI(:,:,40)))])
title('CNN-BTD'); set(gca,'FontName','Times','FontSize',16)
subplot(1,6,5); imagesc(Zhat10(:,:,40))
colorbar; caxis([min(min(SRI(:,:,40))) max(max(SRI(:,:,40)))])
title('CT-STAR'); set(gca,'FontName','Times','FontSize',16)
subplot(1,6,6); imagesc(Zhat11(:,:,40))
colorbar; caxis([min(min(SRI(:,:,40))) max(max(SRI(:,:,40)))])
title('CB-STAR'); set(gca,'FontName','Times','FontSize',16)
saveas(gcf,'figures/sri_ip.fig')   
