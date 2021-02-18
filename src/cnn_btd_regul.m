function [ZS,ZC,ZCbar,cost] = cnn_btd_regul(Z,HSI,MSI,B_hat,C_hat,Cbar_hat,P1,P2,L,R,nIter,innerIter,rho,lambda,options)

if ~exist('options','var')
    options = struct();
end
if ~isfield(options,'strategy') || isempty(options.strategy)
    options.strategy = true;
end

ZS = zeros(size(Z,1)*size(Z,2),R); US = ZS;
ZC = zeros(size(C_hat)); UC = ZC;
ZCbar = zeros(size(Cbar_hat)); UCbar = ZCbar;

cost(1) = Inf; diff_cost(1) = Inf; i = 1;

while i<nIter %&& cost(i)>eps && diff_cost(i)>eps
    
    i = i+1;
    
     % Update A
    temp1h=pw_kron(C_hat,P2*B_hat,R);
    temp1m=pw_kron(Cbar_hat,B_hat,R);
    H1 = (P1'*P1);
    H2 = pw_kron(C_hat,P2*B_hat,R)'*pw_kron(C_hat,P2*B_hat,R);
    H4 = pw_kron(Cbar_hat,B_hat,R)'*pw_kron(Cbar_hat,B_hat,R);
    inv_H2=pinv(H2);
    Bs = H4*inv_H2;
    H5 = P1'*tens2mat(HSI,[],1)'*temp1h + tens2mat(MSI,[],1)'*temp1m;
    A_hat = sylvester(full(H1),Bs,H5*inv_H2);

    % Update B
    temp1h=pw_kron(C_hat,P1*A_hat,R);
    temp1m=pw_kron(Cbar_hat,A_hat,R);
    H1 = (P2'*P2);
    H2 = pw_kron(C_hat,P1*A_hat,R)'*pw_kron(C_hat,P1*A_hat,R);
    H4 = pw_kron(Cbar_hat,A_hat,R)'*pw_kron(Cbar_hat,A_hat,R);
    inv_H2=pinv(H2);
    Bs = H4*inv_H2;
    H5 = P2'*tens2mat(HSI,[],2)'*temp1h + tens2mat(MSI,[],2)'*temp1m;
    B_hat = sylvester(full(H1),Bs,H5*inv_H2);

     % Update S
     I = size(A_hat,1); J = size(B_hat,1);
     rho = sum(diag(eye(I)))/(I);
     %rho = 1e3;
     Hhs = eye(J) + diag(-ones(J-1,1),-1); Hhs(1,end) = -1;
     Hvs = eye(I) + diag(-ones(I-1,1),1); Hvs(end,1) = -1;
     n=1;
    while n<innerIter
        S_hat = [];
        n = n+1;
         for r=1:R
             tmp = A_hat(:,(r-1)*L+1:r*L)*B_hat(:,(r-1)*L+1:r*L)';
             res = sylvester(lambda*Hvs'*Hvs,(1+rho)*eye(J)+lambda*Hhs*Hhs',tmp+rho*reshape(ZS(:,r)+US(:,r),I,J));
             S_hat = [S_hat res(:)];
         end
         ZS = max(0,S_hat-US);
         US = US + (ZS - S_hat);
    end
    
    % Update C
    temp = kron(P2,P1)*S_hat;
    rho = sum(diag(temp'*temp))/size(Z,3);
    %rho = 1e-3;
    H2 = temp'*temp + rho*eye(R);
    H5_tilde = tens2mat(HSI,[],3)'*temp;
    n=1;
    while n<innerIter
        n = n+1;
        H5 = H5_tilde + rho*(ZC+UC);
        C_hat = (pinv(H2)*H5')';
        if options.strategy
            for r=1:R
                C_hat(:,r) = C_hat(:,r)/norm(C_hat(:,r));
            end
        end
        ZC = max(0,C_hat-UC);
        UC = UC + (ZC - C_hat);
    end
    
    % Update Cbar
    temp = S_hat;
    rho = sum(diag(temp'*temp))/size(MSI,3);
    H2 = temp'*temp + rho*eye(R);
    H5_tilde = tens2mat(MSI,[],3)'*temp;
    n=1;
    while n<innerIter
        n = n+1;
        H5 = H5_tilde + rho*(ZCbar+UCbar);
        Cbar_hat = (pinv(H2)*H5')';
        if options.strategy
            for r=1:R
                Cbar_hat(:,r) = Cbar_hat(:,r)/norm(Cbar_hat(:,r));
            end
        end
        ZCbar = max(0,Cbar_hat-UCbar);
        UCbar = UCbar + (ZCbar - Cbar_hat);
    end

    
    cost(i) = frob(HSI - generate_ll1(P1*A_hat,P2*B_hat,C_hat,L,R),'squared') + ...
        frob(MSI - generate_ll1(A_hat,B_hat,Cbar_hat,L,R),'squared');
    diff_cost(i) = cost(i-1)-cost(i);

end

end