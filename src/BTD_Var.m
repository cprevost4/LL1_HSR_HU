function [A,B,S,C,C_tilde,cost,valid] = BTD_Var(SRI,HSI,MSI,P1,P2,Pm,R,B,C,C_tilde,Niter,lambda)

m=1; 
cost(1) = inf; diff_cost(1) = inf; tol = 1e-8;
valid(1) = inf; diff_valid(1) = inf;

while m<Niter %&& cost(m)>tol && diff_cost(m)>tol %&& valid(m)>tol && diff_valid(m)>tol
    
    m=m+1;
    
% Update A
    H1 = P1'*P1; H3 = lambda*eye(size(MSI,1));
    H2 = pw_kron(C,P2*B,R)'*pw_kron(C,P2*B,R);
    H4 = pw_kron(C_tilde,B,R)'*pw_kron(C_tilde,B,R);
    H5 = P1'*tens2mat(HSI,[],1)'*pw_kron(C,P2*B,R) + lambda*tens2mat(MSI,[],1)'*pw_kron(C_tilde,B,R);
    A = real(bartelsStewart(H1,H2,H3,H4,H5));
    
% Update B
    H1 = P2'*P2; H3 = lambda*eye(size(MSI,2));
    H2 = pw_kron(C,P1*A,R)'*pw_kron(C,P1*A,R);
    H4 = pw_kron(C_tilde,A,R)'*pw_kron(C_tilde,A,R);
    H5 = P2'*tens2mat(HSI,[],2)'*pw_kron(C,P1*A,R) + lambda*tens2mat(MSI,[],2)'*pw_kron(C_tilde,A,R);
    B = real(bartelsStewart(H1,H2,H3,H4,H5));
    %B = qr(B,0);
     
% Update S
    S = pw_vec(B,A,R);
    S_tilde = kron(P2,P1)*S;
%     for r=1:R
%         S(:,r) = S(:,r)/norm(S(:,r));
%         S_tilde(:,r) = S_tilde(:,r)/norm(S_tilde(:,r));
%     end

% Update C_tilde
    C_tilde = (pinv(S)*tens2mat(MSI,[],3))';
%     for r=1:R
%         C_tilde(:,r) = C_tilde(:,r)/norm(C_tilde(:,r));
%     end
    
% Update C
    C = (pinv(S_tilde)*tens2mat(HSI,[],3))';
%     for r=1:R
%         C(:,r) = C(:,r)/norm(C(:,r));
%     end
    

    
    cost(m) = frob(tens2mat(HSI,[],3) - kron(P2,P1)*pw_vec(A,B,R)*C','squared') + ...
    lambda*frob(tens2mat(MSI,[],3) - pw_vec(A,B,R)*C_tilde','squared');
    diff_cost(m) = cost(m-1) - cost(m);
    valid(m) = frob(tens2mat(SRI,[],3) - pw_vec(A,B,R)*C','squared');
    diff_valid(m) = valid(m-1) - valid(m);
    
    
    
    
    
end




end

