function [A_0,B_0,C_0,Cbar_0] = BTDRec2(HSI,MSI,P1,P2,R,L)

% Initialization function for CNN_BTD

options.MaxIter = 25; options.Initialization = @ll1_rnd;
%options.Initialization = @ll1_gevd;
obj = 10^50;

%Ctilde_0  = []; 
%A_0 = []; B_0 = [];

for i = 1:1
    U = ll1(MSI,L*ones(1,R),options,'OutputFormat','CPD');
    for r=1:R
        Ai = U{1};
        Bi = U{2};
        Ctildei = U{3};
    end
    err = frob(MSI - ll1gen({Ai,Bi,Ctildei},L*ones(1,R)),'squared');
    if err<obj
        obj = err;
        A_0 = Ai; B_0 = Bi;
    end
end


A_0 = abs(A_0); B_0 = abs(B_0);
temp = pw_vec2(P1*A_0,P2*B_0,R);
% temp = pw_vec(P1*A_0,P2*B_0,R);
C_0=(temp\tens2mat(HSI,[],3))';
temp = pw_vec2(A_0,B_0,R);
% temp = pw_vec(A_0,B_0,R);
Cbar_0=(temp\tens2mat(MSI,[],3))';

end

