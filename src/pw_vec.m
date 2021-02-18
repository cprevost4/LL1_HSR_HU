function out = pw_vec(X,Y,R)

L = size(X,2)/R;
% Verify if you can actually perform this operation
if L ~= size(Y,2)/R
    fprintf("You cannot perform operation !")
    return
end

out = [];

for r=1:R
    out = [out kr(X(:,(r-1)*L+1:r*L),Y(:,(r-1)*L+1:r*L))*ones(L,1)];
end

end

