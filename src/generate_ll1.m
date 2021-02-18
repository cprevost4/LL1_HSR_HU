function tens = generate_ll1(A,B,C,L,R)

tens = zeros(size(A,1),size(B,1),size(C,1));

for r=1:R
    tens = tens + outprod(A(:,(r-1)*L+1:r*L)*B(:,(r-1)*L+1:r*L)',C(:,r));
end

end

