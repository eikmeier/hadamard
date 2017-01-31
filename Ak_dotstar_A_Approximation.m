function[C] = Ak_dotstar_A_Approximation(A)
% COMPUTE APPROXIMATION takes a graph and computes a rank 1 factorization of 
% the adjacency matrix of the graph A^(k-1) where x and then computes
% a Hadamard product.
%
% 


[m,n] = size(A);
k = floor(log(n));
if mod(k,2) ~= 0
    k = k+1;
end

    
%%Compute Aproximation

        %randomized range finder
l = 20;
Om = randn(n,l);
Y = sparse(m,l);
for i = 1:l
    Y(:,i) = A*Om(:,i);
end 
Y = full(Y);
[Q,~] = qr(Y,0);
Q = sparse(Q);

        %direct svd
     %do a rank q approximation
q = 10;
D = Q'*A;
[Utemp,S,V] = svds(D,q);
U = zeros(m,q);
for i = 1:q
    U(:,i) = Q*Utemp(:,i);
end

    %form approximation C
V = V(:,1:q);
S = S(1:q,1:q);
S = S.^(k-1);
Vt = S*V';
Ut = U';

%Composes C ~ A^(k-1)

[ei,ej] = find(A); % give us all the non-zeros
cvals = zeros(numel(ei),1);
for nzi=1:numel(ei) % for each non-zero/edge
    i = ei(nzi);
    j = ej(nzi);
    cvals(nzi) = Ut(:,i)'*(Vt(:,j));
end

C = sparse(ei,ej,cvals,size(A,1),size(A,2));

%Composes C = C.*A
for i = 1:n
    C(:,i) = C(:,i).*A(:,i);
end

% Save C
%save('comyoutube-ungraph-approxrank1.mat','C');


