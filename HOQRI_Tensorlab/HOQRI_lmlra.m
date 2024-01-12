clear
close all
% clc

i=100;
J = [i,i,i];
k = 10;
K = [k,k,k];

sparsity = 10*i;
spt = sptenrand(J, sparsity); 
normX2 = norm(spt)^2;


U0 = orth(rand(J(1),K(1)));
V0 = orth(rand(J(2),K(2)));
W0 = orth(rand(J(3),K(3)));


init = {U0,V0,W0};
X=double(spt);
options.MaxIter=2

fprintf('hooi\n');
fprintf('itr    time     error\n');
[U_hooi,S_hooi,output_hooi] = lmlra_hooi_time(X,init,options);

fprintf('nls\n');
fprintf('itr    time     error\n');
[U_nls,S_nls,output_nls] = lmlra_nls(X,init,rand(K),options);


fprintf('minf\n');
fprintf('itr    time     error\n');
[U_minf,S_minf,output_minf] = lmlra_minf(X,init,rand(K),options);



T.val = spt.vals;
T.sub = spt.subs;
T.size = J;
% T.ind = spt.subs(:,1)+ J(1)*spt.subs(:,2) + J(1)*J(2)*(spt.subs(:,3)-1);
T.sparse  = true;
X=fmt(T);

U = U0;
V = V0;
W = W0;

normG2 = norm( mtkronprod( X, {U,V,W},0) )^2;
fit = normG2;
tim = [];
tic
fprintf('hoqri\n');
fprintf('itr    time     error\n');
for itr = 1:10

    A = 0;
    v23 = full(mtkronprod(X,{U,V,W},1));
    for k=1:K(2)*K(3)
        A = A + v23(:, k) * (v23(:, k)'*U);
    end
    [U,~] = qr(A,0);

    B = 0;
    v13 = full(mtkronprod(X,{U,V,W},2));
    for k=1:K(1)*K(3)
        B = B + v13(:, k) * (v13(:, k)'*V);
    end
    [V,~] = qr(B,0);    

    C = 0;
    v12 = full(mtkronprod(X,{U,V,W},3));
    for k=1:K(1)*K(2)
        C = C + v12(:, k) * (v12(:, k)'*W);
    end
    [W,~] = qr(C,0);

    newnormG2 = norm( mtkronprod( X, {U,V,W}, 0) )^2;
    if abs(   normG2 - newnormG2 )<=1e-5
	break
    end
    normG2 = newnormG2;
    fit = [ fit, normG2];
    tim = [ tim, toc ];
    fprintf('%d %f %f\n',itr, toc, sqrt(normX2 - normG2));
    
end


% disp('	nls		minf		hooi		hoqri');
% disp([output_nls.time(1:10)',output_minf.time(1:10)',output_hooi.time(1:10)',tim(1:10)'])
