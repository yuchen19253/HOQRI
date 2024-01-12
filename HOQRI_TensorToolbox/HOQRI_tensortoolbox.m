clear;
% clc;

i = 100;
J = [i,i,i];
k = 10;
K = [k,k,k];

spt = sptenrand(J,100*i); 
normX2 = norm(spt)^2; 


[~,init,obj,tim] = tucker_als( spt, K, 'tol', 1e-8, 'maxiters',10, 'printitn',0);


U = orth(init{1});
V = orth(init{2});
W = orth(init{3});

normG = norm( ttm( spt, {U',V',W'}, [1,2,3] ) )^2;
fit = normG;
tim = [];
tic
fprintf('hoqri\n');
for itr = 1:10
    A=0;
    for j=1:K(2)
        for k=1:K(3)
            v = double(ttv(spt,{V(:,j),W(:,k)},[2,3]));
            A = A + v * (v' * U);
        end
    end
    [U,~] = qr(A,0);
    
    B=0;
    for i=1:K(1)
        for k=1:K(3)
            v = double(ttv(spt,{U(:,i),W(:,k)},[1,3]));
            B = B + v * (v' * V);
        end
    end
    [V,~] = qr(B,0);
    
    C=0;
    for i=1:K(1)
        for j=1:K(2)
            v = double(ttv(spt,{U(:,i),V(:,j)},[1,2]));
            C = C + v * (v' * W);
        end
    end
    [W,~] = qr(C,0);

    G = ttm( spt, {U',V',W'}, [1,2,3] );
    newnormG = norm( G )^2;
    loss=abs(   normG - newnormG );
    if loss<=1e-8
	break
    end
    normG = newnormG;
    fit = [ fit, normG];
    tim = [ tim, toc ];
    fprintf('itr = %d time = %f loss = %f normG = %f\n',itr, toc, loss, sqrt(normG) );
    
end

