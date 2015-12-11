function out = performance_test()

thesum = 100;

NN = [10, 50, 100, 200, 300, 400, 500, 750, 1000, 2000, 5000, 10000];
t_Ainv = zeros(length(NN), 1);
t_AB   = zeros(length(NN), 1);
t_ABC  = zeros(length(NN), 1);
t_Asum = zeros(length(NN), 1);

figure(1);

mytimer = timer('Period',0.01);

for i = 1:length(NN)
    N = NN(i);
    A = rand(N,N);
    B = rand(N,N);
    C = rand(N,N);
    Ainv = zeros(N,N);
    AB   = zeros(N,N);
    ABC  = zeros(N,N);
    Asum = zeros(N,N);
    
    fprintf(1, 'N=%d:\n',N);
    
    fprintf(1, '   matrix inversion...\n');
    tic;
    Ainv = inv(A);
    t_Ainv(i) = toc;
    
    fprintf(1, '   single matrix multiplication...\n');
    tic;
    AB = A*B;
    t_AB(i) = toc;
    
    fprintf(1, '   double matrix multiplication...\n');
    tic;
    ABC = A*B*C;
    t_ABC(i) = toc;
    
    fprintf(1, '   sum of %d matrices...\n',thesum);
    tic;
    for j=1:thesum
        Asum = Asum + A;
    end
    t_Asum(i) = toc;
    
    loglog(NN,t_Ainv,'blue',NN,t_AB,'red',NN,t_ABC,'yellow',NN,t_Asum,'black');
    legend('inversion','multiplication','double mult','addition');
    if N>999
        fprintf(1, '   press any key to continue...');
        pause;
    end
end
