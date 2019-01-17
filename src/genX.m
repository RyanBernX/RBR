n = 100000;
p = 50;
K = 20;

X = randn(n, p);
blk = n / K;
for i=1:K
    X(1+(i-1)*blk:i*blk, :) = X(1+(i-1)*blk:i*blk, :) + i;
end
Xt = X';
fid = fopen('X.dat', 'w');
fprintf(fid, '%e\n', Xt(:));
fclose(fid);