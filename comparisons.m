T = 20; % number of iterations
[n,m] = size(fea); % input data
k = length(unique(gnd)); % number of clusters

% sequential version
[Ac1,W0] = ours_s(fea,gnd,k,50,.6,rand_ortho(n,k));

% batch-mode-1
[Ac2] = ours_p(fea,gnd,k,50,.6,rand_ortho(n,k));

% batch-mode-2
[Ac3] = ours_p1(fea,gnd,k,50,.6,rand_ortho(n,k));


