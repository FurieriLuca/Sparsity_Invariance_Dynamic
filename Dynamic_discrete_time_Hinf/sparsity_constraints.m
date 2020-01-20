



CQYv = sdpvar(m,n*N);          % decision variables
DQYv = sdpvar(m,n);



T_bin_aug = kron(ones(1,N),Tbin);
fprintf('Encoding the sparsity constraints ...\n')
CQYv = CQYv.*T_bin_aug;
DQYv = DQYv.*Tbin;
