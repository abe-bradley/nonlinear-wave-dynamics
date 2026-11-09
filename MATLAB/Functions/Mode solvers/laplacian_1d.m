function A = laplacian_1d(N)
% Builds dx^2 times the Laplacian in 1D with Dirichlet boundary condition
% N = number of grid points

% note this A is stored as a full matrix
A = diag(-2*ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);

end
