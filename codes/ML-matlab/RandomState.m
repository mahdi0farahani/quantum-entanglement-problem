function rho = RandomState(lambda, dim)
	addpath(genpath('C:\Users\nazanin1\Desktop\QETLAB-0.9'));
	d = diag(drchrnd(lambda, dim));
	U = RandomUnitary(dim);
	rho = U * d * U';
end

function r = drchrnd(a,n)
	% take a sample from a dirichlet distribution
	r = gamrnd(repmat(a,n,1),1,n,1);
	r = r / sum(r,1);
end




