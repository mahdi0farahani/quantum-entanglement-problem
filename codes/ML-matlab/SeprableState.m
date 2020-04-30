function [ S ] = SeprableState( dim1 , dim2 , index )

%dim = 1;
%lambda = ones(1,index);
addpath(genpath('C:\Users\nazanin1\Desktop\QETLAB-0.9'));
coef = drchrnd(1, index);

%dim1 = length( density_matrix_list1(1,:) );
%dim2  = length( density_matrix_list2(1,:) );

S = zeros(dim1*dim2,dim1*dim2);

for i = 1:index
    %x = randsample(index,1);
    %y = randsample(index,1);
    DM1 = RandomDensityMatrix(dim1);
    DM2 = RandomDensityMatrix(dim2);
    DM = kron(DM1,DM2);
    DM = coef(i) .* DM;
    S = S + DM;
end
end

function r = drchrnd(a,n)
	% take a sample from a dirichlet distribution
	r = gamrnd(repmat(a,n,1),1,n,1);
	r = r / sum(r,1);
end

