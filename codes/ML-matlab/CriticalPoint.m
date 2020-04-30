function points = CriticalPoint(rho, d1, d2, maxiter, rate)
	vrho = DMToVec(rho);
	dim = d1*d1*d2*d2 - 1;
	init = dim * 100;
	kets = zeros(d1 * d2, init);
	points = zeros(init, dim);
	parfor i = 1 : init
		kets(:, i) = kron(RandomStateVector(d1), RandomStateVector(d2));
		points(i, :) = DMToVec(kets(:, i) * kets(:, i)');
	end
	
	[ans1, x] = CompAlpha(vrho, points); 
	kets = kets(:, x);
	points = points(x, :);

	eps = 1;
	%disp(ans);
	iter = 0;
	while 1
		if ans1 >= 1
			break;
		end
		if iter == maxiter
			break;
		end
		iter = iter + 1;
		n = size(points, 1);
		for i = 1 : n
			for j = 1: 10
				U1 = Random_Unitary(d1, rand * eps);
				U2 = Random_Unitary(d2, rand * eps);
				ket = kron(U1, U2) * kets(:, i);
				kets = [kets ket];
				points = [points; transpose(DMToVec(ket * ket'))];
			end
		end
		[ans2, x] = CompAlpha(vrho, points);
		kets = kets(:, x);
		points = points(x, :);
		%msg = sprintf('p_c = %.8f', ans);
	%	disp(msg);
		eps = eps * rate;
	end
	disp(ans2);
end

function U = Random_Unitary(d, eps)
	H = Random_Hermitian(d);
	U = expm(eps * 1j * H);
end

function H = Random_Hermitian(d)
	vec = randn(d*d);
	vec = vec / norm(vec);
	H = zeros(d);
	for j = 1: d-1
		a = vec(j) * sqrt(2.0 / double(j*(j+1)));
		for k = 1 : j
			H(k, k) = H(k, k) + a;
		end
		H(j+1, j+1) = -a * double(j);
	end
	tot = d;
	for j = 1:d-1
		for k = j+1:d
			H(j,k) = vec(tot) + vec(tot+1) * 1j;
			H(k,j) = vec(tot) - vec(tot+1) * 1j;
			tot = tot+2;
		end
	end
	H = H / sqrt(2);
	for i = 1 : d
		H(i, i) = H(i, i) + vec(tot) / sqrt(d);
	end
end
