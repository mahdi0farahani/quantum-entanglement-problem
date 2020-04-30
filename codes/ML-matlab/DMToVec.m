function vec = DMToVec(H)
	d = size(H,1);
	vec = [DiagTerm(H,d); NonDiagTerm(H,d)];
	vec = vec * sqrt(d / (d+d-2));
end

function vec = DiagTerm(H,d)
	vec = zeros(d-1,1);
	tmp = 0;
	for j = 1:d-1
		tmp = tmp + H(j,j);
		vec(j) = sqrt(2 / (j*(j+1))) * real(tmp - j * H(j+1,j+1));
	end
end

function vec = NonDiagTerm(H,d)
	vec = zeros(d*(d-1), 1);
	tot = 0;
	for j = 1:d-1
		for k = j+1:d
			tot = tot+1;
			vec(tot) = real(H(j,k) + H(k,j));
			tot = tot+1;
			vec(tot) = imag(H(j,k) - H(k,j));
		end
	end
end
