function BCHM(dim, m)
	if dim == 2
		file_rdm	= '2x2rdm.mat';
		file_alpha 	= sprintf('qubit%d.mat', m);
	else
		file_rdm	= '3x3rdm.mat';
		file_alpha 	= sprintf('qutrit%d.mat', m);
	end
	load(file_rdm,	'points', 'label');
	load(file_alpha, 'alpha');

	n = size(points, 1) / 2;
	X_train = [points(1:n,:) alpha(1:n)];
	X_test  = [points(n+1:2*n,:), alpha(n+1:2*n)];
	y_train = label(1:n);
	y_test	= label(n+1:2*n);

	y_fit = alpha >= 1;
	y_fit = y_fit(n+1:2*n);
	acc = sum(y_fit == y_test) / n;
	disp(acc);

	B = TreeBagger(100, X_train, y_train);
	y_fit = cell2mat(predict(B, X_test)) - 48;  %Character to Integer
	acc = sum(y_fit == y_test) / n;
	disp(acc);
end
