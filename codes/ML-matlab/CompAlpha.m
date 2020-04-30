function [a, x] = CompAlpha(ray, extreme)
    d = size(extreme, 2);%dim of feature vector
    n = size(extreme, 1);%num of points
  	ray = reshape(ray, 1, d);
    %parameteres for linprog
    f = zeros(n+1,1);
    f(1) = -1;
    Aeq = [transpose([-ray; extreme]); 0 ones(1,n)];%ino nemifahmim
    beq = [zeros(d,1); 1];
    lb = zeros(n+1,1);

    opt = optimoptions(@linprog, 'Algorithm', 'dual-simplex', 'Display','final');
    [x, fval] = linprog(f,[],[],Aeq,beq,lb,[],[],opt);
	a = -fval;
	x = find(x > 0) - 1;%nmidunim chera
	x = x(2 : size(x,1));
end
