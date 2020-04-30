clc
clear

%Matlab test
addpath(genpath('C:\Users\nazanin1\Desktop\QETLAB-0.9'));
d1 = 3;
d2 = 3;
n= 0;
%ket1 = kron(RandomStateVector(d1), RandomStateVector(d2));
%ket1 = eye(3,3)./2;
%rho = ket1 * ket1';
%rho = RandomDensityMatrix(9);
%list1 = generate_rand_den_list(d1,10);
%list2 = generate_rand_den_list(d2,10);
for z = 1:500
    rho  = SeprableState(d1,d2,8);
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
    maxiter = 10;
    rate = 1;
    while 1
        if ans1 >= 1
            break;
        end
        if iter == maxiter
            n = n+1;
            break;
        end
        iter = iter + 1;
        n = size(points, 1);
        tic
        for i = 1 : n
            for j = 1: 10
                U1 = Random_Unitary(d1, rand * eps);
                U2 = Random_Unitary(d2, rand * eps);
                ket = kron(U1, U2) * kets(:, i);
                kets = [kets ket];
                points = [points; transpose(DMToVec(ket * ket'))];
            end
        end
        toc
        %fprintf('WOW Its rotating');
        %disp(iter);
        [ans1, x] = CompAlpha(vrho, points);
        kets = kets(:, x);
        points = points(x, :);
        %msg = sprintf('p_c = %.8f', ans);
        %	disp(msg);
        eps = eps * rate;
    end
    %disp(ans1);
end
disp(n/500);



