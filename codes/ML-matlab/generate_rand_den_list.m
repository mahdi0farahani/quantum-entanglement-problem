function [ den_list ] = generate_rand_den_list( dim , lenght )

den_list = {};
for i = 1:lenght
    R = RandomDensityMatrix(dim);
    den_list = [den_list,R];
end

end
