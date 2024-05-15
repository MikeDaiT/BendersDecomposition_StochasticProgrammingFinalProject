dim_x = [3; 0; 0];
dim_y = [10; 0; 0];
vec_c = [150; 230; 260];
vec_f = [238; 210; -170; -150; -36; -10; zeros(4, 1)];
cell_f = {};
for i = 1:3
    cell_f{end+1}=vec_f;
end
vec_p = 1/3*ones(3, 1);
mat_A = -ones(1, 3);
vec_b = -500;
mat_D = zeros(4, 6);
mat_D(1, 1) = 1;
mat_D(2, 2) = 1;
mat_D(1, 3) = -1;
mat_D(2, 4) = -1;
mat_D(3, 5) = -1;
mat_D(3, 6) = -1;
mat_D(4, 5) = -1;
cell_D = {};
for i = 1:3
    cell_D{end+1}=[mat_D, -eye(4)];
end
cell_B = {[-3 0 0; 0 -3.6 0; 0 0 -24; 0 0 0], [-2.5 0 0; 0 -3 0; 0 0 -20; 0 0 0], [-2 0 0; 0 -2.4 0; 0 0 -16; 0 0 0]};
cell_d = {[200; 240; 0; -6000], [200; 240; 0; -6000], [200; 240; 0; -6000]};
% solve
[Opt_value, x_star, time] = Benders(dim_x, dim_y, vec_c, cell_f, vec_p, mat_A, vec_b, cell_B, cell_D, cell_d, 0.00000001, 15);
disp(time);