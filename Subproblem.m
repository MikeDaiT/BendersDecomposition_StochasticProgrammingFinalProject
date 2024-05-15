function [status, y, pi] = Subproblem(x, dim_y, vec_f, mat_B, mat_D, vec_d)
%% Solve the subproblem
% Inputs:
% x, solution of the master problem;
% dim_y, 1*3 vector, dimentions of corresponding variable types of second stage decision variables, quantities for rational, integer, 0-1 variables respectively;
% vec_f, second stage cost vectors;
% mat_B, second stage constraint matrix for first stage decision variables;
% mat_D, second stage constraint matrices for second stage decision variables;
% vec_d, second stage constraint vectors;
% Return:
% status, value 1 means subproblem is not feasible, value 0 means subproblem is feasible;
% y, optimal solution of the subproblem. If the subproblem is not feasible, we set y as an all zero vector;
% pi, if the subproblem is feasible, it's the dual solution of the subproblem; If the subproblem is not feasible, pi is the dual solution of the problem SUB-PI.
%% Data Organization
% Second stage decision variables
if dim_y(1) ~= 0; y_r = sdpvar(dim_y(1), 1); end
if dim_y(2) ~= 0; y_i = intvar(dim_y(2), 1); end
if dim_y(3) ~= 0; y_b = binvar(dim_y(3), 1); end 

%% Solve the SUB-PI
s_p = sdpvar(size(vec_d, 1), 1);
s_n = sdpvar(size(vec_d, 1), 1);
e = ones(size(vec_d));

Obj_PI = e'*(s_p + s_n); % SUB-PI objective
Ctr = 0;
if dim_y(1) ~= 0; Ctr = Ctr + mat_D(:, 1:dim_y(1))*y_r; end
if dim_y(2) ~= 0; Ctr = Ctr + mat_D(:, dim_y(1)+1:dim_y(1)+dim_y(2))*y_i; end
if dim_y(3) ~= 0; Ctr = Ctr + mat_D(:, dim_y(1)+dim_y(2)+1:end)*y_b; end
F_set_PI = [Ctr + s_p - s_n == mat_B*x + vec_d, s_p >= 0, s_n >= 0];
if dim_y(1) ~= 0; F_set_PI = [F_set_PI, y_r >= 0]; end
if dim_y(2) ~= 0; F_set_PI = [F_set_PI, y_i >= 0]; end
if dim_y(3) ~= 0; F_set_PI = [F_set_PI, y_b >= 0]; end % SUB-PI constraint
options = sdpsettings('verbose', 0, 'solver', 'gurobi');
result_PI = optimize(F_set_PI, Obj_PI, options);
if result_PI.problem == 0
    if value(Obj_PI) > 0
        status = 1; % subproblem is infeasible
        y = zeros(dim_y(1)+dim_y(2)+dim_y(3), 1);
        pi = dual(F_set_PI(1));
        if isnan(pi) ~= zeros(size(pi))
            disp('aka');
            pi = zeros(size(pi));
            if abs(pi'*(mat_B*x + vec_d) + value(Obj_PI)) < 0.001*abs(value(Obj_PI))
                pi = -pi;
            end
        end
    else
        status = 0; % subproblem is feasible

%% Solve the Subproblem
        Obj = 0;
        if dim_y(1) ~= 0; Obj = Obj + vec_f(1:dim_y(1))'*y_r; end
        if dim_y(2) ~= 0; Obj = Obj + vec_f(dim_y(1)+1:dim_y(1)+dim_y(2))'*y_i; end
        if dim_y(3) ~= 0; Obj = Obj + vec_f(dim_y(1)+dim_y(2)+1:end)'*y_b; end
        F_set = [Ctr == mat_B*x + vec_d];
        if dim_y(1) ~= 0; F_set = [F_set, y_r >= 0]; end
        if dim_y(2) ~= 0; F_set = [F_set, y_i >= 0]; end
        if dim_y(3) ~= 0; F_set = [F_set, y_b >= 0]; end
        result = optimize(F_set, Obj, options);
        if result.problem == 0
            y = [];
            if dim_y(1) ~= 0; y = [y; value(y_r)]; end 
            if dim_y(2) ~= 0; y = [y; value(y_i)]; end
            if dim_y(3) ~= 0; y = [y; value(y_b)]; end
            pi = dual(F_set(1));
            if abs(pi'*(mat_B*x + vec_d) + value(Obj)) < 0.001*abs(value(Obj))
                pi = -pi;
            end
            if isnan(pi) ~= zeros(size(pi))
                disp('aka');
                pi = zeros(size(pi));
            end
        else
            disp('Subproblem failed.');
            disp(result.info);
        end
    end
else
    disp('SUB-PI failed.');
    disp(result_PI.info);
end

end