function [Opt_value, x_star, time] = Benders(dim_x, dim_y, vec_c, cell_f, vec_p, mat_A, vec_b, cell_B, cell_D, cell_d, eps, IteMax)
%% Benders decomposition ultra ultimate script;
% Problem Form:
%   min c'x+E[f'y]
%   s.t.Ax>=b
%       Dy==Bx+d
%       x>=0, y>=0
% Inputs: 
%   the parameters of the problems with corresponding vector or matrix forms;
%   first stage and second stage decision variables with specified varaible catagories.
% Return: 
%   optimal value; 
%   optimal first stage variables.
%% Data Organization
% Problem Data
%   dim_x, dim_y, 1*3 vector, dimentions of corresponding variable types of first stage and second stage decision variables, quantities for rational, integer, 0-1 variables respectively;
%   NOTICE: all constraints should arrange their coefficients conformed the order of the decision variables, i.e. rational, integer, 0-1 variables from left to right;
%   vec_c, first stage cost vector;
%   cell_f, 1*s cell of second stage cost vectors;
%   vec_p, s-dimention senarios' distribution;
%   mat_A, first stage constraint matrix;
%   vec_b, first stage constraint vector;
%   cell_B, 1*s cell of second stage constraint matrices for first stage decision variables;
%   cell_D, 1*s cell of second stage constraint matrices for second stage decision variables;
%   cell_d, 1*s cell of second stage constraint vectors;

% Computing Parameter
%   eps, computing accuracy 
%   IteMax, maximal iteration time
scr_n = size(vec_p, 1); % Scenario quantities

% Decision Variables
%   first stage decision variables
if dim_x(1) ~= 0; x_r = sdpvar(dim_x(1), 1); end
if dim_x(2) ~= 0; x_i = intvar(dim_x(2), 1); end
if dim_x(3) ~= 0; x_b = binvar(dim_x(3), 1); end
% %   second stage decision variables
% if dim_y(1) ~= 0; y_r = sdpvar(dim_y(1), 1); end
% if dim_y(2) ~= 0; y_i = intvar(dim_y(2), 1); end
% if dim_y(3) ~= 0; y_b = binvar(dim_y(3), 1); end


%% Benders Decomposition
% We generally use multi-cut here. 
% We try 2 modes: 
%   Mode 1: calculate all subproblems with one iteration first stage variable;
%   Mode 2: generate one cut then calculate one new first stage variable immediately.
mode = 1; 

max_cuts = 0; % 先不设置这个



% Master Problem
upp_z = 10000000; % track the upper bound of optimal value
low_z = -10000000; % the lower bound of optimal value
theta = sdpvar(scr_n, 1); % master second stage optimal value upper bound variable

Obj = 0; 
if dim_x(1) ~= 0; Obj = Obj + vec_c(1:dim_x(1))'*x_r; end
if dim_x(2) ~= 0; Obj = Obj + vec_c(dim_x(1)+1:dim_x(1)+dim_x(2))'*x_i; end
if dim_x(3) ~= 0; Obj = Obj + vec_c(dim_x(1)+dim_x(2)+1:end)'*x_b; end
Obj_m = Obj + vec_p'*theta; % objective function of the master problem

Ctr = 0;
if dim_x(1) ~= 0; Ctr = Ctr + mat_A(:, 1:dim_x(1))*x_r; end
if dim_x(2) ~= 0; Ctr = Ctr + mat_A(:, dim_x(1)+1:dim_x(1)+dim_x(2))*x_i; end
if dim_x(3) ~= 0; Ctr = Ctr + mat_A(:, dim_x(1)+dim_x(2)+1:end)*x_b; end

% original master problem constraints
F_set_M_O = [Ctr >= vec_b]; 
if dim_x(1) ~= 0; F_set_M_O = [F_set_M_O, x_r >= 0]; end
if dim_x(2) ~= 0; F_set_M_O = [F_set_M_O, x_i >= 0]; end
if dim_x(3) ~= 0; F_set_M_O = [F_set_M_O, x_b >= 0]; end
F_set_M = F_set_M_O; % master problem constraints
Cut_set = []; % initial cut set

tstart=tic; % solving starts
% Solve
Ite = 0;
while 1
    Ite = Ite + 1;
    options = sdpsettings('verbose', 0, 'solver', 'gurobi');
    % solve the master problem
    result = optimize(F_set_M, Obj_m, options); 
    z = value(Obj_m);
    if result.problem == 12 % if the master problem is infeasible or unbounded
        result = optimize(F_set_M_O, Obj, options);
        z = low_z;
    end
    low_z = z;
    
    if result.problem == 0
        x_mid = []; % master decision variable
        if dim_x(1) ~= 0; x_r_mid = value(x_r); x_mid = [x_mid, x_r_mid]; end
        if dim_x(2) ~= 0; x_i_mid = value(x_i); x_mid = [x_mid, x_i_mid]; end
        if dim_x(3) ~= 0; x_b_mid = value(x_b); x_mid = [x_mid, x_b_mid]; end
        z_hat = vec_c'*x_mid; % current iteration solution
        trigger = 0; % update upper bound trigger
            
        % solve subproblems
        for i = 1:scr_n
            [status, y, pi] = Subproblem(x_mid, dim_y, cell_f{i}, cell_B{i}, cell_D{i}, cell_d{i});
            pi_B = pi'*cell_B{i}; pi_B_x = 0;
            if dim_x(1) ~= 0; pi_B_x = pi_B_x + pi_B(1:dim_x(1))*x_r; end
            if dim_x(2) ~= 0; pi_B_x = pi_B_x + pi_B(dim_x(1)+1:dim_x(1)+dim_x(2))*x_i; end
            if dim_x(3) ~= 0; pi_B_x = pi_B_x + pi_B(dim_x(1)+dim_x(2)+1:end)*x_b; end
            if status == 0 % generate optimality cut
                Cut_set = [Cut_set, theta(i) >= pi_B_x + pi'*cell_d{i}];
                if pi_B*x_mid + pi'*cell_d{i} ~= cell_f{i}'*y
                    disp('Cut is wrong.');
                end
                if trigger == 0
                    z_hat = z_hat + vec_p(i)*cell_f{i}'*y;
                end
            else % generate feasible cut
                Cut_set = [Cut_set, 0 >= pi_B_x + pi'*cell_d{i}];
                trigger = 1; % some subproblem is infeasible, no upper bound update at this master iteration
            end
        end
        if trigger == 0 && z_hat < upp_z % update upper bound
            upp_z = z_hat;
        end
        x_star = x_mid;
        Opt_value = low_z; 
        if max_cuts > 0
            % cuts dropping, 之后再说
        end
        % iteration accuracy
        % disp('---------------------------------');
        % disp(upp_z);
        % disp(low_z);
        % disp((upp_z - low_z)/min(abs(upp_z), abs(low_z)));
        % add cuts
        F_set_M = [F_set_M_O, Cut_set];
    else
        disp('Master problem failed.');
        disp(Ite);
        disp(result.info);
        break;
    end

    % iteration stopping criteria
    if upp_z - low_z <= eps*min(abs(upp_z), abs(low_z)) || Ite == IteMax
        % debug
        if upp_z - low_z < -0.001
            disp('Something wrong, upper bound is smaller than lower bound.');
        end
        % debug
        if upp_z - low_z <= eps*min(abs(upp_z), abs(low_z))
            disp('Obtain the accuracy.');
        elseif Ite == IteMax
            disp('Get the iteration threshold.');
        end
        disp(Ite);
        break;
    end

end
time=toc(tstart); % solving ends


end