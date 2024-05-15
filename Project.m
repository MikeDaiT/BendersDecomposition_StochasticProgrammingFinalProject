close all; clear; clc;

%% Project problem 

%% Data organization
addpath('Data/');
load('IPV1');  load('IPV2'); load('IPV3'); load('IPV4'); load('IPV5'); load('IPV6');   
load('Eload1');  load('Eload2'); load('Eload3'); load('Eload4'); load('Eload5'); load('Eload6');
It = [IPV1;IPV2;IPV3;IPV4;IPV5;IPV6];
Edemand = [Eload1;Eload2;Eload3;Eload4;Eload5;Eload6]; % electric demand
load('Nday'); 
vec_p = Nday / 365;
E_edemand = vec_p*Edemand;

%% Computing Parameter
%   eps, computing accuracy 
%   IteMax, maximal iteration time
scr_n = size(vec_p, 2); % Scenario quantities

% Decision Variables
%   first stage decision variables
A_PV = sdpvar(1, 1); % m^2, area of PV panel, in PV model (P)
E_EL_S =  sdpvar(1, 1); % kW, power of Electrolyzer, in P2G model (P)
E_CH4_S = sdpvar(1, 1); % kW, power of CH4 reaction tank, in P2G model (P)
V_H2_S = sdpvar(1, 1); % N*m^3, volume of H2 storage tank, in H2 storage model (P)
V_CO2_S = sdpvar(1, 1); % N*m^3, volume of CO2 storage tank, in CO2 storage model (P)
x = [A_PV, E_EL_S, E_CH4_S, V_H2_S, V_CO2_S];
 

%   Important second stage decision variables
%       PV model
E_PV = sdpvar(24, 1); % kW, real power of PV (P)
%       CCS model
E_PGU = sdpvar(24, 1); % kW, power of thermal power (P)
E_CCS = sdpvar(24, 1); % kW, power to capture CO2 (P)
V_CO2_CCS = sdpvar(24, 1); % N*m^3, true volume of captured CO2, (P)
%       P2G model
V_H2_EL = sdpvar(24, 1); % N*m^3, volume of generated H2, (P)
E_EL = sdpvar(24, 1); % kW, power used to generate H2, (P)
E_CH4 = sdpvar(24, 1); % kW, power used to generate CH4, (P)
V_CO2_P2G = sdpvar(24, 1); % N*m^3, volume of CO2 used to generate CH4, (P)
V_H2_P2G = sdpvar(24, 1); % N*m^3, volume of H2 used to generate CH4, (P)
%       H2 storage model
V_H2 = sdpvar(24, 1); % N*m^3, H2 quantity of H2 storage equipment, (P)
V_H2_in = sdpvar(24, 1); % N*m^3, H2 in, (P)
V_H2_out = sdpvar(24, 1); % N*m^3, H2 out, (P)
%       CO2 storage model
V_CO2 = sdpvar(24, 1); % N*m^3, CO2 quantity of CO2 storage equipment, (P)
V_CO2_in = sdpvar(24, 1); % N*m^3, CO2 in, (P)
V_CO2_out = sdpvar(24, 1); % N*m^3, CO2 out, (P)






%% Objectives and constraints of master problem and subproblems

max_cuts = 0; % 先不设置这个

% original master problem constraints
F_set_M_O = []; 


%% PV model
ita_PV = 0.200/1000; 
k = 0.200;
E_PVr_max = 300000; % kW, rated power generation of photovoltaic panels

% second stage decision variables
E_PV_cur = sdpvar(24, 1); % kW, discarded power of PV, (P)


F_set_M_O = [F_set_M_O, A_PV*k <= E_PVr_max];


% SUB-PI constraints
s_p_1 = sdpvar(24, 1); 
s_n_1 = sdpvar(24, 1);


%% CCS model
e_PGU = 0.46;
ita_CCS_max = 0.65;
lamdaCO2 = 0.1937;
E_PGUmax = 300000;  %kW
E_PGUmin = 90000;   %kW
E_CCSmax = 4.4551e+03; %kW
dita_E_PGU = 50000;


% second stage decision variables
V_CO2_cur = sdpvar(24, 1); % N*m^3, discarded volume of CO2, (P)
V_CO2_PGU = sdpvar(24, 1); % N*m^3, emission volume of CO2

s_E_PGU_1 = sdpvar(24, 1); % slackness (P)
s_E_PGU_2 = sdpvar(24, 1); % slackness (P)
s_E_PGU_3 = sdpvar(23, 1); % slackness (P)
s_E_PGU_4 = sdpvar(23, 1); % slackness (P)
s_E_CCS = sdpvar(24, 1); % slackness (P)


% SUB-PI constraints
s_p_2 = sdpvar(24, 1); 
s_n_2 = sdpvar(24, 1);
s_p_3 = sdpvar(24, 1);
s_n_3 = sdpvar(24, 1);
s_p_4 = sdpvar(24, 1);
s_n_4 = sdpvar(24, 1);
s_p_5 = sdpvar(24, 1);
s_n_5 = sdpvar(24, 1);
s_p_6 = sdpvar(24, 1);
s_n_6 = sdpvar(24, 1);
s_p_7 = sdpvar(24, 1);
s_n_7 = sdpvar(24, 1);
s_p_8 = sdpvar(23, 1);
s_n_8 = sdpvar(23, 1);
s_p_9 = sdpvar(23, 1);
s_n_9 = sdpvar(23, 1);



%% P2G model
lamda_H2  = 4.2;
lamda_CH4 = 0.3;
E_EL_Smax = 100000;
E_EL_Smin = 10000;



% second stage decision variables
V_CH4_P2G = sdpvar(24, 1); % volume of generating CH4


s_E_EL_1 = sdpvar(24, 1); % slackness (P)
s_E_EL_2 = sdpvar(24, 1); % slackness (P)
s_E_EL_3 = sdpvar(23, 1); % slackness (P)
s_E_EL_4 = sdpvar(23, 1); % slackness (P)
s_E_EL_5 = sdpvar(23, 1); % slackness (P)
s_E_EL_6 = sdpvar(23, 1); % slackness (P)
s_E_CH4_1 = sdpvar(24, 1); % slackness (P)
s_E_CH4_2 = sdpvar(23, 1); % slackness (P)
s_E_CH4_3 = sdpvar(23, 1); % slackness (P)


F_set_M_O = [F_set_M_O, E_EL_Smin <= E_EL_S, E_EL_S <= E_EL_Smax, E_CH4_S <= 2000];



% SUB-PI constraints
s_p_10 = sdpvar(24, 1); 
s_n_10 = sdpvar(24, 1);
s_p_11 = sdpvar(24, 1); 
s_n_11 = sdpvar(24, 1);
s_p_12 = sdpvar(24, 1); 
s_n_12 = sdpvar(24, 1);
s_p_13 = sdpvar(24, 1); 
s_n_13 = sdpvar(24, 1);
s_p_14 = sdpvar(24, 1); 
s_n_14 = sdpvar(24, 1);
s_p_15 = sdpvar(24, 1); 
s_n_15 = sdpvar(24, 1);
s_p_16 = sdpvar(23, 1); 
s_n_16 = sdpvar(23, 1);
s_p_17 = sdpvar(23, 1); 
s_n_17 = sdpvar(23, 1);
s_p_18 = sdpvar(23, 1); 
s_n_18 = sdpvar(23, 1);
s_p_19 = sdpvar(23, 1); 
s_n_19 = sdpvar(23, 1);
s_p_20 = sdpvar(24, 1); 
s_n_20 = sdpvar(24, 1);
s_p_21 = sdpvar(23, 1); 
s_n_21 = sdpvar(23, 1);
s_p_22 = sdpvar(23, 1); 
s_n_22 = sdpvar(23, 1);



%% H2 storage model

% second stage decision variables
s_H2_in = sdpvar(24, 1); % slackness (P)
s_H2_out = sdpvar(24, 1); % slackness (P)
s_H2 = sdpvar(24, 1); % slackness (P)


F_set_M_O = [F_set_M_O, V_H2_S<=20000];


% SUB-PI constraints
s_p_23 = sdpvar(23, 1); 
s_n_23 = sdpvar(23, 1);
s_p_24 = sdpvar(24, 1); 
s_n_24 = sdpvar(24, 1);
s_p_25 = sdpvar(24, 1); 
s_n_25 = sdpvar(24, 1);
s_p_26 = sdpvar(24, 1); 
s_n_26 = sdpvar(24, 1);
s_p_27 = sdpvar(1, 1); 
s_n_27 = sdpvar(1, 1);



%% CO2 storage model

% second stage decision variables
s_CO2_in = sdpvar(24, 1); % slackness (P)
s_CO2_out = sdpvar(24, 1); % slackness (P)
s_CO2 = sdpvar(24, 1); % slackness (P)


F_set_M_O = [F_set_M_O, V_CO2_S<=20000];                  


% SUB-PI constraints
s_p_28 = sdpvar(23, 1); 
s_n_28 = sdpvar(23, 1);
s_p_29 = sdpvar(24, 1); 
s_n_29 = sdpvar(24, 1);
s_p_30 = sdpvar(24, 1); 
s_n_30 = sdpvar(24, 1);
s_p_31 = sdpvar(24, 1); 
s_n_31 = sdpvar(24, 1);
s_p_32 = sdpvar(1, 1); 
s_n_32 = sdpvar(1, 1);



%% Energy balance


% SUB-PI constraints
s_p_33 = sdpvar(24, 1); 
s_n_33 = sdpvar(24, 1);


%% H2 balance


% SUB-PI constraints
s_p_34 = sdpvar(24, 1); 
s_n_34 = sdpvar(24, 1);


%% C balance


% SUB-PI constraints
s_p_35 = sdpvar(24, 1); 
s_n_35 = sdpvar(24, 1);


%% Constraints complement

F_set_PoADD = [A_PV>=0, E_EL_S>=0, E_CH4_S>=0, V_H2_S>=0, V_CO2_S>=0];
F_set_M_O = [F_set_M_O, F_set_PoADD];


C_sub_PoADD = [E_PV>=0, E_PGU>=0, E_CCS>=0, V_CO2_CCS>=0,...
                V_H2_EL>=0, E_EL>=0, E_CH4>=0, V_CO2_P2G>=0, V_H2_P2G>=0,...
                V_H2>=0, V_H2_in>=0, V_H2_out>=0, V_CO2>=0, V_CO2_in>=0, V_CO2_out>=0,...
                E_PV_cur>=0, V_CO2_cur>=0, V_CO2_PGU>=0,...
                s_E_PGU_1>=0, s_E_PGU_2>=0, s_E_PGU_3>=0, s_E_PGU_4>=0, s_E_CCS>=0,...
                V_CH4_P2G>=0, s_E_EL_1>=0, s_E_EL_2>=0, s_E_EL_3>=0, s_E_EL_4>=0,...
                s_E_EL_5>=0, s_E_EL_6>=0, s_E_CH4_1>=0, s_E_CH4_2>=0, s_E_CH4_3>=0,...
                s_H2_in>=0, s_H2_out>=0, s_H2>=0, s_CO2_in>=0, s_CO2_out>=0, s_CO2>=0,...
                ];



s_p = [s_p_1; s_p_2; s_p_3; s_p_4; s_p_5; s_p_6; s_p_7; s_p_8; s_p_9; s_p_10; ...
                    s_p_11; s_p_12; s_p_13; s_p_14; s_p_15; s_p_16; s_p_17; s_p_18; s_p_19; s_p_20; ...
                    s_p_21; s_p_22; s_p_23; s_p_24; s_p_25; s_p_26; s_p_27; s_p_28; s_p_29; s_p_30; ...
                    s_p_31; s_p_32; s_p_33; s_p_34; s_p_35];
s_n = [s_n_1; s_n_2; s_n_3; s_n_4; s_n_5; s_n_6; s_n_7; s_n_8; s_n_9; s_n_10; ...
                    s_n_11; s_n_12; s_n_13; s_n_14; s_n_15; s_n_16; s_n_17; s_n_18; s_n_19; s_n_20; ...
                    s_n_21; s_n_22; s_n_23; s_n_24; s_n_25; s_n_26; s_n_27; s_n_28; s_n_29; s_n_30; ...
                    s_n_31; s_n_32; s_n_33; s_n_34; s_n_35];


%% Objectives
dita_PV = 0.1*1.1^20/(1.1^20-1);
dita_EL = 0.1*1.1^10/(1.1^10-1);
dita_CH4 = 0.1*1.1^20/(1.1^20-1);
dita_H2 = 0.1*1.1^25/(1.1^25-1);
dita_CO2 = 0.1*1.1^15/(1.1^15-1);

cinv_PV = 500; % yuan/kW
cinv_EL = 100; % yuan/kW
cinv_CH4 = 300; % yuan/kW
cinv_H2  = 7.76; % yuan/N.m3
cinv_CO2 = 0.776; % yuan/N.m3

% electric price
E_Pri = [0.382*ones(8, 1); 0.54*ones(4, 1); 0.922*ones(4, 1); 0.54*ones(3, 1); 0.922*ones(3, 1); 0.54*ones(2, 1)]; 


theta = sdpvar(scr_n, 1); % master second stage optimal value upper bound variable

% objective function of the master problem
C_inv = cinv_PV*A_PV*k*dita_PV+cinv_EL*E_EL_S*dita_EL+cinv_CH4*E_CH4_S*dita_CH4...
                 +cinv_H2*V_H2_S*dita_H2+cinv_CO2*V_CO2_S*dita_CO2;
C_OM = cinv_PV*A_PV*k*0.03+cinv_EL*E_EL_S*0.04+cinv_CH4*E_CH4_S*0.05...
                 +cinv_H2*V_H2_S*0.02+cinv_CO2*V_CO2_S*0.02;  
Obj = C_inv + C_OM;
Obj_m = Obj + vec_p*theta; 



%% Benders Decomposition
% We generally use multi-cut here. 

upp_z = 1.00e+18; % track the upper bound of optimal value
low_z = -1.00e+18; % the lower bound of optimal value

F_set_M = F_set_M_O; % master problem constraints
Cut_set = []; % initial cut set

tstart=tic; % solving starts
% Solve
Ite = 0;
IteMax = 40;
eps = 0.00001;
while 1
    Ite = Ite + 1;
    options = sdpsettings('verbose', 0, 'solver', 'gurobi');
    % solve the master problem
    result = optimize(F_set_M, Obj_m, options); 
    z = value(Obj_m) - E_edemand*E_Pri;
    if result.problem == 12 % if the master problem is infeasible or unbounded
        result = optimize(F_set_M, Obj, options);
        z = low_z;
    end
    low_z = z;
    
    if result.problem == 0
        z_hat = value(Obj); 
        trigger = 0; % update upper bound trigger

        % first stage solution
        x_mid = [value(A_PV), value(E_EL_S), value(E_CH4_S), value(V_H2_S), value(V_CO2_S)];
        
        status = 0;
        %% solve subproblems
        for i = 1:scr_n
            it = It(i, :)'; edemand = Edemand(i, :)';
            A_PV_v = value(A_PV); E_EL_S_v = value(E_EL_S); E_CH4_S_v = value(E_CH4_S);
            V_H2_S_v = value(V_H2_S); V_CO2_S_v = value(V_CO2_S);

            % objective function of subproblem
            C_tari = sum(350*E_PGU*320/10000-E_Pri.*edemand- 500*V_CH4_P2G);
            
            C_tax = sum((V_CO2_PGU-V_CO2_CCS)*0.6);
            
            Obj_sub = C_tari + C_tax;

            % subproblem constraints
            C_sub = [];

            C_sub = [C_sub, E_PV + E_PV_cur == ita_PV*it*A_PV_v];

            C_sub = [C_sub, e_PGU*E_PGU - V_CO2_PGU == 0,...
                V_CO2_CCS + V_CO2_cur - ita_CCS_max*V_CO2_PGU == 0,...
                E_CCS - lamdaCO2*V_CO2_CCS == 0,...
                E_PGU - s_E_PGU_1 == E_PGUmin,...
                E_PGU + s_E_PGU_2 == E_PGUmax,...
                E_CCS + s_E_CCS == E_CCSmax,...
                E_PGU(2:24)-E_PGU(1:23) - s_E_PGU_3 == -dita_E_PGU,...
                E_PGU(2:24)-E_PGU(1:23) + s_E_PGU_4 == dita_E_PGU
                ];

            C_sub = [C_sub, E_EL - lamda_H2*V_H2_EL == 0,...
                E_CH4 - lamda_CH4*V_CH4_P2G == 0,...
                V_CH4_P2G - V_CO2_P2G == 0,...
                V_CH4_P2G - V_H2_P2G/4 == 0,...
                E_EL - s_E_EL_1 == 0.3*E_EL_S_v,...
                E_EL + s_E_EL_2 == E_EL_S_v,...
                E_EL(2:24)-E_EL(1:23)-s_E_EL_3 == -0.4*E_EL_S_v,...
                E_EL(2:24)-E_EL(1:23)+s_E_EL_4 == 0.2*E_EL_S_v,...
                E_EL(1:23)-E_EL(2:24)-s_E_EL_5 == -0.4*E_EL_S_v,...
                E_EL(1:23)-E_EL(2:24)+s_E_EL_6 == 0.2*E_EL_S_v,...
                E_CH4 + s_E_CH4_1 == E_CH4_S_v,...
                E_CH4(1:23)-E_CH4(2:24)-s_E_CH4_2 == -0.3*E_CH4_S_v,...
                E_CH4(1:23)-E_CH4(2:24)+s_E_CH4_3 == 0.3*E_CH4_S_v
                ];

            C_sub = [C_sub, V_H2(2:24)-V_H2(1:23)-V_H2_in(1:23)+V_H2_out(1:23) == 0,...
                V_H2_in+s_H2_in == 0.125*V_H2_S_v,...
                V_H2_out+s_H2_out == 0.125*V_H2_S_v,...
                V_H2+s_H2 == V_H2_S_v,...
                V_H2(1) - V_H2(24) == 0
                ];

            C_sub = [C_sub, V_CO2(2:24)-V_CO2(1:23)-V_CO2_in(1:23)+V_CO2_out(1:23)==0,...
                V_CO2_in+s_CO2_in == 0.125*V_CO2_S_v,...
                V_CO2_out+s_CO2_out == 0.125*V_CO2_S_v,...
                V_CO2+s_CO2 == V_CO2_S_v,...
                V_CO2(1)-V_CO2(24)==0
                ];

            C_sub = [C_sub, E_PV + E_PGU - E_CCS - E_EL - E_CH4 == edemand];

            C_sub = [C_sub, V_H2_EL+V_H2_out-V_H2_P2G-V_H2_in == 0];

            C_sub = [C_sub, V_CO2_CCS+V_CO2_out-V_CO2_P2G-V_CO2_in == 0];

            C_sub = [C_sub, C_sub_PoADD];

            % SUB-PI constraints
            C_sub_PI = [];

            C_sub_PI = [C_sub_PI, E_PV + E_PV_cur + s_p_1 - s_n_1 == ita_PV*it*A_PV_v];

            C_sub_PI = [C_sub_PI, e_PGU*E_PGU - V_CO2_PGU + s_p_2 - s_n_2 == 0,...
                V_CO2_CCS + V_CO2_cur - ita_CCS_max*V_CO2_PGU + s_p_3 - s_n_3 == 0,...
                E_CCS - lamdaCO2*V_CO2_CCS + s_p_4 - s_n_4 == 0,...
                E_PGU - s_E_PGU_1 + s_p_5 - s_n_5 == E_PGUmin,...
                E_PGU + s_E_PGU_2 + s_p_6 - s_n_6 == E_PGUmax,...
                E_CCS + s_E_CCS + s_p_7 - s_n_7 == E_CCSmax,...
                E_PGU(2:24)-E_PGU(1:23) - s_E_PGU_3 + s_p_8 - s_n_8 == -dita_E_PGU,...
                E_PGU(2:24)-E_PGU(1:23) + s_E_PGU_4 + s_p_9 - s_n_9 == dita_E_PGU
                ];

            C_sub_PI = [C_sub_PI, E_EL - lamda_H2*V_H2_EL + s_p_10 - s_n_10 == 0,...
                E_CH4 - lamda_CH4*V_CH4_P2G + s_p_11 - s_n_11 == 0,...
                V_CH4_P2G - V_CO2_P2G + s_p_12 - s_n_12 == 0,...
                V_CH4_P2G - V_H2_P2G/4 + s_p_13 - s_n_13 == 0,...
                E_EL - s_E_EL_1 + s_p_14 - s_n_14 == 0.3*E_EL_S_v,...
                E_EL + s_E_EL_2 + s_p_15 - s_n_15 == E_EL_S_v,...
                E_EL(2:24)-E_EL(1:23)-s_E_EL_3 + s_p_16 - s_n_16 == -0.4*E_EL_S_v,...
                E_EL(2:24)-E_EL(1:23)+s_E_EL_4 + s_p_17 - s_n_17 == 0.2*E_EL_S_v,...
                E_EL(1:23)-E_EL(2:24)-s_E_EL_5 + s_p_18 - s_n_18 == -0.4*E_EL_S_v,...
                E_EL(1:23)-E_EL(2:24)+s_E_EL_6 + s_p_19 - s_n_19 == 0.2*E_EL_S_v,...
                E_CH4 + s_E_CH4_1 + s_p_20 - s_n_20 == E_CH4_S_v,...
                E_CH4(1:23)-E_CH4(2:24)-s_E_CH4_2 + s_p_21 - s_n_21 == -0.3*E_CH4_S_v,...
                E_CH4(1:23)-E_CH4(2:24)+s_E_CH4_3 + s_p_22 - s_n_22 == 0.3*E_CH4_S_v
                ];

            C_sub_PI = [C_sub_PI, V_H2(2:24)-V_H2(1:23)-V_H2_in(1:23)+V_H2_out(1:23) + s_p_23 - s_n_23 == 0,...
                V_H2_in+s_H2_in + s_p_24 - s_n_24 == 0.125*V_H2_S_v,...
                V_H2_out+s_H2_out + s_p_25 - s_n_25 == 0.125*V_H2_S_v,...
                V_H2+s_H2 + s_p_26 - s_n_26 == V_H2_S_v,...
                V_H2(1) - V_H2(24) + s_p_27 - s_n_27 == 0
                ];

            C_sub_PI = [C_sub_PI, V_CO2(2:24)-V_CO2(1:23)-V_CO2_in(1:23)+V_CO2_out(1:23) + s_p_28 - s_n_28 == 0,...
                V_CO2_in+s_CO2_in + s_p_29 - s_n_29 == 0.125*V_CO2_S_v,...
                V_CO2_out+s_CO2_out + s_p_30 - s_n_30 == 0.125*V_CO2_S_v,...
                V_CO2+s_CO2 + s_p_31 - s_n_31 == V_CO2_S_v,...
                V_CO2(1)-V_CO2(24) + s_p_32 - s_n_32 ==0
                ];

            C_sub_PI = [C_sub_PI, E_PV + E_PGU - E_CCS - E_EL - E_CH4 + s_p_33 - s_n_33 == edemand];

            C_sub_PI = [C_sub_PI, V_H2_EL+V_H2_out-V_H2_P2G-V_H2_in + s_p_34 - s_n_34 == 0];

            C_sub_PI = [C_sub_PI, V_CO2_CCS+V_CO2_out-V_CO2_P2G-V_CO2_in + s_p_35 - s_n_35 == 0];

            C_sub_PI = [C_sub_PI, C_sub_PoADD, s_p >= 0, s_n >= 0];

            % subproblem mat_B and vec_d
            mat_B = [ita_PV*it, zeros(24, 4);...
                    zeros(286, 5);...
                    zeros(24, 1), 0.3*ones(24, 1), zeros(24, 3);...
                    zeros(24, 1), ones(24, 1), zeros(24, 3);...
                    zeros(23, 1), -0.4*ones(23, 1), zeros(23, 3);...
                    zeros(23, 1), 0.2*ones(23, 1), zeros(23, 3);...
                    zeros(23, 1), -0.4*ones(23, 1), zeros(23, 3);...
                    zeros(23, 1), 0.2*ones(23, 1), zeros(23, 3);...
                    zeros(24, 2), ones(24, 1), zeros(24, 2);...
                    zeros(23, 2), -0.3*ones(23, 1), zeros(23, 2);...
                    zeros(23, 2), 0.3*ones(23, 1), zeros(23, 2);...
                    zeros(23, 5);...
                    zeros(48, 3), 0.125*ones(48, 1), zeros(48, 1);...
                    zeros(24, 3), ones(24, 1), zeros(24, 1);...
                    zeros(24, 5);...
                    zeros(48, 4), 0.125*ones(48, 1);...
                    zeros(24, 4), ones(24, 1);...
                    zeros(73, 5)
                    ];

            vec_d = [zeros(96, 1); E_PGUmin*ones(24, 1); E_PGUmax*ones(24, 1); E_CCSmax*ones(24, 1);...
                    -dita_E_PGU*ones(23, 1); dita_E_PGU*ones(23, 1); zeros(498, 1); edemand; zeros(48, 1)
                ];

            % Solve the SUB-PI
            e = ones(size(s_p));
            
            Obj_PI = e'*(s_p + s_n); % SUB-PI objective
            
            options = sdpsettings('verbose', 0, 'solver', 'gurobi', 'gurobi.Method', 0);
            result_PI = optimize(C_sub_PI, Obj_PI, options);
            if result_PI.problem == 0
                if value(Obj_PI) > 1.00e-6
                    % subproblem is infeasible
                    status = 1;
                    pi = dual(C_sub_PI(1:35));
                    if isnan(pi) ~= zeros(size(pi))
                        disp('aka');
                        pi = zeros(size(pi));
                    end
                    if abs(pi'*(mat_B*x_mid' + vec_d) + value(Obj_PI)) < 0.001*abs(value(Obj_PI))
                        pi = -pi;
                    end
                else
                    % subproblem is feasible
                    status = 0;
            % Solve the Subproblem
                    options = sdpsettings('verbose', 0, 'solver', 'gurobi', 'gurobi.Method', 0);
                    result = optimize(C_sub, Obj_sub, options);
                    if result.problem == 0
                        pi = dual(C_sub(1:35));
                        if abs(-pi'*(mat_B*x_mid' + vec_d) - sum(E_Pri.*edemand) - value(Obj_sub)) < 0.001*abs(value(Obj_sub))
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
            %%
            
            if status == 0 % generate optimality cut
                Cut_set = [Cut_set, theta(i) >= pi'*mat_B*x' + pi'*vec_d];
                if abs(pi'*(mat_B*x_mid' + vec_d) - sum(E_Pri.*edemand) - value(Obj_sub)) > 0.001*abs(value(Obj_sub))
                    disp('Optimality cut is wrong.');
                end
                if trigger == 0
                    z_hat = z_hat + vec_p(i)*value(Obj_sub);
                end
            else % generate feasible cut
                Cut_set = [Cut_set, 0 >= pi'*mat_B*x' + pi'*vec_d];
                if abs(pi'*(mat_B*x_mid' + vec_d) - value(Obj_PI)) > 0.001*abs(value(Obj_PI))
                    disp('Feasible cut is wrong.');
                end
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

    disp(upp_z - low_z);
    % disp(low_z);
    disp(x_star);
    % iteration stopping criteria
    if abs(upp_z - low_z) <= eps*max(abs(upp_z), abs(low_z)) || Ite == IteMax
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


%% Solve all solutions
E_PV_deter = [];
E_PGU_deter = [];
E_CCS_deter = [];
E_EL_deter = [];
E_CH4_deter = [];

V_H2_EL_deter = [];
V_H2_out_deter = [];
V_H2_P2G_deter = [];
V_H2_in_deter = [];
V_H2_deter = [];

V_CO2_CCS_deter = [];
V_CO2_out_deter = [];
V_CO2_P2G_deter = [];
V_CO2_in_deter = [];
V_CO2_deter = [];


for i = 1:scr_n
            it = It(i, :)'; edemand = Edemand(i, :)';
            A_PV_v = value(A_PV); E_EL_S_v = value(E_EL_S); E_CH4_S_v = value(E_CH4_S);
            V_H2_S_v = value(V_H2_S); V_CO2_S_v = value(V_CO2_S);

            % objective function of subproblem
            C_tari = sum(350*E_PGU*320/10000-E_Pri.*edemand- 500*V_CH4_P2G);
            
            C_tax = sum((V_CO2_PGU-V_CO2_CCS)*0.6);
            
            Obj_sub = C_tari + C_tax;

            % subproblem constraints
            C_sub = [];

            C_sub = [C_sub, E_PV + E_PV_cur == ita_PV*it*A_PV_v];

            C_sub = [C_sub, e_PGU*E_PGU - V_CO2_PGU == 0,...
                V_CO2_CCS + V_CO2_cur - ita_CCS_max*V_CO2_PGU == 0,...
                E_CCS - lamdaCO2*V_CO2_CCS == 0,...
                E_PGU - s_E_PGU_1 == E_PGUmin,...
                E_PGU + s_E_PGU_2 == E_PGUmax,...
                E_CCS + s_E_CCS == E_CCSmax,...
                E_PGU(2:24)-E_PGU(1:23) - s_E_PGU_3 == -dita_E_PGU,...
                E_PGU(2:24)-E_PGU(1:23) + s_E_PGU_4 == dita_E_PGU
                ];

            C_sub = [C_sub, E_EL - lamda_H2*V_H2_EL == 0,...
                E_CH4 - lamda_CH4*V_CH4_P2G == 0,...
                V_CH4_P2G - V_CO2_P2G == 0,...
                V_CH4_P2G - V_H2_P2G/4 == 0,...
                E_EL - s_E_EL_1 == 0.3*E_EL_S_v,...
                E_EL + s_E_EL_2 == E_EL_S_v,...
                E_EL(2:24)-E_EL(1:23)-s_E_EL_3 == -0.4*E_EL_S_v,...
                E_EL(2:24)-E_EL(1:23)+s_E_EL_4 == 0.2*E_EL_S_v,...
                E_EL(1:23)-E_EL(2:24)-s_E_EL_5 == -0.4*E_EL_S_v,...
                E_EL(1:23)-E_EL(2:24)+s_E_EL_6 == 0.2*E_EL_S_v,...
                E_CH4 + s_E_CH4_1 == E_CH4_S_v,...
                E_CH4(1:23)-E_CH4(2:24)-s_E_CH4_2 == -0.3*E_CH4_S_v,...
                E_CH4(1:23)-E_CH4(2:24)+s_E_CH4_3 == 0.3*E_CH4_S_v
                ];

            C_sub = [C_sub, V_H2(2:24)-V_H2(1:23)-V_H2_in(1:23)+V_H2_out(1:23) == 0,...
                V_H2_in+s_H2_in == 0.125*V_H2_S_v,...
                V_H2_out+s_H2_out == 0.125*V_H2_S_v,...
                V_H2+s_H2 == V_H2_S_v,...
                V_H2(1) - V_H2(24) == 0
                ];

            C_sub = [C_sub, V_CO2(2:24)-V_CO2(1:23)-V_CO2_in(1:23)+V_CO2_out(1:23)==0,...
                V_CO2_in+s_CO2_in == 0.125*V_CO2_S_v,...
                V_CO2_out+s_CO2_out == 0.125*V_CO2_S_v,...
                V_CO2+s_CO2 == V_CO2_S_v,...
                V_CO2(1)-V_CO2(24)==0
                ];

            C_sub = [C_sub, E_PV + E_PGU - E_CCS - E_EL - E_CH4 == edemand];

            C_sub = [C_sub, V_H2_EL+V_H2_out-V_H2_P2G-V_H2_in == 0];

            C_sub = [C_sub, V_CO2_CCS+V_CO2_out-V_CO2_P2G-V_CO2_in == 0];

            C_sub = [C_sub, C_sub_PoADD];

            options = sdpsettings('verbose', 0, 'solver', 'gurobi', 'gurobi.Method', 0);
            result_deter = optimize(C_sub, Obj_sub, options);

            E_PV_deter = [E_PV_deter; value(E_PV)];
            E_PGU_deter = [E_PGU_deter; value(E_PGU)];
            E_CCS_deter = [E_CCS_deter; value(E_CCS)];
            E_EL_deter = [E_EL_deter; value(E_EL)];
            E_CH4_deter = [E_CH4_deter; value(E_CH4)];

            V_H2_EL_deter = [V_H2_EL_deter; value(V_H2_EL)];
            V_H2_out_deter = [V_H2_out_deter; value(V_H2_out)];
            V_H2_P2G_deter = [V_H2_P2G_deter; value(V_H2_P2G)];
            V_H2_in_deter = [V_H2_in_deter; value(V_H2_in)];
            V_H2_deter = [V_H2_deter; value(V_H2)];

            V_CO2_CCS_deter = [V_CO2_CCS_deter; value(V_CO2_CCS)];
            V_CO2_out_deter = [V_CO2_out_deter; value(V_CO2_out)];
            V_CO2_P2G_deter = [V_CO2_P2G_deter; value(V_CO2_P2G)];
            V_CO2_in_deter = [V_CO2_in_deter; value(V_CO2_in)];
            V_CO2_deter = [V_CO2_deter; value(V_CO2)];
end




%% Plot
% Power balance
figure  
Plot_E_yuan=zeros(2,144);       
  
Plot_E_yuan(1, :)=E_PV_deter'/1000;      
Plot_E_yuan(2, :)=E_PGU_deter'/1000;

Plot_E_fuhe =zeros(4,144); 

Plot_E_fuhe(1, :) = -E_CCS_deter'/1000;
Plot_E_fuhe(2, :) = -E_EL_deter'/1000;
Plot_E_fuhe(3, :) = -E_CH4_deter'/1000;
Plot_E_fuhe(4, :) = -[Edemand(1, :), Edemand(2, :), Edemand(3, :), Edemand(4, :), Edemand(5, :), Edemand(6, :)]/1000;

bar(Plot_E_yuan','stacked');
hold on
bar(Plot_E_fuhe','stacked');
xlabel('Time/h');
ylabel('Power/MW');
title('Electric Power Balance');
legend('PV','DG', ...
       'Carbon Device','Hydrogen Device','CH4 Device','Load'); 
xlim([0 150   ]);
%set(gca,'FontSize',12 );         

%% H2 balance
    
figure  
% yyaxis left
Plot_H_yuan=zeros(2,144);       
       
Plot_H_yuan(1, :)=V_H2_EL_deter';      
Plot_H_yuan(2, :)=V_H2_out_deter';

Plot_H_fuhe =zeros(2,144); 

Plot_H_fuhe(1, :) = -V_H2_P2G_deter';
Plot_H_fuhe(2, :) = -V_H2_in_deter';

bar(Plot_H_yuan','stacked');
hold on
bar(Plot_H_fuhe','stacked');
xlabel('Time/h');
ylabel('Power/N.m3');
% yyaxis right
% hold on
plot(V_H2_deter,'k','LineWidth',2);
ylabel('H2 Volume/N.m3');
title('H2 Balance');
legend('H2 Generation','H2 Out','H2 for P2G','H2 In','H2 Volume/N.m3'); 
xlim([0 150   ]);
%set(gca,'FontSize',12 );  

%% CO2 balance
  
figure  
Plot_C_yuan=zeros(2,144);       
     
Plot_C_yuan(1, :)=V_CO2_CCS_deter';      
Plot_C_yuan(2, :)=V_CO2_out_deter';

Plot_C_fuhe =zeros(2,144); 

Plot_C_fuhe(1, :) = -V_CO2_P2G_deter';
Plot_C_fuhe(2, :) = -V_CO2_in_deter';

bar(Plot_C_yuan','stacked');
hold on
bar(Plot_C_fuhe','stacked');
hold on
plot(V_CO2_deter,'-k','LineWidth',2)
xlabel('Time/h');
ylabel('Power/N.m3');
title('CO2 Balance');
legend('CCS','Out','P2G','In','CO2 Vloume/N.m3'); 
xlim([0 150   ]);
%set(gca,'FontSize',12 );  

%%  容量配置结果

fprintf('配置光伏容量%dMW\n',k*A_PV_v/1000 );
fprintf('配置电制氢容量%d/MW\n',E_EL_S_v/1000 );
fprintf('配置甲烷化容量%d/MW\n',E_CH4_S_v/1000 );
fprintf('配置氢气存储容量%d/N.m3\n',V_H2_S_v );
fprintf('配置CO2存储容量%d/N.m3\n',V_CO2_S_v );
