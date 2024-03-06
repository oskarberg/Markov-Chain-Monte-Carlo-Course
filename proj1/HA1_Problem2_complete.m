%% Calc Probability of (V in I=(a,b))
load powercurve_D240.mat


% Parameters for Weibull distribution
lambda =  [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10.0 10.9 11.7 11.7];
k =       [2.0  2.0  2.0  1.9 1.9 1.9 1.9 1.9 2.0  1.9  2.0  2.0];
a = 4; % Lower bound
b = 25; % Upper bound

PinI_save = zeros(12,1);

for i=1:12

    total_probability = wblcdf(b, lambda(i), k(i)) - wblcdf(a, lambda(i), k(i));
    PinI_save(i) = total_probability;
end

PinI_save

%% Choosing the instrumental function

k = 1.95;
lambda = 10.05;
v = 4:0.1:25;
P_v = P(v); 
WeibullPDF_v = wblpdf(v, lambda, k); 
GammaPDF_v = gampdf(v, 10, 11/9); 


scaled_P_v = P_v / max(P_v);


expression = scaled_P_v' .* WeibullPDF_v; %./ GammaPDF_v; 

plot(v, expression)
hold on
plot(v, GammaPDF_v)



%% a) Crude monte carlo and b) covariate method



lambda =  [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10.0 10.9 11.7 11.7];
k =       [2.0  2.0  2.0  1.9 1.9 1.9 1.9 1.9 2.0  1.9  2.0  2.0];

N = 10000;

Monthly_expected_power = zeros(12,1);
Monthly_std = zeros(12,1);
Monthly_conf_int = zeros(12,2);

Monthly_expected_Z = zeros(12,1);
Monthly_std_Z = zeros(12,1);
Monthly_conf_int_Z = zeros(12,2);

for i = 1: 12

    
    V = wblrnd(lambda(i), k(i), N, 1);
    
    Power = P(V);
    mean_power = mean(Power);


    % ------------------START Control Variate--------------%
    control_variate = -V;
    m = -gamma(1 + 1/k(i)) * (lambda(i)^1);

    cov_covariate = mean((Power - mean_power) .* (-V - m));
    var_covariate = mean((-V - m).^2);


    beta = -cov_covariate/var_covariate;
    % New process
    Z = Power + beta*(-V - m);

    Monthly_std_Z(i) = std(Z);
    Monthly_expected_Z(i) = mean(Z);

    Monthly_conf_int_Z(i,1) = mean(Z) + norminv(0.005) * std(Z)/sqrt(N);
    Monthly_conf_int_Z(i,2) = mean(Z) + norminv(0.995) * std(Z)/sqrt(N);
    % ----------------END Control Variate------------------%


    Monthly_expected_power(i) = mean_power;

    std_power = std(Power);
    Monthly_std(i) = std_power;
   
    Monthly_conf_int(i,1) = mean_power + norminv(0.005) * std_power/sqrt(N);
    Monthly_conf_int(i,2) = mean_power + norminv(0.995) * std_power/sqrt(N);
end



Monthly_expected_power
%Monthly_std
Monthly_conf_int
Monthly_conf_int(:,2) - Monthly_conf_int(:,1)


%% b) Print values for covariate method

Monthly_expected_Z
%Monthly_std_Z
Monthly_conf_int_Z
Monthly_conf_int_Z(:,2) - Monthly_conf_int_Z(:,1)


%% a) Truncated Weibull sampling using inverse method


a = 4;
b = 25;

N = 10000; 

lambda = [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10.0 10.9 11.7 11.7];
k = [2.0  2.0  2.0  1.9 1.9 1.9 1.9 1.9 2.0  1.9  2.0  2.0];

Monthly_expected_power = zeros(12,1);
Monthly_std = zeros(12,1);
Monthly_conf_int = zeros(12,2);


%Monthly_expected_Z = zeros(12,1);
%Monthly_std_Z = zeros(12,1);

for i = 1:12

    F_a = wblcdf(a, lambda(i), k(i));
    F_b = wblcdf(b, lambda(i), k(i));

    U = rand(N, 1);

    adjusted_U = F_a + U * (F_b - F_a);
    
    V_conditional = wblinv(adjusted_U, lambda(i), k(i));


    % pdf_weibull = wblpdf(V_conditional, lambda(i), k(i));
    % pdf_importance = wblpdf(V_conditional, lambda(i), k(i)) / (F_b - F_a);
    % Power = P(V_conditional).* (pdf_weibull./ pdf_importance);

    Power = P(V_conditional);
    
    mean_power = mean(Power);


    %----------------- START Control Variate--------------
    % control_variate = -V_conditional;
    % m = -m_save(i);
    % 
    % cov_covariate = mean((Power - mean_power) .* (-V_conditional - m));
    % var_covariate = mean((-V_conditional - m).^2);
    % 
    % beta = -cov_covariate/var_covariate;
    % % New process
    % Z = Power + beta*(-V_conditional - m);
    % 
    % Monthly_expected_Z(i) = mean(Z);
    % Monthly_std_Z(i) = std(Z);
     %---------------- END Control Variate--------------------

 
    std_power = std(Power);
    
    Monthly_expected_power(i) = mean_power;
    Monthly_std(i) = std_power;
    
    Monthly_conf_int(i,1) = mean_power + norminv(0.005) * std_power / sqrt(N);
    Monthly_conf_int(i,2) = mean_power + norminv(0.995) * std_power / sqrt(N);
end

Monthly_expected_power .* PinI_save
%Monthly_std
%Monthly_conf_int
Monthly_conf_int .* PinI_save
Monthly_conf_int(:,2) - Monthly_conf_int(:,1)



%% c) Importance sampling


const1 = 10;
const2 = 11/9; 

a = 4; 
b = 25; 
N = 10000; 


Monthly_expected_power = zeros(12, 1);
Monthly_std = zeros(12, 1);

Monthly_conf_int = zeros(12,2);

for i = 1:12


    U = rand(N, 1); 
    F_a_g = gamcdf(a, const1, const2); 
    F_b_g = gamcdf(b, const1, const2);  
    adjusted_U = F_a_g + U * (F_b_g - F_a_g);
    V_g = gaminv(adjusted_U, const1, const2);

    % Calculating importance weights
    pdf_weibull = wblpdf(V_g, lambda(i), k(i));
    pdf_importance = gampdf(V_g, const1, const2) / (F_b_g - F_a_g);
    weights = pdf_weibull ./ pdf_importance;

    Power_g = P(V_g);
    weighted_power = Power_g .* weights;

    weighted_mean_power = mean(weighted_power);
    weighted_std_power = std(weighted_power);

    % Store results
    Monthly_expected_power(i) = weighted_mean_power;
    Monthly_std(i) = weighted_std_power;

    Monthly_conf_int(i,1) = weighted_mean_power + norminv(0.005) * weighted_std_power / sqrt(N);
    Monthly_conf_int(i,2) = weighted_mean_power + norminv(0.995) * weighted_std_power / sqrt(N);
end

Monthly_expected_power
Monthly_conf_int
Monthly_conf_int(:,2) - Monthly_conf_int(:,1)



%% d) Antithectic Varaible

a = 4;
b = 25;

N = 10000; 

lambda = [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10.0 10.9 11.7 11.7];
k =      [2.0  2.0  2.0  1.9 1.9 1.9 1.9 1.9 2.0  1.9  2.0  2.0];

Monthly_expected_power = zeros(12,1);
Monthly_std = zeros(12,1);
Monthly_conf_int = zeros(12,2);

for i = 1:12
    
    F_a = wblcdf(a, lambda(i), k(i));
    F_b = wblcdf(b, lambda(i), k(i));

    U1 = rand(N/2, 1);
    U2 = 1 - U1;

    adjusted_U1 = F_a + U1 * (F_b - F_a);
    adjusted_U2 = F_a + U2 * (F_b - F_a);

    V_conditional1 = wblinv(adjusted_U1, lambda(i), k(i));
    V_conditional2 = wblinv(adjusted_U2, lambda(i), k(i));

    Power1 = P(V_conditional1);
    Power2 = P(V_conditional2);

    Average_Power = (Power1 + Power2) / 2;

    mean_power = mean(Average_Power);
    std_power = std(Average_Power);

    Monthly_expected_power(i) = mean_power;
    Monthly_std(i) = std_power;

    Monthly_conf_int(i,1) = mean_power + norminv(0.005) * std_power / sqrt(N/2);
    Monthly_conf_int(i,2) = mean_power + norminv(0.995) * std_power / sqrt(N/2);
end


Monthly_expected_power .* PinI_save
Monthly_conf_int.*PinI_save
%Monthly_std

Monthly_conf_int(:,2) - Monthly_conf_int(:,1)



%% e)

Probability_value = PinI_save % PinI_save from previous script

%% f) Dont forget to load the importance sampling section before


lambda = [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10.0 10.9 11.7 11.7];
k =      [2.0  2.0  2.0  1.9 1.9 1.9 1.9 1.9 2.0  1.9  2.0  2.0];

rho = 1.225; 
d = 240;  

E_Ptot = zeros(1, 12);

for i = 1:12
    E_V_cubed = gamma(1 + 3/k(i)) * lambda(i)^3;
    E_Ptot(i) = 0.5 * rho * pi * (d/2)^2 * E_V_cubed;
end

disp(E_Ptot);


%Monthyl expected power from previous best. Importance sampling the best.
Avg_Power_Ratio = Monthly_expected_power' ./ E_Ptot;

Avg_Power_Ratio= Avg_Power_Ratio' %per month

std_Avg_Power_Ratio = (1./E_Ptot)' .* Monthly_std;

Monthly_conf_int_Avg_Ratio = zeros(12,2);
Monthly_conf_int_Avg_Ratio(:,1) = Avg_Power_Ratio + norminv(0.005) * std_Avg_Power_Ratio/sqrt(N);
Monthly_conf_int_Avg_Ratio(:,2) = Avg_Power_Ratio + norminv(0.995) * std_Avg_Power_Ratio/sqrt(N);

Monthly_conf_int_Avg_Ratio
Monthly_conf_int_Avg_Ratio(:,2) - Monthly_conf_int_Avg_Ratio(:,1)


%% g)

hours_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]' * 24; % Hours in each month
rated_power = 15000000; % Rated power in W

% Capacity Factor Calculation
Max_Possible_Output = rated_power*hours_per_month; % Max possible output per month

Availability_Factor = PinI_save; % PinI_save from previous script

Real_output = Monthly_expected_power.* hours_per_month.* Availability_Factor;

Capacity_Factor = Real_output ./ Max_Possible_Output


mean(Capacity_Factor)
mean(Availability_Factor)



