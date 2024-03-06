%% 3a)

load('powercurve_D240.mat')
clc

k = 1.95;
lambda = 10.05;



const1 = 10;
const2 = 11/9;

a = 4; 
b = 25; 
N = 10000; 



% Truncated Gamma sampling
U = rand(N, 1);  % Reusing U for Gamma sampling
F_a_g = gamcdf(a, const1, const2); 
F_b_g = gamcdf(b, const1, const2);  
adjusted_U = F_a_g + U * (F_b_g - F_a_g);
V_g = gaminv(adjusted_U, const1, const2);

% Calculate the importance weights
pdf_weibull = wblpdf(V_g, lambda, k);
pdf_importance = gampdf(V_g, const1, const2) / (F_b_g - F_a_g);
weights = pdf_weibull ./ pdf_importance;



% Calculate weighted power
Power_g = P(V_g); % P is the power curve function
weighted_power = Power_g .* weights;

% Calculate mean and standard deviation
weighted_mean_power = mean(weighted_power)
weighted_std_power = std(weighted_power)


Monthly_conf_int_LB = weighted_mean_power + norminv(0.005) * weighted_std_power/sqrt(N)
Monthly_conf_int_UB = weighted_mean_power + norminv(0.995) * weighted_std_power/sqrt(N)


Epv1 = weighted_mean_power;

%% Crude monte carlo for 3c), estiamte for marginal variances 
V = wblrnd(lambda, k, N, 1);

Power = P(V);
mean_power_crude = mean(Power)
Marginal_var_power_crude = var(Power)




%% Motivation for choice of importance function (Plotting of 3d fucntions using surf)

% Define the parameters
k = 1.95;
lambda = 10.05;
alpha = 0.638;
p = 3;
q = 1.5;


% definitions for the univariate Weibull PDF and CDF functions 
f = @(v) wblpdf(v, lambda, k); % PDF
F = @(v) wblcdf(v, lambda, k); % CDF


% Create a meshgrid for v1 and v2
[v1, v2] = meshgrid(0:0.5:35, 0:0.5:35); % Adjust range as necessary

% Correcting the Compute the joint PDF code
jointPDF = f(v1) .* f(v2) .* (  1 + alpha * (1 - F(v1).^p).^(q-1) .* (1 - F(v2).^p).^(q-1) .* ...
             ( F(v1).^p * (1 + p*q) - 1 ).* ( F(v2).^p* (1 + p*q) - 1 )  );



%total_power = P(v1).*P(v2);
total_power = P(v1) + P(v2);
total_power = reshape(total_power, size(v1,1), size(v1,1) );


Indicator_power = total_power > 15000000;
%Indicator_power(Indicator_power == 0) = NaN;

figure(1)

% Plot the joint PDF using surf
surf(v1, v2, jointPDF);
xlabel('v1');
ylabel('v2');
zlabel('Joint PDF');
title('Bivariate Weibull Distribution');
xlim([0 30]);
ylim([0 30]);


% Indicator_power = total_power < 15000000;
% 
% figure(3)
% 
% % Plot the joint PDF using surf
% surf(v1, v2, Indicator_power.*jointPDF);
% xlabel('v1');
% ylabel('v2');
% zlabel('Joint PDF');
% title('Bivariate Weibull Distribution');
% xlim([0 35]);
% ylim([0 35]);

% Indicator_power = total_power == 15000000;
% 
% figure(4)
% 
% % Plot the joint PDF using surf
% surf(v1, v2, Indicator_power.*jointPDF);
% xlabel('v1');
% ylabel('v2');
% zlabel('Joint PDF');
% title('Bivariate Weibull Distribution');
% xlim([0 35]);
% ylim([0 35]);



figure(2)

%k1 = 6; theta1 = 4/5;
%k1 = 7; theta1 =9/6;


%g = @(v1, v2) gampdf(v1, k1, theta1) .* gampdf(v2, k1, theta1);
% Define the bivariate independent weibull PDF
g = @(v1,v2) wblpdf(v1, lambda, k) .* wblpdf(v2, lambda, k) 

pdfValues = g(v1, v2);

surf(v1, v2, pdfValues);
xlabel('v1');
ylabel('v2');
zlabel('Weibull PDF');
title('Bivariate Independent Weibull PDF');
xlim([0 30]);
ylim([0 30]);


%% Importance sampling 3b)

clc
const1 = 10;
const2 = 11/9; 

a = 4; 
b = 25; 
N = 10000; 



% Truncated Gamma sampling, Sample from 
U1 = rand(N, 1);  
U2 = rand(N, 1);
F_a_g = gamcdf(a, const1, const2); 
F_b_g = gamcdf(b, const1, const2);  
adjusted_U1 = F_a_g + U1 * (F_b_g - F_a_g);
adjusted_U2 = F_a_g + U2 * (F_b_g - F_a_g);
V1_g = gaminv(adjusted_U1, const1, const2);
V2_g = gaminv(adjusted_U2, const1, const2);


% joint bivariate Weibull PDF
jointPDF_Vg = f(V1_g) .* f(V2_g) .* (1 + alpha * (1 - F(V1_g).^p).^(q-1) .* (1 - F(V2_g).^p).^(q-1) .* ...
             ( F(V1_g).^p * (1 + p*q) - 1 ).* ( F(V2_g).^p* (1 + p*q) - 1 ));


% Calculate the PDF of the instrumental gamma distribution for each sample
instrumentalPDF_V1_g = gampdf(V1_g, const1, const2);
instrumentalPDF_V2_g = gampdf(V2_g, const1, const2);

% Adjust for truncation
%F_b_g = gamcdf(b, const1, const2);  
%F_a_g = gamcdf(a, const1, const2);  
adjustment_factor = F_b_g - F_a_g;

% Compute the importance weights
weights = jointPDF_Vg ./ (instrumentalPDF_V1_g .* instrumentalPDF_V2_g / adjustment_factor^2);


weighted_power = P(V1_g).*P(V2_g) .* weights;

% Calculate mean and standard deviation
weighted_mean_power = mean(weighted_power)
weighted_std_power = std(weighted_power)/sqrt(N)

Epv1pv2 = weighted_mean_power;

%% Some calculations for tasks 3b and 3c

Covariance = Epv1pv2 - Epv1^2

Variance = 2*Marginal_var_power_crude + 2*Covariance


%% Importance sampling 3d)
clc


a = 0; % Untruncated 
b = inf; % Untruncated 
N = 10000; 


% Truncated Gamma sampling, Sample from 
U1 = rand(N, 1);  
U2 = rand(N, 1);
F_a_g = wblcdf(a, const1, const2); 
F_b_g = wblcdf(b, const1, const2);  
adjusted_U1 = F_a_g + U1 * (F_b_g - F_a_g);
adjusted_U2 = F_a_g + U2 * (F_b_g - F_a_g);
V1_g = wblinv(adjusted_U1, lambda, k);
V2_g = wblinv(adjusted_U2, lambda, k);



% joint bivariate Weibull PDF
jointPDF_Vg = f(V1_g) .* f(V2_g) .* (1 + alpha * (1 - F(V1_g).^p).^(q-1) .* (1 - F(V2_g).^p).^(q-1) .* ...
             ( F(V1_g).^p * (1 + p*q) - 1 ).* ( F(V2_g).^p* (1 + p*q) - 1 ));


% Calculate the PDF of the instrumental gamma distribution for each sample
adjustment_factor = F_b_g - F_a_g; % = 1
instrumentalPDF_V1_g = wblpdf(V1_g, lambda, k) / adjustment_factor;
instrumentalPDF_V2_g = wblpdf(V2_g, lambda, k) / adjustment_factor;


% Compute the importance weights
weights = jointPDF_Vg ./ (instrumentalPDF_V1_g .* instrumentalPDF_V2_g );

condition_met =  (P(V1_g) + P(V2_g)) > 15000000; 
condition_not_met =  (P(V1_g) + P(V2_g)) < 15000000;

% Apply weights to the condition (indicator function)
weighted_condition = condition_met .* weights;
weighted_condition_not = condition_not_met .* weights;

% Estimate the probability
probability_estimate = sum(weighted_condition) / sum(weights)
probability_estimate_not = sum(weighted_condition_not) / sum(weights)

sum(probability_estimate + probability_estimate_not)

% Confidence Intervals
conf_int_LB = probability_estimate + norminv(0.005) * std(weighted_condition)/sqrt(N)
conf_int_UB = probability_estimate + norminv(0.995) * std(weighted_condition)/sqrt(N)

conf_int_not_LB = probability_estimate_not + norminv(0.005) * std(weighted_condition_not)/sqrt(N)
conf_int__not_UB = probability_estimate_not + norminv(0.995) * std(weighted_condition_not)/sqrt(N)

%std(weighted_condition)/sqrt(N)