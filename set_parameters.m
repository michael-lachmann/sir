% function set_parameters()
T_Y_TO_R_PARA = 22.6; % // median of  [21.1, 22.6, 24.4]
T_EXPOSED_PARA = 7.0; % // median of  [5.6, 7, 8.2]
ASYMP_RATE = 0.179;
PROP_TRANS_IN_E = 0.126;
T_ONSET_TO_H = 5.9;
T_H_TO_R = 14.0;
T_H_TO_D = 14.0;
symp_h_ratio_overall = [0.00048721, 0.00048721, 0.03287572, 0.11337395, 0.17733063] ;
hosp_f_ratio = [0.04 , 0.12365475, 0.03122403, 0.10744644, 0.23157691] ;
symp_h_ratio = [2.79135866e-04, 2.14621858e-04, 1.32154040e-02, 2.85633688e-02, 3.38733218e-02, ...
    2.79135866e-03, 2.14621858e-03, 1.32154040e-01, 2.85633688e-01, 3.38733218e-01  ] ; % // for now only first row is used because there is only one risk group.
symp_h_ratio_corrrect = [ 1.610326595443765, 1.924464960134284, 2.31133016137442, 3.724051596082457, 4.95257504190157]; % // hospitalization is 10 times higher for high risk individuals. I made an average of 1*[num low]+10*[num high] for each age class.

g.gamma_h = 1.0 / T_H_TO_R ;
gamma_y_c = 1.0 / T_Y_TO_R_PARA ;
g.gamma_y = gamma_y_c ;
g.gamma_a = g.gamma_y ;

sigma_c = 1.0 / T_EXPOSED_PARA ;
g.sigma = sigma_c ;

g.eta = 1.0 / T_ONSET_TO_H ;

g.mu = 1.0 / T_H_TO_D ;

g.tau = 1.0 - ASYMP_RATE ;

g.omega_y = 1.0 ;
g.omega_h = 0.0 ;

for i=1:g.ngroup
    symp_h_ratio(i) = symp_h_ratio(i)*symp_h_ratio_corrrect(i) ;
    
    for i=1:g.ngroup
        g.nu(i) = hosp_f_ratio(i) * g.gamma_h / (g.mu + (g.gamma_h - g.mu  ) * hosp_f_ratio(i) ) ; % // hosp_f_ratio is an array of size age.
        g.pi(i) = symp_h_ratio(i) * g.gamma_y / (g.eta + (g.gamma_y - g.eta) * symp_h_ratio(i) ) ; % // symp_h_ratio is an array
        % // omega_e - relative infectiousness in E, IY, IA % // symp_h_ratio_overall length age.
        g.omega_e(i) = ((symp_h_ratio_overall(i) / g.eta) + ((1 - symp_h_ratio_overall(i) ) / g.gamma_y)) * ...
            g.omega_y * g.sigma * PROP_TRANS_IN_E / (1 - PROP_TRANS_IN_E) ;
        g.omega_a(i) = ((symp_h_ratio_overall(i) / g.eta) + ((1 - symp_h_ratio_overall(i)) / gamma_y_c)) * ...
            g.omega_y * sigma_c * PROP_TRANS_IN_E / (1 - PROP_TRANS_IN_E) ;
    end
end

SN = {'S', 'E', 'A', 'Y', 'H', 'R', 'D'};




