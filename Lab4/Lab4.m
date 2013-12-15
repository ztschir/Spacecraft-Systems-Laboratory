%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Lab 4 Top Level Script%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Lab4

    clear all;
    clc;

%%% Constants %%% 
    global rho0;
    global beta;
    global g_mars;
    global r_mars;
    global H;
    global q_dot_int;
    
    rho0         = 0.0525;           % kg/m^3
    g_mars       = 3.69;             % m/s^2
    r_mars       = 3397000;          % m
    g_earth      = 9.81;             % m/s^2
    R            = 8314.472;         % J/K-kmol 
    molar_weight = 43.27;            % kg/kmol
    R_bar        = R / molar_weight; % J/K-kg
    T_init       = 150;              % K
    H            = (R_bar * T_init) / g_mars

    mass         = 603;              % kg
    c_d          = 1.65; 
    diameter     = 2.55;             % m
    ref_area     = pi*(diameter/2)^2; 

    bolt_const   = 5.67E-8;          % W/m^2-K^4
    graves_const = 1.9027E-4; 

%%% Problem 1 %%%

%%% Part a %%%
    beta      = mass / (c_d * ref_area)


%%% Part b %%%
    gamma = 13;         % deg
    C = - (rho0 * H) / (2 * beta * sind(gamma))
               

%%% Part c %%%

    height = 125000:-1:0;  % m
    vel_atm = 5.6 * 1000; % m/s
    vel = vel_atm * exp(C * exp(-height/H));

    %plot(vel, height), grid on;
    %xlabel('Velocity (m/s)'), ylabel('Altitude (m)')


%%% Part d %%%

    max_deacc_height = H * log(-2 * single(C))
    max_magnitude    = (vel_atm^2 * sind(gamma)) / (53.3 * H)


%%% Part e %%%
    tspan = 1:10000;
    v0 = vel_atm;
    h0 = max(height);
    s0 = [v0; h0; gamma];
    
    options = odeset('RelTol',1e-8, 'AbsTol', 1e-8, 'Events', @eventfun);
    [t,s]   = ode45(@odefun,tspan,s0,options);
    %plot(vel, height, s(:,1), s(:,2)), grid on
    %xlabel('Velocity (m/s)'), ylabel('Altitude (m)')
   
   
   
%%% Part f %%%
    k_mars = 1.9027E-4; 
    r_n    = 0.665; % m
    
    q_dot_allen = k_mars .* (((rho0 .* exp(- (height ./ H))) ./ r_n).^0.5) ...
                  .* vel.^3;
    q_dot_int   = k_mars .* (((rho0 .* exp(- (s(:,2) ./ H))) ./ r_n).^0.5) ...
                  .* s(:,1).^3;
    
    %plot(q_dot_allen, height, q_dot_int, s(:,2)), grid on
    %xlabel('Heatrate (W/m^2)'), ylabel('Altitude (m)')
    
    
%%% Part g %%%
    peak_analytic_q_dot = k_mars * ((beta * sind(gamma)) / (3 * r_n * H))^0.5 ...
                        * ((vel_atm^3)/(exp(0.5)))
    peak_allen_q_dot    = max(q_dot_allen)
    peak_int_q_dot      = max(q_dot_int)
    
    J_s_int             = sum(q_dot_int)
    J_s_allen           = sum(q_dot_allen)
    
    J_s_analytic        = k_mars * vel_atm^2 * ...
                          ((pi * H * beta) / r_n * sind(gamma))^0.5
    
    
%%% Part h %%%
    global emissivity;
    global boltzmann_const;
    emissivity = 0.7; 
    boltzmann_const = 5.670373E-8;
    T_surf = nthroot(q_dot_int ./ (emissivity .* boltzmann_const), 4);
    T_surf_cel = T_surf - 273;     % C
    %plot(t,T_surf_cel), grid on
    %xlabel('Time (s)'), ylabel('Temperture (deg C)')
    
%%% Part i %%%
    int_J_s       = J_s_int / 10000;             % J/cm^2
    TPS_mass_frac = 0.091 * int_J_s^0.51575
    TPS_mass      = (TPS_mass_frac / 100) * mass % kg
    TPS_struc     = TPS_mass * 0.15              % kg
    TPS_tot       = TPS_mass + TPS_struc         % kg
    


%%% Part j %%%
    int_peak_heat_rate   = peak_int_q_dot / 10000   % W/cm^2
    allen_peak_heat_rate = peak_allen_q_dot / 10000 % W/cm^2
    int_J_s                                         % J/cm^2
    allen_J_s            = J_s_allen / 10000        % J/cm^2
    analytic_J_s         = J_s_analytic / 10000     % J/cm^2


%%% Extra Credit Problem 1 %%%
%     global rho_sc;
%     global cp_sc;
%     global k_sc;
%     
%     T0              = 25 + 273;               % K
%     max_bond_surf_T = 200 + 273;              % K
%     rho_sc          = 0.288 * (100^3 / 1000); % kg/m^3
%     cp_sc           = 1120;                   % J/kg-K
%     k_sc            = 0.0576;                 % W/m-K
% 
%     tspan = linspace(0,100,100);
%     xmesh = linspace(0,1,100);
%     m     = 0;
% 
%     sol = pdepe(m,@pdefun,@icfun,@bcfun,xmesh,tspan);
% 
%     surf(xmesh,tspan,sol)
    
    
%%% Problem 2 %%%


%%% Part a %%%
    cd_p      = 0.62;
    V_t       = 57;                              % m/s
    T_surf    = 200;                             % K
    H_spec    = (R_bar * T_surf) / g_mars;       % m
    rho0_surf = 0.0159                           % kg/m^3   
    rho_spec  = rho0_surf * exp(-(1 / H_spec));  % kg/m^3
    diameter  = sqrt((8 * mass * g_mars) / ...
               (pi * cd_p * rho_spec * V_t^2))   % m
    
%%% Part b %%%
    diff_percent = (diameter - 11.7) / 11.7

    

%%% Problem 3 %%%

    mass_lander = 382;  % kg
    V_t_land    = 2.4;  % m/s
    max_load    = 10;   % g's
    load_frac   = 0.7; 
    peak_stress = 1.62 * 1000; % kPa
    max_strain  = 0.66;
    num_legs    = 3;
    
%%% Part a %%%
    S_f = (V_t_land^2) / ...
          (2 * (g_earth * max_load * load_frac - g_mars))     % m
    S_o = S_f / max_strain                                    % m
    A_t = (max_load * mass_lander * g_earth) / peak_stress;   % m^2
    A_l = A_t / num_legs                                      % m^2
    R_l = sqrt( A_l / pi)
    
    V_i = 1:10;                                               % m/s
    S_f_i = (V_i.^2) ./ ...
          (2 .* (g_earth .* max_load .* load_frac - g_mars)); % m
    S_o_i = S_f_i ./ max_strain; 
    plot(V_i, S_o_i), grid on
    ylabel('S0 (m)'), xlabel('Velocity (m/s)')
    
end   

%% Define the ODEs
function [sdot] = odefun(t,s)
    global rho0;
    global beta;
    global g_mars;
    global r_mars;
    global H;
    
    sdot = zeros(3,1);

    V     = s(1);
    h     = s(2);
    gamma = s(3);
    rho   = (rho0 * exp(-h/H));


    sdot(1) = -((rho*(V^2)) / (2 * beta)) + g_mars * sind(gamma);
    sdot(2) = -V * sind(gamma);

    sdot(3) = -(V * cosd(gamma) / (r_mars + s(2))) + ...
               ((g_mars * cosd(gamma)) / V);
end

%% Define the event function to end at 0 altitude
function [height,isterminal,direction] = eventfun(t,s)
    height     = s(2);
    isterminal = 1;
    direction  = 0;
end


%% Define the pde function terms
function [c,f,s] = pdefun(x,t,T,DTDx)
    global rho_sc;
    global cp_sc;
    global k_sc;    

    c = 1 / (k_sc / (cp_sc * rho_sc));
    f = DTDx;
    s = 0;
end

%% Define the initial condition function
function T0 = icfun(x)
    global q_dot_int;
    global emissivity;
    global boltzmann_const;
    T0 = nthroot(q_dot_int ./ (emissivity .* boltzmann_const), 4);

end
%% Define the boundary condition function
function [pl,ql,pr,qr] = bcfun(xl,Tl,xr,Tr,t)
    global emissivity;
    global boltzmann_const;
    
    pl = t;
    %ql = emissivity * boltzmann_const * Tl^4;
    ql = 0;
    pr = 0;
    qr = 1;
    
end



