global m1 m2 l1 l2 g 

m1 = 0.15; m2 = 0.36;l1 = 0.3;l2 = 0.2; g  = 9.81;

dist_target = 0.8; 
v_target = 0.2; 
Ts = 0.01;

[t_vec, pos_prof, vel_prof] = generate_profile_IS(dist_target, v_target, 1, l1, l2, m1, m2, Ts);
% หรือ [t_vec, pos_prof, vel_prof] = generate_profile_ACCRT(dist_target, v_target, 1, l1, l2, m1, m2, Ts);

accel_prof = gradient(vel_prof) ./ Ts;

ode_func = @(t, Y) crane_dynamics(t, Y, t_vec, accel_prof, m1, m2, l1, l2, g);
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
[t_ode, Y] = ode45(ode_func, [0 25], [0; 0; 0; 0], options);

th1_sim = rad2deg(Y(:,1));
th2_sim = rad2deg(Y(:,2));

% MATHEMATICAL MODELING
function dYdt = crane_dynamics(t, Y, time_vec, accel_profile, m1, m2, l1, l2, g)
    a_trolley = interp1(time_vec, accel_profile, t, 'linear', 0);
    M = [(m1+m2)*l1^2,   m2*l1*l2; 
         m2*l1*l2,       m2*l2^2];
    K = [(m1+m2)*g*l1,   0; 
         0,              m2*g*l2];
    B = [-(m1+m2)*l1; 
         -m2*l2];

    theta = [Y(1); Y(2)]; 
    dtheta = [Y(3); Y(4)];

    ddtheta = M \ (B * a_trolley - K * theta);
    
    dYdt = [dtheta; ddtheta];
end