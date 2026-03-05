function [time_vec, pos_profile, vel_profile] = generate_profile_IS(dist_target, v_target, selected_shaper, lh, ll, mh, ml, Ts)
    % 1. System Identification
    g = 9.81;
    a = mh * lh * ll;
    b = (mh * g * ll) + (ml * g * ll) + (mh * g * lh) + (ml * g * lh);
    c = (mh * g * g) + (ml * g * g);
    discriminant = sqrt(b * b - 4 * a * c);
    wn1 = sqrt((b - discriminant) / (2 * a));
    wn2 = sqrt((b + discriminant) / (2 * a));
    T1 = 2 * pi / wn1;
    T2 = 2 * pi / wn2;

    % 2. Generate FIR Filter
    N = 9;
    n_filter = 0:N-1;
    sigma = (N-1)/6;
    g_gauss = exp(-0.5*((n_filter-(N-1)/2)/sigma).^2);
    h = g_gauss / sum(g_gauss);

    % 3. Call Shaper Function (ต้องมีไฟล์ shaper_xx.m อยู่ในโฟลเดอร์เดียวกัน)
    switch selected_shaper
        case 1, shaper_vec = shaper_01_zv_zv(T1, T2, Ts);
        case 2, shaper_vec = shaper_02_zv_zvd(T1, T2, Ts);
        case 3, shaper_vec = shaper_03_zv_zvdd(T1, T2, Ts);
        case 4, shaper_vec = shaper_04_zvd_zv(T1, T2, Ts);
        case 5, shaper_vec = shaper_05_zvd_zvd(T1, T2, Ts);
        case 6, shaper_vec = shaper_06_zvd_zvdd(T1, T2, Ts);
        case 7, shaper_vec = shaper_07_zvdd_zv(T1, T2, Ts);
        case 8, shaper_vec = shaper_08_zvdd_zvd(T1, T2, Ts);
        case 9, shaper_vec = shaper_09_zvdd_zvdd(T1, T2, Ts);
        otherwise, error('Invalid Shaper Number! (Use 1-9)');
    end

    % 4. Automatic Pulse Generator
    start_time = 0.0;
    pulse_duration = dist_target / v_target;
    shaper_duration = length(shaper_vec) * Ts;
    filter_duration = length(h) * Ts;
    end_time = start_time + pulse_duration;
    t_final = end_time + shaper_duration + filter_duration + 1.0; 
    
    time_vec = 0:Ts:t_final;
    v_command = zeros(size(time_vec));
    v_command(time_vec >= start_time & time_vec < end_time) = v_target;

    % 5. Convolution & Integration
    v_shaped = conv(v_command, shaper_vec);
    v_shaped = v_shaped(1:length(time_vec)); 
    v_filtered = conv(v_shaped, h);
    vel_profile = v_filtered(1:length(time_vec)); 
    pos_profile = cumtrapz(time_vec, vel_profile); 
end