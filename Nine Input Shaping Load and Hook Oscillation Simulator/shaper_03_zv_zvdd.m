function shaper_vec = shaper_03_zv_zvdd(T1, T2, Ts)
    % Mode 1: ZV
    t1 = [0, 0.5*T1];
    a1 = [0.5, 0.5];
    
    % Mode 2: ZVDD (4 impulses)
    t2 = [0, 0.5*T2, T2, 1.5*T2];
    a2 = [0.125, 0.375, 0.375, 0.125];
    
    v1 = times2vec(t1, a1, Ts);
    v2 = times2vec(t2, a2, Ts);
    shaper_vec = conv(v1, v2);
end

function vec = times2vec(times, amps, Ts)
    idx = round(times/Ts) + 1;
    vec = zeros(1, max(idx));
    vec(idx) = amps;
end