function shaper_vec = shaper_07_zvdd_zv(T1, T2, Ts)
    % Mode 1: ZVDD
    t1 = [0, 0.5*T1, T1, 1.5*T1];
    a1 = [0.125, 0.375, 0.375, 0.125];
    
    % Mode 2: ZV
    t2 = [0, 0.5*T2];
    a2 = [0.5, 0.5];
    
    v1 = times2vec(t1, a1, Ts);
    v2 = times2vec(t2, a2, Ts);
    shaper_vec = conv(v1, v2);
end

function vec = times2vec(times, amps, Ts)
    idx = round(times/Ts) + 1;
    vec = zeros(1, max(idx));
    vec(idx) = amps;
end