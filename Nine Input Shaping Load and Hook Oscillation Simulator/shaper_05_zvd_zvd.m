function shaper_vec = shaper_05_zvd_zvd(T1, T2, Ts)
    % Mode 1: ZVD
    t1 = [0, 0.5*T1, T1];
    a1 = [0.25, 0.5, 0.25];
    
    % Mode 2: ZVD
    t2 = [0, 0.5*T2, T2];
    a2 = [0.25, 0.5, 0.25];
    
    v1 = times2vec(t1, a1, Ts);
    v2 = times2vec(t2, a2, Ts);
    shaper_vec = conv(v1, v2);
end

function vec = times2vec(times, amps, Ts)
    idx = round(times/Ts) + 1;
    vec = zeros(1, max(idx));
    vec(idx) = amps;
end