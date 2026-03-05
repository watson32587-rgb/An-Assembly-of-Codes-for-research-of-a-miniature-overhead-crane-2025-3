function [time_vec, pos_profile, vel_profile] = generate_profile_ACCRT(dist_target, v_target, St, lh, ll, mh, ml, Ts)
    % 1. System Identification
    [wn1, wn2, b1, b2] = ModePeriodize(lh, ll, mh, ml);
    
    % 2. Solver & Jacobian
    Tc = St * 2 * pi / wn1;
    Cq = AcJacobian(wn1, wn2, b1, b2, Tc);
    BC = AcBC(v_target);
    C = Cq \ BC;
    
    % 3. Calculate Timing
    ta = findTa(v_target, dist_target, C, b1, b2, wn1, wn2, Tc);
    if ta < 0
        warning('Distance too short for this velocity! ta is set to 0.');
        ta = 0; 
    end
    
    % 4. Generate Trajectory
    d1_idx = round(Tc / Ts);
    d2_idx = round(ta / Ts);
    total_samples = 2*d1_idx + d2_idx; 
    
    time_vec = (0:total_samples-1) * Ts;
    accel_profile = zeros(size(time_vec));
    vel_profile = zeros(size(time_vec));
    
    for k = 1:total_samples
        t_sim = (k-1) * Ts;
        if k <= d1_idx
            accel_profile(k) = ACshaped_Eq(t_sim, C, b1, b2, wn1, wn2);
            vel_profile(k) = ACshaped_V(t_sim, C, b1, b2, wn1, wn2); 
        elseif k > d1_idx && k <= (d1_idx + d2_idx)
            accel_profile(k) = 0;
            vel_profile(k) = v_target;
        elseif k > (d1_idx + d2_idx)
            k_ref = (2*d1_idx + d2_idx) - k + 1;
            if k_ref >= 1 && k_ref <= d1_idx
                accel_profile(k) = -accel_profile(k_ref);
                vel_profile(k) = vel_profile(k_ref); 
            end
        end
    end
    % Integrate Velocity -> Position
    pos_profile = cumtrapz(time_vec, vel_profile);
end

%% ================= LOCAL FUNCTIONS =================
function [wn1, wn2, b1, b2] = ModePeriodize(lh, ll, mh, ml)
    g = 9.81;
    M_s = [(mh+ml)*lh^2, ml*lh*ll; ml*lh*ll, ml*ll^2];
    K_s = [(mh+ml)*g*lh, 0; 0, ml*g*ll];
    B_s = [-(mh+ml)*lh; -ml*ll];
    [mds_s, lambda_s] = eig(K_s, M_s);
    omg_s = sqrt(diag(lambda_s));
    b_nu_s = mds_s.' * B_s; 
    wn1 = real(omg_s(1)); wn2 = real(omg_s(2));
    b1 = real(b_nu_s(1)); b2 = real(b_nu_s(2));
end

function BC = AcBC(u)
    BC = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; u];
end

function yy = ACshaped_Eq(t, C, b1, b2, wn1, wn2)
    a1=C(1); a2=C(2); a3=C(3); a4=C(4); a5=C(5); a6=C(6);
    yy = a6-t.*(-a5+a2.*b1.*wn1+a4.*b2.*wn2)-a1.*b1-a3.*b2+1.0./wn1.^2.* ...
         (a1.*b1.*wn1.^2.*cos(wn1.*t)+a2.*b1.*wn1.^2.*sin(wn1.*t))+1.0./wn2.^2.* ...
         (a3.*b2.*wn2.^2.*cos(wn2.*t)+a4.*b2.*wn2.^2.*sin(wn2.*t)) ; 
end

function y = ACshaped_V(t, C, b1, b2, w1, w2)
    a1=C(1); a2=C(2); a3=C(3); a4=C(4); a5=C(5); a6=C(6); a7=C(7);
    y = -1.0./w1.^2.*(a2.*b1.*w1.*cos(w1.*t) - a1.*b1.*w1.*sin(w1.*t)) ...
      -1.0./w2.^2.*(a4.*b2.*w2.*cos(w2.*t) - a3.*b2.*w2.*sin(w2.*t)) ...
      -t.*(-a6 + a1.*b1 + a3.*b2) ...
      -(t.^2.*(-a5 + a2.*b1.*w1 + a4.*b2.*w2))./2.0 ...
      +(a2.*b1.*w2 + a4.*b2.*w1 + a7.*w1.*w2)./(w1.*w2);
end

function ta = findTa(u, x, C, b1, b2, wn1, wn2, Tc)
    a1=C(1); a2=C(2); a3=C(3); a4=C(4); a5=C(5); a6=C(6); a7=C(7); a8=C(8); d1 = Tc;
    y = -1.0/wn1^2 * (a1*b1*cos(wn1*d1) + a2*b1*sin(wn1*d1)) ...
        -1.0/wn2^2 * (a3*b2*cos(wn2*d1) + a4*b2*sin(wn2*d1)) ...
        -(d1^3 * (-a5 + a2*b1*wn1 + a4*b2*wn2))/6.0 ...
        -(d1^2 * (-a6 + a1*b1 + a3*b2))/2.0 ...
        +1.0/wn1^2 * 1.0/wn2^2 * (a8*wn1^2*wn2^2 + a1*b1*wn2^2 + a3*b2*wn1^2) ...
        +(d1 * (a2*b1*wn2 + a4*b2*wn1 + a7*wn1*wn2))/(wn1*wn2);
    ta = (x - 2*y) / u;
end

function Cq = AcJacobian(wn1,wn2,b1,b2,Tc)
    mt1 = [-b1.^2.*1.0./wn1.^2+(1.0./wn1.^3.*(b1.^2.*wn1.^3-b1.^2.*wn1.*wn2.^2))./((wn1+wn2).*(wn1-wn2)),-b1.*b2.*1.0./wn2.^2-(b1.*b2)./((wn1+wn2).*(wn1-wn2))+(b1.*b2.*wn1.^2.*1.0./wn2.^2)./((wn1+wn2).*(wn1-wn2)),0.0,0.0,0.0,0.0,-b1.^2.*1.0./wn1.^2+(1.0./wn1.^3.*(b1.^2.*wn1.^3.*cos(Tc.*wn1)-b1.^2.*wn1.*wn2.^2.*cos(Tc.*wn1)))./((wn1+wn2).*(wn1-wn2))+(Tc.*b1.^2.*sin(Tc.*wn1))./(wn1.*2.0)];
    mt2 = [-b1.*b2.*1.0./wn2.^2-(b1.*b2.*cos(Tc.*wn1))./((wn1+wn2).*(wn1-wn2))+(b1.*b2.*wn1.^2.*1.0./wn2.^2.*cos(Tc.*wn2))./((wn1+wn2).*(wn1-wn2)),0.0,(b1.^2.*1.0./wn1.^2.*(wn1.*sin(Tc.*wn1)+Tc.*wn1.^2.*cos(Tc.*wn1)))./2.0-(1.0./wn1.^3.*(b1.^2.*wn1.^4.*sin(Tc.*wn1)-b1.^2.*wn1.^2.*wn2.^2.*sin(Tc.*wn1)))./((wn1+wn2).*(wn1-wn2)),(b1.*b2.*wn1.*sin(Tc.*wn1))./((wn1+wn2).*(wn1-wn2))-(b1.*b2.*wn1.^2.*sin(Tc.*wn2))./(wn2.*(wn1+wn2).*(wn1-wn2))];
    mt3 = [-Tc.*b1+(b1.*sin(Tc.*wn1))./wn1,0.0,0.0,0.0,-b1.^2./wn1+(1.0./wn1.^3.*(b1.^2.*wn1.^4-b1.^2.*wn1.^2.*wn2.^2))./((wn1+wn2).*(wn1-wn2)),-b1.*b2.*wn1.*1.0./wn2.^2-(b1.*b2.*wn1)./((wn1+wn2).*(wn1-wn2))+(b1.*b2.*wn1.^3.*1.0./wn2.^2)./((wn1+wn2).*(wn1-wn2)),0.0];
    mt4 = [-(Tc.*b1.^2)./wn1+(b1.^2.*1.0./wn1.^2.*(sin(Tc.*wn1)-Tc.*wn1.*cos(Tc.*wn1)))./2.0+(1.0./wn1.^3.*(b1.^2.*wn1.^3.*sin(Tc.*wn1)-b1.^2.*wn1.*wn2.^2.*sin(Tc.*wn1)))./((wn1+wn2).*(wn1-wn2)),-(b1.*b2.*sin(Tc.*wn1))./((wn1+wn2).*(wn1-wn2))-Tc.*b1.*b2.*wn1.*1.0./wn2.^2+(b1.*b2.*wn1.^3.*1.0./wn2.^3.*sin(Tc.*wn2))./((wn1+wn2).*(wn1-wn2)),b1.*wn1];
    mt5 = [-b1.^2./wn1+(Tc.*b1.^2.*sin(Tc.*wn1))./2.0+(1.0./wn1.^3.*(b1.^2.*wn1.^4.*cos(Tc.*wn1)-b1.^2.*wn1.^2.*wn2.^2.*cos(Tc.*wn1)))./((wn1+wn2).*(wn1-wn2)),-b1.*b2.*wn1.*1.0./wn2.^2-(b1.*b2.*wn1.*cos(Tc.*wn1))./((wn1+wn2).*(wn1-wn2))+(b1.*b2.*wn1.^3.*1.0./wn2.^2.*cos(Tc.*wn2))./((wn1+wn2).*(wn1-wn2)),b1./wn1-(b1.*cos(Tc.*wn1))./wn1-(Tc.^2.*b1.*wn1)./2.0];
    mt6 = [-b1.*b2.*1.0./wn1.^2+(b1.*b2)./((wn1+wn2).*(wn1-wn2))-(b1.*b2.*1.0./wn1.^2.*wn2.^2)./((wn1+wn2).*(wn1-wn2)),-b2.^2.*1.0./wn2.^2-(1.0./wn2.^3.*(b2.^2.*wn2.^3-b2.^2.*wn1.^2.*wn2))./((wn1+wn2).*(wn1-wn2)),0.0,0.0,0.0,0.0,-b1.*b2.*1.0./wn1.^2+(b1.*b2.*cos(Tc.*wn2))./((wn1+wn2).*(wn1-wn2))-(b1.*b2.*1.0./wn1.^2.*wn2.^2.*cos(Tc.*wn1))./((wn1+wn2).*(wn1-wn2))];
    mt7 = [-b2.^2.*1.0./wn2.^2-(1.0./wn2.^3.*(b2.^2.*wn2.^3.*cos(Tc.*wn2)-b2.^2.*wn1.^2.*wn2.*cos(Tc.*wn2)))./((wn1+wn2).*(wn1-wn2))+(Tc.*b2.^2.*sin(Tc.*wn2))./(wn2.*2.0),0.0,-(b1.*b2.*wn2.*sin(Tc.*wn2))./((wn1+wn2).*(wn1-wn2))+(b1.*b2.*wn2.^2.*sin(Tc.*wn1))./(wn1.*(wn1+wn2).*(wn1-wn2))];
    mt8 = [(b2.^2.*1.0./wn2.^2.*(wn2.*sin(Tc.*wn2)+Tc.*wn2.^2.*cos(Tc.*wn2)))./2.0+(1.0./wn2.^3.*(b2.^2.*wn2.^4.*sin(Tc.*wn2)-b2.^2.*wn1.^2.*wn2.^2.*sin(Tc.*wn2)))./((wn1+wn2).*(wn1-wn2)),-Tc.*b2+(b2.*sin(Tc.*wn2))./wn2,0.0,0.0,0.0,-b1.*b2.*1.0./wn1.^2.*wn2+(b1.*b2.*wn2)./((wn1+wn2).*(wn1-wn2))-(b1.*b2.*1.0./wn1.^2.*wn2.^3)./((wn1+wn2).*(wn1-wn2))];
    mt9 = [-b2.^2./wn2-(1.0./wn2.^3.*(b2.^2.*wn2.^4-b2.^2.*wn1.^2.*wn2.^2))./((wn1+wn2).*(wn1-wn2)),0.0,(b1.*b2.*sin(Tc.*wn2))./((wn1+wn2).*(wn1-wn2))-Tc.*b1.*b2.*1.0./wn1.^2.*wn2-(b1.*b2.*1.0./wn1.^3.*wn2.^3.*sin(Tc.*wn1))./((wn1+wn2).*(wn1-wn2))];
    mt10 = [-(Tc.*b2.^2)./wn2+(b2.^2.*1.0./wn2.^2.*(sin(Tc.*wn2)-Tc.*wn2.*cos(Tc.*wn2)))./2.0-(1.0./wn2.^3.*(b2.^2.*wn2.^3.*sin(Tc.*wn2)-b2.^2.*wn1.^2.*wn2.*sin(Tc.*wn2)))./((wn1+wn2).*(wn1-wn2)),b2.*wn2,-b1.*b2.*1.0./wn1.^2.*wn2+(b1.*b2.*wn2.*cos(Tc.*wn2))./((wn1+wn2).*(wn1-wn2))-(b1.*b2.*1.0./wn1.^2.*wn2.^3.*cos(Tc.*wn1))./((wn1+wn2).*(wn1-wn2))];
    mt11 = [-b2.^2./wn2+(Tc.*b2.^2.*sin(Tc.*wn2))./2.0-(1.0./wn2.^3.*(b2.^2.*wn2.^4.*cos(Tc.*wn2)-b2.^2.*wn1.^2.*wn2.^2.*cos(Tc.*wn2)))./((wn1+wn2).*(wn1-wn2)),b2./wn2-(b2.*cos(Tc.*wn2))./wn2-(Tc.^2.*b2.*wn2)./2.0,0.0,0.0,0.0,b1.*1.0./wn1.^2-(1.0./wn1.^3.*(b1.*wn1.^3-b1.*wn1.*wn2.^2))./((wn1+wn2).*(wn1-wn2)),b2.*1.0./wn2.^2+(1.0./wn2.^3.*(b2.*wn2.^3-b2.*wn1.^2.*wn2))./((wn1+wn2).*(wn1-wn2)),0.0];
    mt12 = [Tc.*b1.*1.0./wn1.^2-(1.0./wn1.^3.*(b1.*wn1.^2.*sin(Tc.*wn1)-b1.*wn2.^2.*sin(Tc.*wn1)))./((wn1+wn2).*(wn1-wn2)),Tc.*b2.*1.0./wn2.^2-(1.0./wn2.^3.*(b2.*wn1.^2.*sin(Tc.*wn2)-b2.*wn2.^2.*sin(Tc.*wn2)))./((wn1+wn2).*(wn1-wn2)),-1.0,b1.*1.0./wn1.^2-(1.0./wn1.^3.*(b1.*wn1.^3.*cos(Tc.*wn1)-b1.*wn1.*wn2.^2.*cos(Tc.*wn1)))./((wn1+wn2).*(wn1-wn2))];
    mt13 = [b2.*1.0./wn2.^2+(1.0./wn2.^3.*(b2.*wn2.^3.*cos(Tc.*wn2)-b2.*wn1.^2.*wn2.*cos(Tc.*wn2)))./((wn1+wn2).*(wn1-wn2)),Tc.^2./2.0,b1.*1.0./wn1.^2-(1.0./wn1.^3.*(b1.*wn1.^3-b1.*wn1.*wn2.^2))./((wn1+wn2).*(wn1-wn2)),b2.*1.0./wn2.^2+(1.0./wn2.^3.*(b2.*wn2.^3-b2.*wn1.^2.*wn2))./((wn1+wn2).*(wn1-wn2)),0.0,0.0,0.0,0.0,b1.*1.0./wn1.^2-(1.0./wn1.^3.*(b1.*wn1.^3.*cos(Tc.*wn1)-b1.*wn1.*wn2.^2.*cos(Tc.*wn1)))./((wn1+wn2).*(wn1-wn2))];
    mt14 = [b2.*1.0./wn2.^2+(1.0./wn2.^3.*(b2.*wn2.^3.*cos(Tc.*wn2)-b2.*wn1.^2.*wn2.*cos(Tc.*wn2)))./((wn1+wn2).*(wn1-wn2)),0.0,(1.0./wn1.^3.*(b1.*wn1.^4.*sin(Tc.*wn1)-b1.*wn1.^2.*wn2.^2.*sin(Tc.*wn1)))./((wn1+wn2).*(wn1-wn2)),-(1.0./wn2.^3.*(b2.*wn2.^4.*sin(Tc.*wn2)-b2.*wn1.^2.*wn2.^2.*sin(Tc.*wn2)))./((wn1+wn2).*(wn1-wn2)),Tc,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
    mt15 = [(1.0./wn1.^3.*(wn1.^5-wn1.^3.*wn2.^2))./((wn1+wn2).*(wn1-wn2)),0.0,0.0,(1.0./wn1.^3.*(wn1.^4.*sin(Tc.*wn1)-wn1.^2.*wn2.^2.*sin(Tc.*wn1)))./((wn1+wn2).*(wn1-wn2)),0.0,0.0,(1.0./wn1.^3.*(wn1.^5.*cos(Tc.*wn1)-wn1.^3.*wn2.^2.*cos(Tc.*wn1)))./((wn1+wn2).*(wn1-wn2)),0.0,0.0,(1.0./wn1.^3.*(wn1.^5-wn1.^3.*wn2.^2))./((wn1+wn2).*(wn1-wn2)),0.0,0.0,0.0,0.0,0.0];
    mt16 = [(1.0./wn1.^3.*(wn1.^5.*cos(Tc.*wn1)-wn1.^3.*wn2.^2.*cos(Tc.*wn1)))./((wn1+wn2).*(wn1-wn2)),0.0,0.0,-(1.0./wn1.^3.*(wn1.^6.*sin(Tc.*wn1)-wn1.^4.*wn2.^2.*sin(Tc.*wn1)))./((wn1+wn2).*(wn1-wn2)),0.0,0.0,0.0,0.0,0.0,0.0,-(1.0./wn2.^3.*(wn2.^5-wn1.^2.*wn2.^3))./((wn1+wn2).*(wn1-wn2)),0.0,0.0,-(1.0./wn2.^3.*(wn2.^4.*sin(Tc.*wn2)-wn1.^2.*wn2.^2.*sin(Tc.*wn2)))./((wn1+wn2).*(wn1-wn2)),0.0,0.0];
    mt17 = [-(1.0./wn2.^3.*(wn2.^5.*cos(Tc.*wn2)-wn1.^2.*wn2.^3.*cos(Tc.*wn2)))./((wn1+wn2).*(wn1-wn2)),0.0,0.0,-(1.0./wn2.^3.*(wn2.^5-wn1.^2.*wn2.^3))./((wn1+wn2).*(wn1-wn2)),0.0,0.0,0.0,0.0,0.0,-(1.0./wn2.^3.*(wn2.^5.*cos(Tc.*wn2)-wn1.^2.*wn2.^3.*cos(Tc.*wn2)))./((wn1+wn2).*(wn1-wn2)),0.0,0.0,(1.0./wn2.^3.*(wn2.^6.*sin(Tc.*wn2)-wn1.^2.*wn2.^4.*sin(Tc.*wn2)))./((wn1+wn2).*(wn1-wn2)),0.0];
    Cq = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8,mt9,mt10,mt11,mt12,mt13,mt14,mt15,mt16,mt17],12,12);
end