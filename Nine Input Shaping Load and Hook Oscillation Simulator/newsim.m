%=========================================================================
% THESIS PLOT: ANGLE & VELOCITY ONLY (PNG EXPORT, NO POSITION)
%=========================================================================
clc; clearvars -except data01 data02 data03 data04 data05 data06 data07 data08 data09 data10 data11 data12 data13 data14 data15 data16; 
close all;

fprintf('========================================================\n');
fprintf('STARTING SIMULATION & EXPERIMENTAL PLOTTING (PNG FORMAT)\n');
fprintf('========================================================\n');

global m1 m2 l1 l2 g 
m1 = 0.15;  m2 = 0.36;  l1 = 0.3;   l2 = 0.2;   g  = 9.81;
dist_target = 0.8; v_target = 0.2; Ts = 0.01;

% Format settings
font_name = 'TH Saraban New';
font_size = 16;
text_color = [0 0 0]; 

% Colors & LineWidth
col_vel = 'b';      % ความเร็ว = น้ำเงิน
col_ang = 'g';      % มุมแกว่ง = เขียว
lw = 1;             % ความหนาเส้น = 1

% Folders Setup
desktopPath = fullfile(char(java.lang.System.getProperty('user.home')), 'Desktop', 'New Sim PNG');
f_sim_hook = fullfile(desktopPath, '1_Sim_Hook');
f_sim_load = fullfile(desktopPath, '2_Sim_Load');
f_exp_hook = fullfile(desktopPath, '3_Exp_Hook');
f_exp_load = fullfile(desktopPath, '4_Exp_Load');

folders = {desktopPath, f_sim_hook, f_sim_load, f_exp_hook, f_exp_load};
for f = 1:length(folders)
    if ~exist(folders{f}, 'dir'), mkdir(folders{f}); end
end

% Labels
file_names = {'IS01','IS02','IS03','IS04','IS05','IS06','IS07','IS08','IS09', ...
             'AC05','AC06','AC07','AC08','AC09','AC10'};

N = 15;
stop_time_exp = zeros(1, N);
stop_time_sim = zeros(1, N);
rms_exp_hook = zeros(1, N); rms_exp_load = zeros(1, N);
rms_sim_hook = zeros(1, N); rms_sim_load = zeros(1, N);
AllData = struct();

% Fixed Axis Limits (จากภาพ AC06.pdf)
xlim_val = [0, 20];
ylim_vel = [-0.3, 0.3];   % แกนซ้าย ความเร็ว
ylim_ang = [-5, 5];   % แกนขวา มุมแกว่ง

%% =========================================================================
%% STEP 1: DATA PROCESSING
%% =========================================================================
for k = 1:N
    fprintf('Processing [%02d/%d]: %s... ', k, N, file_names{k});
    
    %% --- Experimental Data ---
    data = eval(sprintf('data%02d', k));
    dt_exp = mean(diff(data.t));
    v_raw = [0; diff(data.x) ./ dt_exp];
    
    % Shift Time
    idx_start = find(abs(v_raw) > 0.005 & abs(data.x) > 0.001, 1, 'first');
    if isempty(idx_start), idx_start = 1; end
    t_exp_shifted = data.t(idx_start:end) - data.t(idx_start);
    v_raw = v_raw(idx_start:end);
    
    % Exp Stop Time
    [~, pk_idx] = max(v_raw);
    st_offset = find(abs(v_raw(pk_idx:end)) < 0.001, 1, 'first');
    if isempty(st_offset), idx_stop = pk_idx; else, idx_stop = pk_idx + st_offset - 1; end
    stop_time_exp(k) = t_exp_shifted(idx_stop);
    
    % Exp Angle Calculation
    raw_th1 = atan2d((data.y2(idx_start:end)), (data.x2(idx_start:end) - data.x(idx_start:end)));
    raw_th2 = atan2d((data.y1(idx_start:end)), (data.x1(idx_start:end) - data.x(idx_start:end)));
    th1_exp = -(raw_th1 - mean(raw_th1(idx_stop:end)));
    th2_exp = -(raw_th2 - mean(raw_th2(idx_stop:end)));
    
    rms_exp_hook(k) = sqrt(mean(th1_exp(idx_stop:end).^2));
    rms_exp_load(k) = sqrt(mean(th2_exp(idx_stop:end).^2));
    
    v_exp_sm10 = smoothdata(v_raw, 'gaussian', 15);
    th1_exp_smooth = smoothdata(th1_exp, 'gaussian', 15);
    th2_exp_smooth = smoothdata(th2_exp, 'gaussian', 15);
    
    %% --- Simulation Data ---
    if k <= 9
        [t_vec, pos_prof, vel_prof] = generate_profile_IS(dist_target, v_target, k, l1, l2, m1, m2, Ts);
    else
        st_coeff = (k - 5) / 10; 
        [t_vec, pos_prof, vel_prof] = generate_profile_ACCRT(dist_target, v_target, st_coeff, l1, l2, m1, m2, Ts);
    end
    
    accel_prof = gradient(vel_prof) ./ Ts;
    ode_func = @(t, Y) crane_dynamics(t, Y, t_vec, accel_prof, m1, m2, l1, l2, g);
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
    [t_ode, Y] = ode45(ode_func, [0 25], [0;0;0;0], options);
    
    t_ode = t_ode(:).'; Y = Y.';
    th1_sim = rad2deg(Y(1,:)); th2_sim = rad2deg(Y(2,:));
    v_sim = interp1(t_vec, vel_prof, t_ode, 'linear', 0);
    
    % Sim Stop Time
    [~, pk_sim] = max(v_sim);
    st_off_sim = find(abs(v_sim(pk_sim:end)) < 0.001, 1, 'first');
    if isempty(st_off_sim), idx_stop_sim = pk_sim; else, idx_stop_sim = pk_sim + st_off_sim - 1; end
    stop_time_sim(k) = t_ode(idx_stop_sim);
    
    rms_sim_hook(k) = sqrt(mean(th1_sim(idx_stop_sim:end).^2));
    rms_sim_load(k) = sqrt(mean(th2_sim(idx_stop_sim:end).^2));
    
    %% --- Store Data ---
    AllData(k).t_exp = t_exp_shifted; AllData(k).v_exp = v_exp_sm10;
    AllData(k).th1_exp = th1_exp_smooth;AllData(k).th2_exp = th2_exp_smooth;
    AllData(k).t_sim = t_ode; AllData(k).v_sim = v_sim;
    AllData(k).th1_sim = th1_sim; AllData(k).th2_sim = th2_sim;
    
    fprintf('Done.\n');
end

%% =========================================================================
%% STEP 2: PLOT GENERATION
%% =========================================================================

% --- 1. SET 1: Sim Hook ---
for k = 1:N
    fig = figure('Name', ['Sim_Hook_' file_names{k}], 'Color', 'w', 'Position', [100 100 800 450], 'Visible', 'off');
    yyaxis left; hold on;
    plot(AllData(k).t_sim, AllData(k).v_sim, '-', 'Color', col_vel, 'LineWidth', lw);
    ylim(ylim_vel); ax = gca; ax.YColor = text_color; 
    
    yyaxis right; hold on;
    plot(AllData(k).t_sim, AllData(k).th1_sim, '-', 'Color', col_ang, 'LineWidth', lw);
    ylim(ylim_ang); ax.YColor = text_color;
    
    format_base_ax(ax, xlim_val, font_name, font_size, text_color);
    set_top_ylabels(ax, 'Velocity (m/s)', 'Angle (deg)', text_color, font_name, font_size);
    exportgraphics(fig, fullfile(f_sim_hook, [file_names{k} '.png']), 'Resolution', 300);
    close(fig);
end

% --- 2. SET 2: Sim Load ---
for k = 1:N
    fig = figure('Name', ['Sim_Load_' file_names{k}], 'Color', 'w', 'Position', [150 150 800 450], 'Visible', 'off');
    yyaxis left; hold on;
    plot(AllData(k).t_sim, AllData(k).v_sim, '-', 'Color', col_vel, 'LineWidth', lw);
    ylim(ylim_vel); ax = gca; ax.YColor = text_color; 
    
    yyaxis right; hold on;
    plot(AllData(k).t_sim, AllData(k).th2_sim, '-', 'Color', col_ang, 'LineWidth', lw);
    ylim(ylim_ang); ax.YColor = text_color;
    
    format_base_ax(ax, xlim_val, font_name, font_size, text_color);
    set_top_ylabels(ax, 'Velocity (m/s)', 'Angle (deg)', text_color, font_name, font_size);
    exportgraphics(fig, fullfile(f_sim_load, [file_names{k} '.png']), 'Resolution', 300);
    close(fig);
end

% --- 3. SET 3: Exp Hook ---
for k = 1:N
    fig = figure('Name', ['Exp_Hook_' file_names{k}], 'Color', 'w', 'Position', [200 200 800 450], 'Visible', 'off');
    yyaxis left; hold on;
    plot(AllData(k).t_exp, AllData(k).v_exp, '-', 'Color', col_vel, 'LineWidth', lw);
    ylim(ylim_vel); ax = gca; ax.YColor = text_color; 
    
    yyaxis right; hold on;
    plot(AllData(k).t_exp, AllData(k).th1_exp, '-', 'Color', col_ang, 'LineWidth', lw);
    ylim(ylim_ang); ax.YColor = text_color;
    
    format_base_ax(ax, xlim_val, font_name, font_size, text_color);
    set_top_ylabels(ax, 'Velocity (m/s)', 'Angle (deg)', text_color, font_name, font_size);
    exportgraphics(fig, fullfile(f_exp_hook, [file_names{k} '.png']), 'Resolution', 300);
    close(fig);
end

% --- 4. SET 4: Exp Load ---
for k = 1:N
    fig = figure('Name', ['Exp_Load_' file_names{k}], 'Color', 'w', 'Position', [250 250 800 450], 'Visible', 'off');
    yyaxis left; hold on;
    plot(AllData(k).t_exp, AllData(k).v_exp, '-', 'Color', col_vel, 'LineWidth', lw);
    ylim(ylim_vel); ax = gca; ax.YColor = text_color; 
    
    yyaxis right; hold on;
    plot(AllData(k).t_exp, AllData(k).th2_exp, '-', 'Color', col_ang, 'LineWidth', lw);
    ylim(ylim_ang); ax.YColor = text_color;
    
    format_base_ax(ax, xlim_val, font_name, font_size, text_color);
    set_top_ylabels(ax, 'Velocity (m/s)', 'Angle (deg)', text_color, font_name, font_size);
    exportgraphics(fig, fullfile(f_exp_load, [file_names{k} '.png']), 'Resolution', 300);
    close(fig);
end

%% =========================================================================
%% EXPORT EXCEL DATA
%% =========================================================================
summary_table = table(file_names', ...
    stop_time_sim', stop_time_exp', ...
    rms_sim_hook', rms_exp_hook', ...
    rms_sim_load', rms_exp_load', ...
    'VariableNames', {'Method_Name', ...
    'Sim_StopTime_sec', 'Exp_StopTime_sec', ...
    'Sim_RMS_Hook_deg', 'Exp_RMS_Hook_deg', ...
    'Sim_RMS_Load_deg', 'Exp_RMS_Load_deg'});

excel_filename = fullfile(desktopPath, 'Results_Summary.xlsx');
writetable(summary_table, excel_filename);

fprintf('========================================================\n');
fprintf('SUCCESS! 60 PNG Plots & Excel exported to:\n%s\n', desktopPath);
fprintf('========================================================\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================
function align_zeros_yyaxis(ax)
    yyaxis(ax,'left'); limL = ylim(ax);
    yyaxis(ax,'right'); limR = ylim(ax);
    if limL(1) >= 0 || limL(2) <= 0 || limR(1) >= 0 || limR(2) <= 0
        return; 
    end
    fracL = limL(2) / (limL(2) - limL(1));
    fracR = limR(2) / (limR(2) - limR(1));
    if fracL > fracR
        limR(2) = (fracL / (1 - fracL)) * (-limR(1));
    else
        limL(2) = (fracR / (1 - fracR)) * (-limL(1));
    end
    yyaxis(ax,'left'); ylim(ax, limL);
    yyaxis(ax,'right'); ylim(ax, limR);
end

function format_base_ax(ax, xlim_val, f_name, f_size, t_color)
    % 1. จัดการแกน X (สั่งให้เดินทีละ 2 วินาที)
    if ~strcmp(xlim_val, 'none')
        xlim(ax, xlim_val); 
        xticks(ax, xlim_val(1):2:xlim_val(2)); % <--- เพิ่มการตั้งค่า Grid แกน X
        xlabel(ax, 'Time (sec)', 'Color', t_color); 
    end
    
    ax.Color = 'w'; 
    ax.XColor = t_color; 
    
    % 2. บังคับสีแกน Y และจัดการสเกลแกนขวา (สั่งให้เดินทีละ 1 องศา)
    if length(ax.YAxis) == 2
        % สลับมาจัดการแกนขวา (Angle)
        yyaxis(ax, 'right');
        yl_r = ylim(ax);
        yticks(ax, floor(yl_r(1)):1:ceil(yl_r(2))); % <--- เพิ่มการตั้งค่า Grid แกน Y ขวา
        
        ax.YAxis(1).Color = t_color;
        ax.YAxis(2).Color = t_color;
    else
        ax.YColor = t_color;
    end
    
    % 3. ตั้งค่าฟอนต์และความโปร่งของกราฟ
    ax.FontName = f_name; 
    ax.FontSize = f_size;
    ax.LineWidth = 0.5; 
    
    % 4. เปิดเส้น Grid (จะอิงตามสเกล 2 วินาที และ 1 องศาที่เราตั้งไว้)
    grid(ax, 'on'); 
    ax.GridColor = [0 0 0];
    ax.GridAlpha = 0.1; 
    
    % 5. วาดกรอบกราฟ
    box(ax, 'on');
    set(ax, 'Box', 'on');
end

function set_top_ylabels(ax, left_txt, right_txt, t_color, f_name, f_size)
    yyaxis(ax, 'left'); ylabel(ax, '');
    % นำความหนา (Bold) ออก ขยับขึ้นนิดหน่อยเป็น 1.03 ให้มีพื้นที่หายใจ
    text(ax, 0, 1.03, left_txt, 'Units', 'normalized', 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', 'FontName', f_name, 'FontSize', f_size, 'Color', t_color);
    
    yyaxis(ax, 'right'); ylabel(ax, '');
    text(ax, 1, 1.03, right_txt, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom', 'FontName', f_name, 'FontSize', f_size, 'Color', t_color);
end

function dYdt = crane_dynamics(t, Y, time_vec, accel_profile, m1, m2, l1, l2, g)
    a_trolley = interp1(time_vec, accel_profile, t, 'linear', 0);
    M = [(m1+m2)*l1^2, m2*l1*l2; m2*l1*l2, m2*l2^2];
    K = [(m1+m2)*g*l1, 0; 0, m2*g*l2];
    B = [-(m1+m2)*l1; -m2*l2];
    theta = [Y(1); Y(2)]; dtheta = [Y(3); Y(4)];
    ddtheta = M \ (B * a_trolley - K * theta);
    dYdt = [dtheta; ddtheta];
end