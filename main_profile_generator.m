clear; clc; close all;

%% --- SECTION 1: CONFIGURATION ---
% เลือก Shaper (1-9)
selected_file_number = 9; 

% Parameters
lh = 0.2; ll = 0.3; mh = 0.01; ml = 1.5;
Ts = 0.01;      % Sample Time

v_target = 0.2;  
dist_target = 0.8;  
start_time = 0.0;

% Motor Config
radius_cm = 4.81;          
dist_per_rev = pi * (radius_cm / 100); 
steps_per_rev = 1600;

%% --- SECTION 2: IDENTIFICATION ---
g = 9.81;
a = mh * lh * ll;
b = (mh * g * ll) + (ml * g * ll) + (mh * g * lh) + (ml * g * lh);
c = (mh * g * g) + (ml * g * g);
discriminant = sqrt(b * b - 4 * a * c);
wn1 = sqrt((b - discriminant) / (2 * a));
wn2 = sqrt((b + discriminant) / (2 * a));
T1 = 2 * pi / wn1;
T2 = 2 * pi / wn2;

fprintf('Generating Profile #0%d\n', selected_file_number);
fprintf('Target: %.4f m @ %.2f m/s\n', dist_target, v_target);

%% --- SECTION 3: GENERATE FIR FILTER ---
N = 9;
n_filter = 0:N-1;
sigma = (N-1)/6;
g_gauss = exp(-0.5*((n_filter-(N-1)/2)/sigma).^2);
h = g_gauss / sum(g_gauss);

%% --- SECTION 4: CALL SHAPER FUNCTION ---
switch selected_file_number
    case 1, shaper_vec = shaper_01_zv_zv(T1, T2, Ts);
    case 2, shaper_vec = shaper_02_zv_zvd(T1, T2, Ts);
    case 3, shaper_vec = shaper_03_zv_zvdd(T1, T2, Ts);
    case 4, shaper_vec = shaper_04_zvd_zv(T1, T2, Ts);
    case 5, shaper_vec = shaper_05_zvd_zvd(T1, T2, Ts);
    case 6, shaper_vec = shaper_06_zvd_zvdd(T1, T2, Ts);
    case 7, shaper_vec = shaper_07_zvdd_zv(T1, T2, Ts);
    case 8, shaper_vec = shaper_08_zvdd_zvd(T1, T2, Ts);
    case 9, shaper_vec = shaper_09_zvdd_zvdd(T1, T2, Ts);
    otherwise, error('Invalid File Number!');
end

%% --- SECTION 5: AUTOMATIC PULSE GENERATOR ---
% 1. คำนวณเวลา Pulse จากระยะทาง
pulse_duration = dist_target / v_target;

% 2. คำนวณเวลาจบ Simulation
shaper_duration = length(shaper_vec) * Ts;
filter_duration = length(h) * Ts;
end_time = start_time + pulse_duration;
t_final = end_time + shaper_duration + filter_duration + 1.0; 

% 3. สร้าง Command
time_vec = 0:Ts:t_final;
v_command = zeros(size(time_vec));
v_command(time_vec >= start_time & time_vec < end_time) = v_target;

%% --- SECTION 6: CONVOLUTION & PROCESS ---
v_shaped = conv(v_command, shaper_vec);
v_shaped = v_shaped(1:length(time_vec)); 

v_filtered = conv(v_shaped, h);
v_filtered = v_filtered(1:length(time_vec)); 

pos_profile = cumtrapz(time_vec, v_filtered); 
freq_profile = (v_filtered / dist_per_rev) * steps_per_rev;
freq_profile_int = round(freq_profile);

%% --- SECTION 7: PLOT (SEPARATE FIGURES) ---

% 7.1 Velocity (เทียบ Command vs Shaped)
figure('Name', '1. Velocity Profile', 'NumberTitle', 'off');
plot(time_vec, v_command, 'k--', 'LineWidth', 1); hold on;
plot(time_vec, v_filtered, 'b', 'LineWidth', 2);
title('Velocity Profile (Command vs Shaped)');
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('Base Pulse', 'Shaped Profile'); grid on;

% 7.2 Displacement (ระยะทาง)
figure('Name', '2. Position Profile', 'NumberTitle', 'off');
plot(time_vec, pos_profile, 'g', 'LineWidth', 2);
yline(dist_target, 'r--', 'Target Distance');
title(sprintf('Position Profile (Target: %.3fm, Actual: %.3fm)', dist_target, pos_profile(end)));
xlabel('Time (s)'); ylabel('Position (m)');
grid on;

% 7.3 Frequency (สำหรับ Motor)
figure('Name', '3. Frequency Command', 'NumberTitle', 'off');
plot(time_vec, freq_profile_int, 'r', 'LineWidth', 1.5);
title('Frequency Command for Stepper Motor');
xlabel('Time (s)'); ylabel('Frequency (Hz)');
max_freq = max(freq_profile_int);
text(start_time, max_freq, sprintf(' Max Freq: %d Hz', max_freq), 'VerticalAlignment', 'bottom');
grid on;

%% --- SECTION 8: EXPORT ---
filename = sprintf('motor_profile_0%d.txt', selected_file_number);
fileID = fopen(filename, 'w');
fprintf(fileID, '// Profile #%d: Dist=%.3fm, V=%.2fm/s\n', selected_file_number, dist_target, v_target);
fprintf(fileID, 'const int speed_profile[] PROGMEM = {\n  ');
for i = 1:length(freq_profile_int)
    fprintf(fileID, '%d', freq_profile_int(i));
    if i < length(freq_profile_int), fprintf(fileID, ', '); end
    if mod(i, 20) == 0, fprintf(fileID, '\n  '); end
end
fprintf(fileID, '\n};\n');
fprintf(fileID, 'const long profile_length = %d;\n', length(freq_profile_int));
fclose(fileID);
fprintf('Exported to %s\n', filename);