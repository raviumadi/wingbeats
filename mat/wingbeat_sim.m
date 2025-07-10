fs = 192e3;
bandwidth = 65e3:-100:35e3;
c = 343;

params = struct();
params.fs = fs;
params.bandwidth = bandwidth;
params.kr = 5;
params.target_distance = 20;  % m
params.initial_velocity = 2;  % m/s
params.initial_call_duration = 0.005;  % s
params.f_wing = 5;  % Hz
params.theta = pi / 4;  % wingbeat excursion (adjust as needed)
params.theta_mode = 'dynamic'; % 'dynamic' or 'fixed'
params.max_wingbeat_freq = 3*params.f_wing;
params.min_theta = 0.3 * params.theta;
params.motile = 0;
params.makeAudio = 1;
params.wing_sync_mode = 'dynamic'; % 'dynamic' or 'fixed'

result = simulateEcholocationWings(params);

%%
% Assumes you already ran: result = simulateEcholocationWings(params);

phi_star = result.phi_star;
wingbeat_phase = result.wingbeat_phase;
synchrony_flag = result.synchrony_flag;
actual_synchrony_flag = result.actual_synchrony_flag;
call_rate = result.call_rate;
call_idx = 1:length(phi_star);
theta = params.theta;

% Align call_rate (one element shorter) for plotting
% call_rate = [call_rate(1), call_rate];
% call_rate(1) = [];

figure('Position', [100, 100, 1200, 500]);

yyaxis left
plot(call_idx, rad2deg(phi_star), 'r-', 'LineWidth', 2); hold on;
plot(call_idx, rad2deg(wingbeat_phase), 'k--', 'LineWidth', 2);
    plot(call_idx, rad2deg(wingbeat_phase), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
yline(rad2deg(theta), 'k:', 'Max Wing Excursion $\theta$', 'LineWidth', 2,'Interpreter', 'latex', 'FontSize', 16, 'Fontweight', 'bold');

% Show *all* predicted asynchronous calls (synchrony_flag==0) as gray dots
async_pred_idx = find(synchrony_flag == 0);
scatter(call_idx(async_pred_idx), rad2deg(phi_star(async_pred_idx)), ...
    20, [0.4 0.4 0.4], 'filled', 'DisplayName', 'Pred. Asynchronous');

% Highlight *actual* synchronous calls as green circles
actual_sync_idx = find(actual_synchrony_flag == 1);
scatter(call_idx(actual_sync_idx), rad2deg(phi_star(actual_sync_idx)), ...
    70, 'g', 'o', 'LineWidth', 2, 'DisplayName', 'Actual Synchrony');

xlabel('Call Number');
ylabel('Phase (deg)');
ylim([min([rad2deg(phi_star), rad2deg(wingbeat_phase),rad2deg(theta)]), max([rad2deg(phi_star), rad2deg(wingbeat_phase),rad2deg(theta)]) * 1.1]);

yyaxis right
plot(call_idx, call_rate, 'm-', 'LineWidth', 1.2, 'DisplayName', 'Call Rate (Hz)');
ylabel('Call Rate (Hz)');
ylim([0, max(call_rate) * 1.1]);

title('Synchrony Limits in Wingbeat–Call Timing During Echolocation');
subtitle('')
legend({'$\phi^*$ (Expected Call Phase)', ...
        'Current Wing phase', ...
        'Wingbeat Phase at Call Time', ...
        '$\theta$ (Wing Excursion Limit)', ...
        'Pred. Asynchronous', ...
        'Synchronous Calls', ...
        'Call Rate (Hz)'}, ...
        'Location', 'best', 'Interpreter', 'latex', 'FontSize', 16, 'Fontweight', 'bold');

formatLatexYY(gca)

%%
% Inputs from simulation result
wingbeat_frequency = result.wingbeat_frequency; % [Hz] per call
wingbeat_amplitude = result.wingbeat_amplitude; % [rad] per call
delta_t = result.delta_t; % time between calls
call_times = cumsum(delta_t); % absolute time of each call
N = length(call_times);

% Generate time vector (for smooth waveform)
T_end = call_times(end);
fs_plot = 5000; % plot sample rate (Hz)
t = linspace(0, T_end, round(T_end*fs_plot));

% For each call, get the corresponding wingbeat frequency and amplitude
% If using dynamic freq/ampl, interpolate to full time vector
% Ensure vectors are the same length for interp1
call_times_interp = call_times(1:length(wingbeat_frequency)); % Drop the last if needed
wb_freq_interp = interp1(call_times_interp, wingbeat_frequency, t, 'previous', 'extrap');
wb_amp_interp  = interp1(call_times_interp, wingbeat_amplitude, t, 'previous', 'extrap');

% Simulate a generic wingbeat as a sine wave (0 = downstroke, pi/2 = upstroke)
% If you want actual amplitude, multiply by wb_amp_interp
wingbeat_waveform = wb_amp_interp .* sin(2*pi*wb_freq_interp.*t);

% Plotting
figure('Position',[100 100 1200 400]);
plot(t, rad2deg(wingbeat_waveform), '-', 'Color', '#1E90FF', 'LineWidth', 1.5); hold on;

xlabel('Time (s)');
ylabel('Wingbeat Position (deg)');
% title('Wingbeat Cycle with Call Emissions Marked');
% legend({'Wingbeat','Call Emission'},'Location','best');

% Optionally, shade synchrony vs asynchrony using result.synchrony_flag:
sync_flag = result.actual_synchrony_flag; % logical vector for synchrony at each call
for i = 2:N-1
    [~, idx] = min(abs(t - call_times(i)));
    if sync_flag(i)
        plot(t(idx), rad2deg(wingbeat_waveform(idx)), 'o', 'Color', '#006400', 'MarkerFaceColor', '#006400', 'MarkerSize', 10, 'LineWidth', 2); % Green for sync
    else
        plot(t(idx), rad2deg(wingbeat_waveform(idx)), 'o', 'Color', '#FF4500', 'MarkerFaceColor', '#FF4500', 'MarkerSize', 8, 'LineWidth', 1);  % Black for async
    end
end
legend({'Wingbeat','Call Emission (Async)','Call Emission (Sync)'},'Location','best', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');

formatLatex(gca)


