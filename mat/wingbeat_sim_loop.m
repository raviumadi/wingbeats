fs = 192e3;
bandwidth = 65e3:-100:35e3;
c = 343;

% Common parameter structure
params = struct();
params.fs = fs;
params.bandwidth = bandwidth;
params.kr = 5;
params.target_distance = 9;  % m
params.max_call_rate = 200;
params.initial_velocity = 2;  % m/s
params.initial_call_duration = 0.005;  % s
params.f_wing = 5;  % Hz
params.theta = pi / 4;
params.max_wingbeat_freq = 3 * params.f_wing;
params.min_theta = 0.3 * params.theta;
params.motile = 0;
params.makeAudio = 0;

% Simulation conditions
conditions = {
    struct('wing_sync_mode','fixed','theta_mode','fixed',  'subtitle','Constant Wingbeat Frequency, Constant Excursion ($f_w$ fixed, $\theta$ fixed)')
    struct('wing_sync_mode','dynamic','theta_mode','fixed',  'subtitle','Dynamic Wingbeat Frequency, Constant Excursion ($f_w$ dynamic, $\theta$ fixed)')
    struct('wing_sync_mode','dynamic','theta_mode','dynamic','subtitle','Dynamic Wingbeat Frequency, Dynamic Excursion ($f_w$ dynamic, $\theta$ dynamic)')
    };

for i = 1:length(conditions)
    % Apply condition-specific parameters
    params.wing_sync_mode = conditions{i}.wing_sync_mode;
    params.theta_mode     = conditions{i}.theta_mode;
    subtitle_str          = conditions{i}.subtitle;

    % Run simulation
    result = simulateEcholocationWings(params);

    % Gather results for plotting
    phi_star = result.phi_star;
    wingbeat_phase = result.wingbeat_phase;
    synchrony_flag = result.synchrony_flag;
    actual_synchrony_flag = result.actual_synchrony_flag;
    call_rate = result.call_rate;
    call_idx = 1:length(phi_star);

    % Use initial theta for yline
    theta = params.theta;

    % Plot
    figure('Position', [100, 100, 1200, 400]);
    yyaxis left
    plot(call_idx, rad2deg(phi_star), '-', 'Color', [0,128,128]./255, 'LineWidth', 2); hold on;
    plot(call_idx, rad2deg(wingbeat_phase), 'k--', 'LineWidth', 2);
    plot(call_idx, rad2deg(wingbeat_phase), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    yline(rad2deg(theta), 'k:', 'Max Wing Excursion $\theta$', 'LineWidth', 2,'Interpreter', 'latex', 'FontSize', 16, 'Fontweight', 'bold');

    % Show predicted asynchronous calls
    async_pred_idx = find(synchrony_flag == 0);
    scatter(call_idx(async_pred_idx), rad2deg(phi_star(async_pred_idx)), ...
        50, [255,69,0]./255, 'filled', 'DisplayName', 'Pred. Asynchronous');

    % Show actual synchronous calls
    actual_sync_idx = find(actual_synchrony_flag == 1);
    scatter(call_idx(actual_sync_idx), rad2deg(phi_star(actual_sync_idx)), ...
        70, 'Marker', 'o', 'MarkerFaceColor', [255,215,0]./255, 'MarkerEdgeColor', [255,140,0]./255, 'LineWidth', 2, 'DisplayName', 'Synchronous Calls');

    xlabel('Call Number');
    ylabel('Phase (deg)');
    ylim([min([rad2deg(phi_star), rad2deg(wingbeat_phase),rad2deg(theta)]), max([rad2deg(phi_star), rad2deg(wingbeat_phase),rad2deg(theta)]) * 1.1]);

    yyaxis right
    plot(call_idx, call_rate, 'r-', 'LineWidth', 1.2, 'DisplayName', 'Call Rate (Hz)');
    ylabel('Call Rate (Hz)');
    ylim([0, max(call_rate) * 1.1]);

    title('Synchrony Limits in Wingbeat-Call Timing During Echolocation', 'Interpreter', 'latex', 'FontWeight','bold');
    subtitle(subtitle_str);
    if i ==1
        legend({'$\phi^*$ (Expected Call Phase)', ...
            'Current Wing phase', ...
            'Wingbeat Phase at Call Time', ...
            '$\theta$ (Wing Excursion Limit)', ...
            'Pred. Asynchronous', ...
            'Synchronous Calls', ...
            'Call Rate (Hz)'}, ...
            'Location', 'best', 'Interpreter', 'latex', 'FontSize', 16, 'Fontweight', 'bold');
    end
    formatLatexYY(gca)

    % Generate a suitable filename for each condition
    savepath = '/Users/ravi/Documents/projects/wingbeat/wingbeat_call_synchrony/fig';
    filename = sprintf('synchrony_condition_%d.pdf', i);
    fullfile_out = fullfile(savepath, filename);

    exportgraphics(gcf, fullfile_out, 'ContentType','vector', 'BackgroundColor','none', 'Resolution',300);

    % For PNG as well, uncomment this:
    % filename_png = sprintf('synchrony_condition_%d.png', i);
    % exportgraphics(gcf, fullfile(savepath, filename_png), 'Resolution', 300);

    %% waveform plots
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
    subtitle(subtitle_str);

    % Optionally, shade synchrony vs asynchrony using result.synchrony_flag:
    sync_flag = result.actual_synchrony_flag; % logical vector for synchrony at each call
    for k = 2:N-1
        [~, idx] = min(abs(t - call_times(k)));
        if sync_flag(k)
            plot(t(idx), rad2deg(wingbeat_waveform(idx)), 'o', 'Color', '#006400', 'MarkerFaceColor', '#006400', 'MarkerSize', 10, 'LineWidth', 2); % Green for sync
        else
            plot(t(idx), rad2deg(wingbeat_waveform(idx)), 'o', 'Color', '#FF4500', 'MarkerFaceColor', '#FF4500', 'MarkerSize', 8, 'LineWidth', 1);  % Black for async
        end
    end
    if i == 1
        legend({'Wingbeat','Call Emission (Async)','Call Emission (Sync)'},'Location','best', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
    end
    formatLatex(gca)

     % Generate a suitable filename for each condition
    savepath = '/Users/ravi/Documents/projects/wingbeat/wingbeat_call_synchrony/fig';
    filename = sprintf('waveform_synchrony_condition_%d.pdf', i);
    fullfile_out = fullfile(savepath, filename);

    exportgraphics(gcf, fullfile_out, 'ContentType','vector', 'BackgroundColor','none', 'Resolution',300);

    % For PNG as well, uncomment this:
    % filename_png = sprintf('synchrony_condition_%d.png', i);
    % exportgraphics(gcf, fullfile(savepath, filename_png), 'Resolution', 300);

end