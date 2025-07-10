function result = simulateEcholocationWings(params)
% Simulates FM bat's echolocation behaviour with capped call rate,
% and with options for dynamic or fixed wingbeat and theta.

fs = 192e3;
c = 343;
alfa = atmatt(25, 50, params.bandwidth);
A = alfa(params.bandwidth == max(params.bandwidth));
target_diameter = freq2wavelen(min(params.bandwidth), c);
source_level = 94;
source_level_minimum = 70;
detection_level = 20;

if ~isfield(params, 'max_call_rate')
    params.max_call_rate = 200; % default in Hz
end
if ~isfield(params, 'max_wingbeat_freq')
    params.max_wingbeat_freq = 12; % default Hz
end
if ~isfield(params, 'min_theta')
    params.min_theta = 0.1; % min excursion
end
if ~isfield(params, 'wing_sync_mode')
    params.wing_sync_mode = 'dynamic';
end
if ~isfield(params, 'theta_mode')
    params.theta_mode = 'fixed';
end

new_target_distance(1) = params.target_distance;
Ta(1) = 2 * new_target_distance(1) / c;
delta_t(1) = 0;
echo_time(1) = delta_t(1) + Ta(1);

call_number = 1;
done = false;

while ~done
    Ta(call_number + 1) = 2 * new_target_distance(call_number) / c;
    % Responsivity model: get IPI
    ipi = params.kr * Ta(call_number + 1);

    % Cap IPI to max_call_rate
    min_ipi = 1 / params.max_call_rate;
    is_call_rate_capped = false;
    if ipi < min_ipi
        ipi = min_ipi;
        is_call_rate_capped = true;
    end
    delta_t(call_number + 1) = ipi;

    echo_time(call_number + 1) = max(cumsum(delta_t(1:call_number + 1))) + Ta(call_number + 1);

    if params.motile
        delta_s(call_number) = delta_t(call_number + 1) * params.initial_velocity + randn() / 10;
    else
        delta_s(call_number) = delta_t(call_number + 1) * params.initial_velocity;
    end
    v(call_number) = delta_s(call_number) / delta_t(call_number + 1);
    new_target_distance(call_number + 1) = new_target_distance(call_number) - delta_s(call_number);

    TS(call_number) = calculateTargetStrength(min(params.bandwidth), target_diameter, new_target_distance(call_number + 1), A);
    TL(call_number) = 2 * calculateGeometricAttenuation(new_target_distance(call_number + 1));
    RL(call_number) = source_level - abs(TL(call_number)) + TS(call_number);

    if RL(call_number) < detection_level
        detection_point = call_number + 1;
    end
    if Ta(call_number) > params.kr * params.initial_call_duration
        decrement_point = call_number + 1;
    end
    if new_target_distance(call_number) < 0.05
        done = true;
    end

    % --- WINGBEAT FREQUENCY/THETA SWITCH -----
    switch lower(params.wing_sync_mode)
        case {'dynamic', 'sync'}
            if call_number > 1
                call_interval = delta_t(call_number); % IPI
                f_w = 1 / call_interval;
            else
                f_w = params.f_wing;
            end
            theta = params.theta;
            % Determine the max freq limit for this mode
            if isfield(params, 'theta_mode') && strcmpi(params.theta_mode, 'dynamic')
                dynamic_max_wingbeat_freq = 2 * params.max_wingbeat_freq;
            else
                dynamic_max_wingbeat_freq = params.max_wingbeat_freq;
            end

            if f_w > dynamic_max_wingbeat_freq
                f_w = dynamic_max_wingbeat_freq; % frequency is clamped
            end

            % Now theta adjustment
            switch lower(params.theta_mode)
                case 'dynamic'
                    if call_number > 1
                        f_w_for_theta = min(f_w, dynamic_max_wingbeat_freq);
                        % theta scales inversely with frequency up to the max
                        theta = params.theta * (params.max_wingbeat_freq / f_w_for_theta);
                        % clamp theta to never be larger than the original value
                        theta = min(theta, params.theta);
                    end
                case 'fixed'
                    theta = params.theta;
                otherwise
                    error('Unknown theta_mode. Use ''dynamic'' or ''fixed''.');
            end
        case 'fixed'
            f_w = params.f_wing;
            theta = params.theta;
        otherwise
            error('Unknown wing_sync_mode. Use ''dynamic'' or ''fixed''.');
    end
    wingbeat_frequency(call_number) = f_w;
    wingbeat_amplitude(call_number) = theta;

    lambda(call_number) = v(call_number) / f_w;
    Td = params.initial_call_duration;
    Tp = 0.001;
    phi_star(call_number) = theta * f_w * (params.kr * Ta(call_number) - Tp - Td);
    synchrony_flag(call_number) = phi_star(call_number) <= theta;
    T_wing = 1 / f_w;
    call_time = sum(delta_t(1:call_number + 1));
    % phase within cycle, between 0 and 2*pi:
    phi_frac = mod(call_time, T_wing) / T_wing; % range [0,1)
    % map to radians
    phi_radians = 2 * pi * phi_frac;
    % Map to physical angular position (from -theta to +theta)
    wingbeat_phase(call_number) = theta * sin(phi_radians);

    % Actual synchrony logic
    T_wing = 1 / f_w;
    this_call_time = sum(delta_t(1:call_number + 1));
    wing_number = floor(this_call_time / T_wing);

    if call_number > 1
        prev_call_time = sum(delta_t(1:call_number));
        prev_wing_number = floor(prev_call_time / T_wing);
    else
        prev_wing_number = -1;
    end
    % Actual synchrony only if not capped at min_theta and not forced async by call rate
    if wing_number ~= prev_wing_number && phi_star(call_number) <= theta && ...
            theta > params.min_theta + eps && ~is_call_rate_capped
        actual_synchrony_flag(call_number) = 1;
    else
        actual_synchrony_flag(call_number) = 0;
    end

    call_number = call_number + 1;
end

call_rate_responsivity = abs(1 ./ diff(delta_t));
call_rate = 1./delta_t(2:end);
new_call_duration = params.initial_call_duration * ones(call_number, 1)';
new_call_duration(decrement_point:end) = Ta(decrement_point:end) ./ params.kr;
new_call_duration(new_call_duration < 0.0005) = 0.0005;

new_amplitude = source_level * ones(call_number - 1, 1);
new_amplitude(detection_point:end) = source_level - abs((detection_level - RL(detection_point:end)));
new_amplitude(new_amplitude < source_level_minimum) = linspace(source_level_minimum, source_level_minimum - 6, length(new_amplitude(new_amplitude < source_level_minimum)));
call_levels = 20e-6 * 10.^(new_amplitude ./ 20);

seq = randn(fs * round(max(cumsum(delta_t))), 1) ./ 10e5;
call_points = round(fs .* cumsum(delta_t));
echo_seq = randn(round(fs * max(echo_time)), 1) ./ 10e5;
echo_points = round(fs .* echo_time);

overlap_point = find(call_points + round(fs * new_call_duration) > echo_points, 1, 'first');
max_call_rate_point = find(1./diff(delta_t) > params.max_call_rate, 1, 'first');
[~, Tb_idx] = min(abs(call_rate_responsivity - params.max_call_rate));
Tb = delta_t(Tb_idx);
ddt = diff(delta_t);
Tb_prime = ddt(Tb_idx);
time_to_contact = new_target_distance / params.initial_velocity;

if params.makeAudio
    for i = 1:call_number - 1
        call = generateVirtualBatCall(min(params.bandwidth), max(params.bandwidth), new_call_duration(i), fs, 0);
        call = call .* call_levels(i);
        [ir, ~] = airAttenuationFilter(new_target_distance(i), fs, 64, 2);
        attenuation = 2 * calculateGeometricAttenuation(new_target_distance(i));
        att_factor = call_levels(i) - 20e-6 * 10^(abs(attenuation) / 20);
        echo = att_factor .* conv(call, ir, "same");
        seq(call_points(i)+1:call_points(i)+length(call)) = call;
        echo_seq(echo_points(i)+1:echo_points(i)+length(echo)) = echo;
    end
    audio_out = [seq(1:min(end, length(echo_seq))), echo_seq(1:min(end, length(seq)))];
else
    audio_out = [];
end

result = struct(...
    'delta_s', delta_s, ...
    'delta_t', delta_t, ...
    'call_rate_responsivity', call_rate_responsivity, ...
    'call_rate', call_rate, ...
    'call_levels', call_levels, ...
    'call_points', call_points, ...
    'decrement_point', decrement_point, ...
    'detection_point', detection_point, ...
    'echo_points', echo_points, ...
    'echo_time', echo_time, ...
    'initial_velocity', params.initial_velocity,...
    'new_amplitude', new_amplitude, ...
    'new_call_duration', new_call_duration, ...
    'initial_target_distance', params.target_distance, ...
    'new_target_distance', new_target_distance, ...
    'Ta', Ta, ...
    'Tb', Tb, ...
    'Tb_prime', Tb_prime, ...
    'Tb_idx', Tb_idx, ...
    'time_to_contact', time_to_contact, ...
    'TS', TS, ...
    'Kr', params.kr,...
    'max_call_rate', params.max_call_rate, ...
    'MaxCallRatePoint', max_call_rate_point,...
    'OverlapPoint', overlap_point,...
    'v', v, ...
    'lambda', lambda, ...
    'phi_star', phi_star, ...
    'synchrony_flag', synchrony_flag, ...
    'wingbeat_phase', wingbeat_phase, ...
    'wingbeat_frequency', wingbeat_frequency, ...
    'wingbeat_amplitude', wingbeat_amplitude, ...
    'actual_synchrony_flag', actual_synchrony_flag, ...
    'audio', audio_out ...
    );
end