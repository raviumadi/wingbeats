clear; clc;

%% Global constants
fs = 192e3;
bandwidth = 65e3:-100:35e3;
c = 343;
n_runs = 200;  % number of repetitions per condition

%% Define simulation conditions
conditions = {
    struct('wing_sync_mode','fixed','theta_mode','fixed',  'subtitle','Constant Wingbeat Frequency, Constant Excursion')
    struct('wing_sync_mode','dynamic','theta_mode','fixed',  'subtitle','Dynamic Wingbeat Frequency, Constant Excursion')
    struct('wing_sync_mode','dynamic','theta_mode','dynamic','subtitle','Dynamic Wingbeat Frequency, Dynamic Excursion')
    };

%% Preallocate full results table with all required columns
results = table([],[],[],[],[],{},[], ...
    'VariableNames', {'cond_id','run_id','total_calls','n_sync','n_async','async_phases','calls_per_wb'});

rng('shuffle'); % for randomness

%% Main simulation loop
for cond_id = 1:length(conditions)
    condition = conditions{cond_id};
    
    for run_id = 1:n_runs
        
        % Random parameters
        params = struct();
        params.fs = fs;
        params.bandwidth = bandwidth;
        params.kr = randi([3,7]);
        params.target_distance = randi([10,15]);
        params.max_call_rate = 200;
        params.initial_velocity = randi([2,5]);
        params.initial_call_duration = 0.005;
        params.f_wing = rand_uniform(3,7);
        params.theta = rand_uniform(pi/6, pi/4);
        params.max_wingbeat_freq = rand_uniform(2,4) * params.f_wing;
        params.min_theta = rand_uniform(0.1,0.5) * params.theta;
        params.motile = 1;
        params.makeAudio = 0;
        
        % Apply condition
        params.wing_sync_mode = condition.wing_sync_mode;
        params.theta_mode = condition.theta_mode;
        
        % Run simulation
        result = simulateEcholocationWings(params);
        
        % Extract key results
        total_calls = length(result.synchrony_flag);
        n_sync = sum(result.actual_synchrony_flag == 1);
        n_async = total_calls - n_sync;
        
        % Async call phases (vectorised)
        rel_phase = mod(result.phi_star, 2*pi);
        is_async = (result.actual_synchrony_flag == 0);
        async_calls_per_wb = rel_phase(is_async) / (2*pi);
        
        % Total wingbeats
        len_dT = length(result.delta_t);
        len_fw = length(result.wingbeat_frequency);
        min_len = min(len_dT, len_fw);
        wingbeats_total = sum(result.delta_t(1:min_len) .* result.wingbeat_frequency(1:min_len));
        calls_per_wb = total_calls / wingbeats_total;
        
        % Append full row safely
        results = [results; {cond_id, run_id, total_calls, n_sync, n_async, async_calls_per_wb, calls_per_wb}];
    end
end


%% boxplot with sync & async significance markers

figure('Position', [100 100 1200 400]); hold on;

colors_sync = [0 114 189]/255;   % blue
colors_async = [217 83 25]/255;  % red

positions = [1 2; 3 4; 5 6]; % sync, async per condition
condition_labels = {'Cond 1','Cond 2','Cond 3'};

p_sync = zeros(3,3); 
p_async = zeros(3,3);
y_all = [];

for cond_id = 1:3
    data_sync = results.n_sync(results.cond_id == cond_id);
    data_async = results.n_async(results.cond_id == cond_id);
    
    % Plot sync
    boxplot(data_sync, 'Positions', positions(cond_id,1), 'Colors', colors_sync, ...
        'Widths', 0.5, 'Symbol','');
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'), get(h(1),'YData'), colors_sync, 'FaceAlpha',0.3, 'EdgeColor',colors_sync);
    
    % Plot async
    boxplot(data_async, 'Positions', positions(cond_id,2), 'Colors', colors_async, ...
        'Widths', 0.5, 'Symbol','');
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'), get(h(1),'YData'), colors_async, 'FaceAlpha',0.3, 'EdgeColor',colors_async);

    y_all = [y_all; data_sync; data_async];
end

xticks(1.5:2:6);
xticklabels(condition_labels);
xlabel('Condition','FontSize',14);
ylabel('Number of Calls','FontSize',14);
title('Synchronous (blue) and Asynchronous (red) Calls','FontSize',16);
set(gca,'FontSize',12,'XTickLabelRotation',0)
formatLatex(gca)

% Calculate stats for both sync & async
fprintf('\\nSynchronous call comparisons:\\n')
for i = 1:3
    for j = i+1:3
        [p,h] = ranksum(results.n_sync(results.cond_id==i), results.n_sync(results.cond_id==j));
        p_sync(i,j) = p;
        fprintf('Cond %d vs %d: p=%.4g\\n', i,j,p);
    end
end

fprintf('\\nAsynchronous call comparisons:\\n')
for i = 1:3
    for j = i+1:3
        [p,h] = ranksum(results.n_async(results.cond_id==i), results.n_async(results.cond_id==j));
        p_async(i,j) = p;
        fprintf('Cond %d vs %d: p=%.4g\\n', i,j,p);
    end
end

% Add significance stars for both sync and async

y_max = max(y_all);
y_start = y_max * 0.7;
y_offset = y_max * 0.05;

star_pairs = [1 3; 1 5; 3 5];
n_pairs = size(star_pairs,1);

% Sync markers
for k = 1:n_pairs
    i = star_pairs(k,1)/2 + 0.5;
    j = star_pairs(k,2)/2 + 0.5;
    p_val = p_sync(i,j);
    
    y_current = y_start + (k-1)*y_offset;
    line([star_pairs(k,1) star_pairs(k,2)], [y_current y_current],'Color','b','LineWidth',1.2);
    text(mean(star_pairs(k,:)), y_current + y_offset*0.1, significanceStars(p_val), ...
        'FontSize',16, 'HorizontalAlignment','center');
end

% Async markers (higher level)
for k = 1:n_pairs
    i = star_pairs(k,1)/2 + 0.5;
    j = star_pairs(k,2)/2 + 0.5;
    p_val = p_async(i,j);
    
    y_current = y_start + (n_pairs+k-1)*y_offset;
    line([star_pairs(k,1)+1 star_pairs(k,2)+1], [y_current y_current],'Color','r','LineWidth',1.2);
    text(mean([star_pairs(k,1)+1, star_pairs(k,2)+1]), y_current + y_offset*0.1, significanceStars(p_val), ...
        'FontSize',16, 'HorizontalAlignment','center');
end

ylim([0 y_start + (2*n_pairs)*y_offset]);

 %% Histograms of async phases
% figure('Position',[100 100 1200 400]);
% for cond_id = 1:3
%     subplot(3,1,cond_id);
%     tmp = results{results.cond_id==cond_id, 'async_phases'};
%     tmp_nonempty = tmp(~cellfun(@isempty,tmp));
%     all_async_phases = [tmp_nonempty{:}];
%     histogram(all_async_phases, 20,'FaceColor','#0072BD','FaceAlpha',0.7);
%     xlabel('Wingbeat Phase (fraction of cycle)'); ylabel('Count');
%     title(conditions{cond_id}.subtitle);
%     xlim([0 1]); grid on;
% end


%% Normalised histograms of calls-per-wingbeat

figure('Position', [100 100 1200 500]);

% Define custom colors per condition
hist_colors = [ ...
    0 114 189;    % blue
    217 83 25;    % red
    119 172 48] / 255;  % greenish

% Collect data to unify bin edges and limits
all_data = [];
for cond_id = 1:3
    all_data = [all_data; results.calls_per_wb(results.cond_id == cond_id)];
end

% Define global bin edges based on pooled data
bin_edges = linspace(min(all_data), max(all_data), 21);

for cond_id = 1:3
    subplot(3,1,cond_id);
    data = results.calls_per_wb(results.cond_id == cond_id);
    
    histogram(data, bin_edges, ...
        'Normalization','probability', ...
        'FaceColor', hist_colors(cond_id,:), ...
        'FaceAlpha', 0.6, 'EdgeColor','k','LineWidth',0.8);
    
    if cond_id == 3
    xlabel('Calls per Wingbeat','FontSize',12);
    end
    ylabel('Proportion','FontSize',12);
    title(conditions{cond_id}.subtitle,'FontSize',13);
    formatLatex(gca)
end

%% Save full results
save('simulation_results.mat','results');

%% save figs
%% List of figure handles and names

figHandles = [figure(1), figure(2)];
fileNames = {'sim_run_boxplots', 'sim_run_hist'};

savePath = '/Users/ravi/Documents/projects/wingbeat/wingbeat_call_synchrony/fig';

if ~exist(savePath, 'dir')
    mkdir(savePath); % create directory if it doesn't exist
end

for k = 1:length(figHandles)
    fig = figHandles(k);
    
    % Make sure figure is active
    figure(fig);
    
    % Save as .fig
    savefig(fig, fullfile(savePath, [fileNames{k} '.fig']));
    
    % Save as .pdf (high-quality vector)
    exportgraphics(fig, fullfile(savePath, [fileNames{k} '.pdf']), 'ContentType','vector');
end

%% Random number helper
function val = rand_uniform(a,b)
val = a + (b-a)*rand;
end

%% Significance star helper function
function stars = significanceStars(p)
    if p < 0.0001
        stars = '****';
    elseif p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    elseif p < 0.05
        stars = '*';
    else
        stars = 'n.s.';
    end
end