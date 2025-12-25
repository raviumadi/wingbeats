%% analyse_simulation_results.m
% Proper analysis for Monte-Carlo simulation outputs from results table.
% Focus: descriptive stats + effect sizes + robustness, not NHST as "proof".
% Optional permutation tests are provided as model-comparison diagnostics.

clear; clc;
saveFigs = 1;
%% Load results
% Either load from .mat, or assume 'results' is in workspace
if ~exist('results','var')
    S = load('simulation_results.mat','results');
    results = S.results;
end

assert(istable(results), 'results must be a MATLAB table');

% Sanity checks for expected columns
reqCols = {'cond_id','run_id','total_calls','n_sync','n_async','async_phases','calls_per_wb'};
for k = 1:numel(reqCols)
    assert(ismember(reqCols{k}, results.Properties.VariableNames), ...
        'Missing column: %s', reqCols{k});
end

%% Add derived variables (normalised + robust)
results.frac_async = results.n_async ./ max(results.total_calls, 1); % avoid divide by zero
results.frac_sync  = results.n_sync  ./ max(results.total_calls, 1);

% Some runs might produce NaN/Inf if wingbeats_total became 0; guard it.
bad_cpw = ~isfinite(results.calls_per_wb) | results.calls_per_wb < 0;
if any(bad_cpw)
    warning('Found %d invalid calls_per_wb values; removing those rows from analyses.', sum(bad_cpw));
end

results_valid = results(~bad_cpw, :);

%% Condition labels (edit to match manuscript)
condNames = containers.Map( ...
    {1,2,3}, ...
    {'$f_w$ fixed, $\theta$ fixed', ...
     '$f_w$ dyn, $\theta$ fixed', ...
     '$f_w$ dyn, $\theta$ dyn'} );

condIDs = unique(results_valid.cond_id(:))';
condIDs = sort(condIDs);

%% ========= 1) DESCRIPTIVE SUMMARY TABLES =========
% Robust summaries per condition for key outputs

metrics = { ...
    'n_sync',      'Synchronous calls (count)'; ...
    'n_async',     'Asynchronous calls (count)'; ...
    'frac_async',  'Asynchronous calls (fraction)'; ...
    'calls_per_wb','Calls per wingbeat'};

summaryTables = struct();

for m = 1:size(metrics,1)
    varName = metrics{m,1};
    pretty  = metrics{m,2};

    Tsum = table();
    for ci = 1:numel(condIDs)
        cid = condIDs(ci);
        x = results_valid.(varName)(results_valid.cond_id == cid);
        x = x(isfinite(x));

        row = table();
        row.cond_id = cid;
        row.cond_label = string(condNames(cid));
        row.N = numel(x);

        row.mean = mean(x);
        row.std  = std(x);

        row.median = median(x);
        row.IQR = iqr(x);
        row.p05 = prctile(x, 5);
        row.p25 = prctile(x, 25);
        row.p75 = prctile(x, 75);
        row.p95 = prctile(x, 95);

        % Robust measure of spread (MAD)
        row.MAD = mad(x, 1); % median absolute deviation

        Tsum = [Tsum; row]; %#ok<AGROW>
    end

    summaryTables.(varName) = Tsum;
    fprintf('\n=== %s ===\n', pretty);
    disp(Tsum);
end

%% ========= 2) EFFECT SIZES + BOOTSTRAP CIs =========
% Pairwise comparisons: differences in medians + Cliff's delta
% (No p-values required; these quantify "how different" conditions are.)

pairList = nchoosek(condIDs,2);

effectReport = table();
B = 5000; % bootstrap resamples

for m = 1:size(metrics,1)
    varName = metrics{m,1};
    pretty  = metrics{m,2};

    for p = 1:size(pairList,1)
        a = pairList(p,1);
        b = pairList(p,2);

        xa = results_valid.(varName)(results_valid.cond_id == a);
        xb = results_valid.(varName)(results_valid.cond_id == b);
        xa = xa(isfinite(xa)); xb = xb(isfinite(xb));

        % Difference in medians (b - a)
        dMed = median(xb) - median(xa);
        ciMed = bootstrap_ci_diffmedian(xa, xb, B, 0.05);

        % Cliff's delta (robust nonparametric effect size)
        delta = cliffs_delta(xb, xa); % "b relative to a"
        ciDelta = bootstrap_ci_cliffsdelta(xa, xb, B, 0.05);

        row = table();
        row.metric = string(pretty);
        row.var = string(varName);
        row.condA = a;
        row.condB = b;
        row.labelA = string(condNames(a));
        row.labelB = string(condNames(b));
        row.diffMedian = dMed;
        row.diffMedian_CI_low = ciMed(1);
        row.diffMedian_CI_high = ciMed(2);
        row.cliffsDelta = delta;
        row.cliffsDelta_CI_low = ciDelta(1);
        row.cliffsDelta_CI_high = ciDelta(2);

        effectReport = [effectReport; row]; %#ok<AGROW>
    end
end

fprintf('\n=== Effect sizes (pairwise) ===\n');
disp(effectReport);

%% ========= 3) OPTIONAL: PERMUTATION TESTS (MODEL COMPARISON) =========
% These do NOT establish biological significance; they only quantify whether
% distributions differ under the SAME parameter ensemble.
doPermutation = false;  % set true if you want these (recommend: Supplement only)
nPerm = 10000;

permReport = table();
if doPermutation
    for m = 1:size(metrics,1)
        varName = metrics{m,1};
        pretty  = metrics{m,2};

        for p = 1:size(pairList,1)
            a = pairList(p,1);
            b = pairList(p,2);

            xa = results_valid.(varName)(results_valid.cond_id == a);
            xb = results_valid.(varName)(results_valid.cond_id == b);
            xa = xa(isfinite(xa)); xb = xb(isfinite(xb));

            % Test statistic: difference in medians
            obs = median(xb) - median(xa);
            pval = permutation_pvalue_diffmedian(xa, xb, nPerm);

            row = table(string(pretty), string(varName), a, b, obs, pval, ...
                'VariableNames', {'metric','var','condA','condB','obsDiffMedian','p_perm'});
            permReport = [permReport; row]; %#ok<AGROW>
        end
    end

    fprintf('\n=== Permutation diagnostics (optional) ===\n');
    disp(permReport);
end

%% ========= 4) CIRCULAR STATS FOR ASYNC PHASES =========
% You stored async phases as fractions of a cycle in [0,1].
% We'll compute per-condition:
% - circular mean angle (rad)
% - mean phase (fraction)
% - resultant length R (0..1): concentration
% - circular std (approx)
circSummary = table();

for ci = 1:numel(condIDs)
    cid = condIDs(ci);

    % Extract cell array of vectors
    tmp = results_valid.async_phases(results_valid.cond_id == cid);

    % Your results table stores vectors directly, not necessarily cells.
    % If it's a cell column, it will be cell; if numeric arrays, handle both.
    phases = [];
    if iscell(tmp)
        tmp_nonempty = tmp(~cellfun(@isempty, tmp));
        if ~isempty(tmp_nonempty)
            phases = [tmp_nonempty{:}];
        end
    else
        % If it is a nested type, adapt here; default:
        phases = tmp(:);
    end

    phases = phases(isfinite(phases));
    phases = phases(phases >= 0 & phases <= 1);

    row = table();
    row.cond_id = cid;
    row.cond_label = string(condNames(cid));
    row.N_async_phases = numel(phases);

    if isempty(phases)
        row.meanPhase_frac = NaN;
        row.meanAngle_rad  = NaN;
        row.R = NaN;
        row.circStd_rad = NaN;
    else
        ang = phases * 2*pi; % radians
        [mu, R] = circ_mean_and_R(ang);
        row.meanAngle_rad = mu;
        row.meanPhase_frac = mod(mu, 2*pi) / (2*pi);
        row.R = R;
        row.circStd_rad = sqrt(-2*log(max(R, eps))); % common approximation
    end

    circSummary = [circSummary; row]; %#ok<AGROW>
end

fprintf('\n=== Circular summaries of async phases ===\n');
disp(circSummary);

%% ========= 5) PLOTS (JTB-appropriate, no significance stars) =========
hgt = 500;
len = 400;
% 5.1 Box + jitter for fraction async
figure('Position',[100 100 hgt len]); hold on;
title('Fraction of asynchronous calls per simulation run', 'Interpreter','latex');
xlabel('Condition'); ylabel('Fraction async');

xpos = [];
yval = [];
grp  = [];
for ci = 1:numel(condIDs)
    cid = condIDs(ci);
    x = results_valid.frac_async(results_valid.cond_id == cid);
    x = x(isfinite(x));
    xpos = [xpos; ci*ones(size(x))]; %#ok<AGROW>
    yval = [yval; x]; %#ok<AGROW>
    grp  = [grp; cid*ones(size(x))]; %#ok<AGROW>
end

boxplot(yval, xpos, 'Labels', values_for_labels(condIDs, condNames), 'Symbol', '');
% Tick labels
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;
% ax.FontName = 'Times';

% Axis labels (separate objects!)
xlabel('Condition','Interpreter','latex','FontSize',14);
ylabel('Calls per wingbeat','Interpreter','latex','FontSize',14);

% Optional: rotate tick labels cleanly
ax.XTickLabelRotation = 0;   % or 15–30 if long labels

% Box aesthetics (optional)
set(findobj(ax,'Tag','Box'),'LineWidth',1.2);
set(findobj(ax,'Tag','Median'),'LineWidth',1.5);
% jitter points
scatter(xpos + 0.08*(rand(size(xpos))-0.5), yval, 12, 'filled', 'MarkerFaceAlpha',0.25);
grid on;
xlim([0.5 3.5])
formatLatex(gca);

% 5.2 ECDF overlay for calls_per_wb
figure('Position',[100 100 hgt len]); hold on;
title('ECDF of calls per wingbeat');
xlabel('Calls per wingbeat'); ylabel('F(x)');
for ci = 1:numel(condIDs)
    cid = condIDs(ci);
    x = results_valid.calls_per_wb(results_valid.cond_id == cid);
    x = x(isfinite(x));
    [f,xx] = ecdf(x);
    plot(xx, f, 'LineWidth', 2);
end
legend(values_for_labels(condIDs, condNames), 'Location','best', 'Interpreter','latex');
grid on;
formatLatex(gca);

% 5.3 Phase density for async phases
figure('Position',[100 100 hgt len]);
tiledlayout(numel(condIDs),1,'TileSpacing','compact','Padding','compact');
for ci = 1:numel(condIDs)
    cid = condIDs(ci);
    nexttile; hold on;

    tmp = results_valid.async_phases(results_valid.cond_id == cid);
    phases = [];
    if iscell(tmp)
        tmp_nonempty = tmp(~cellfun(@isempty, tmp));
        if ~isempty(tmp_nonempty)
            phases = [tmp_nonempty{:}];
        end
    else
        phases = tmp(:);
    end
    phases = phases(isfinite(phases));
    phases = phases(phases >= 0 & phases <= 1);

    if isempty(phases)
        text(0.5, 0.5, 'No async phases', 'HorizontalAlignment','center');
        xlim([0 1]); ylim([0 1]);
    else
        edges = linspace(0,1,25);
        h = histogram(phases, edges, 'Normalization','probability');
        ylabel('Proportion');
        xlim([0 1]);
    end
    title(string(condNames(cid)));
    if ci == numel(condIDs)
        xlabel('Wingbeat phase (fraction of cycle)');
    else
        set(gca,'XTickLabel',[]);
    end
    grid on;
    formatLatex(gca);
end

% Figure: Synchronous vs Asynchronous call counts per run
figure('Position',[100 100 hgt len]);

% Short symbolic labels for x-axis
condNamesSymbol = containers.Map( ...
    {1,2,3}, {'C1','C2','C3'} );

labels = values_for_labels(condIDs, condNamesSymbol);

% --- Subplot 1: Synchronous calls ---
subplot(1,2,1); hold on;

xpos = [];
yval = [];
for ci = 1:numel(condIDs)
    cid = condIDs(ci);
    x = results_valid.n_sync(results_valid.cond_id == cid);
    xpos = [xpos; ci*ones(size(x))]; %#ok<AGROW>
    yval = [yval; x]; %#ok<AGROW>
end

boxplot(yval, xpos, 'Labels', labels, 'Symbol','');
scatter(xpos + 0.08*(rand(size(xpos))-0.5), yval, ...
    10, 'filled', 'MarkerFaceAlpha',0.25);

ylabel('Synchronous calls per run','Interpreter','latex');
title('$n_{\mathrm{sync}}$','Interpreter','latex');
ylim([0 120])
formatLatex(gca);

% --- Subplot 2: Asynchronous calls ---
subplot(1,2,2); hold on;

xpos = [];
yval = [];
for ci = 1:numel(condIDs)
    cid = condIDs(ci);
    x = results_valid.n_async(results_valid.cond_id == cid);
    xpos = [xpos; ci*ones(size(x))]; %#ok<AGROW>
    yval = [yval; x]; %#ok<AGROW>
end

boxplot(yval, xpos, 'Labels', labels, 'Symbol','');
scatter(xpos + 0.08*(rand(size(xpos))-0.5), yval, ...
    10, 'filled', 'MarkerFaceAlpha',0.25);

ylabel('Asynchronous calls per run','Interpreter','latex');
title('$n_{\mathrm{async}}$','Interpreter','latex');
ylim([0 120])
formatLatex(gca);

% --- Add legend-style annotation explaining conditions ---
annotation('textbox',[0.18 0.68 0.35 0.18], ...
    'String',{ ...
        '\textbf{C1}: $f_w$ fixed, $\theta$ fixed', ...
        '\textbf{C2}: $f_w$ dynamic, $\theta$ fixed', ...
        '\textbf{C3}: $f_w$ dynamic, $\theta$ dynamic'}, ...
    'Interpreter','latex', ...
    'FontSize',11, ...
    'EdgeColor','none');

%% save figs
%% List of figure handles and names
if saveFigs
    figHandles = [figure(1), figure(2), figure(3), figure(4)];
    fileNames = {'frac_async_box_jitter', 'calls_per_wb_ecdf', 'async_phase_hist', 'sync_async_counts_box'};

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
end

%% ========= 6) EXPORT TABLES =========
writetable(effectReport, 'effect_sizes_report.csv');
if doPermutation
    writetable(permReport, 'permutation_report.csv');
end
writetable(circSummary, 'circular_phase_summary.csv');

disp('Saved: effect_sizes_report.csv, circular_phase_summary.csv (and permutation_report.csv if enabled).');

%% ===================== Helper functions =====================

function labels = values_for_labels(condIDs, condNames)
% Convert condIDs to a cell array of short labels for plots
labels = cell(numel(condIDs),1);
for i = 1:numel(condIDs)
    cid = condIDs(i);
    % shorten long labels if needed
    s = char(condNames(cid));
    if numel(s) > 28
        s = [s(1:26) '…'];
    end
    labels{i} = s;
end
end

function ci = bootstrap_ci_diffmedian(xa, xb, B, alpha)
% Bootstrap CI for diff in medians: median(xb)-median(xa)
if nargin < 4, alpha = 0.05; end
na = numel(xa); nb = numel(xb);
boot = zeros(B,1);
for b = 1:B
    sa = xa(randi(na, na, 1));
    sb = xb(randi(nb, nb, 1));
    boot(b) = median(sb) - median(sa);
end
ci = prctile(boot, [100*alpha/2, 100*(1-alpha/2)]);
end

function d = cliffs_delta(x, y)
% Cliff's delta: P(x>y)-P(x<y)
% Robust nonparametric effect size in [-1,1].
x = x(:); y = y(:);
nx = numel(x); ny = numel(y);
% Efficient computation using sorting
xs = sort(x);
ys = sort(y);

% Count how many y each x exceeds and is less than
% For each x: number of y < x and number of y > x
lt = zeros(nx,1);
gt = zeros(nx,1);

for i = 1:nx
    % y < x : last index where y < x
    lt(i) = find_last_lt(ys, xs(i));
    % y > x : number greater = ny - last index where y <= x
    gt(i) = ny - find_last_le(ys, xs(i));
end

d = (sum(lt) - sum(gt)) / (nx*ny);
end

function ci = bootstrap_ci_cliffsdelta(xa, xb, B, alpha)
% Bootstrap CI for Cliff's delta comparing xb vs xa
if nargin < 4, alpha = 0.05; end
na = numel(xa); nb = numel(xb);
boot = zeros(B,1);
for b = 1:B
    sa = xa(randi(na, na, 1));
    sb = xb(randi(nb, nb, 1));
    boot(b) = cliffs_delta(sb, sa);
end
ci = prctile(boot, [100*alpha/2, 100*(1-alpha/2)]);
end

function p = permutation_pvalue_diffmedian(xa, xb, nPerm)
% Permutation p-value for diff in medians (two-sided)
xa = xa(:); xb = xb(:);
obs = median(xb) - median(xa);

pool = [xa; xb];
nA = numel(xa);
n = numel(pool);

count = 0;
for k = 1:nPerm
    idx = randperm(n);
    A = pool(idx(1:nA));
    B = pool(idx(nA+1:end));
    stat = median(B) - median(A);
    if abs(stat) >= abs(obs)
        count = count + 1;
    end
end
p = (count + 1) / (nPerm + 1); % add-one smoothing
end

function [mu, R] = circ_mean_and_R(alpha)
% Circular mean and resultant length for angles alpha (rad)
alpha = alpha(:);
C = mean(cos(alpha));
S = mean(sin(alpha));
mu = atan2(S, C);
R = sqrt(C^2 + S^2);
end

function idx = find_last_lt(sortedVec, x)
% last index where sortedVec < x, or 0 if none
idx = find(sortedVec < x, 1, 'last');
if isempty(idx), idx = 0; end
end

function idx = find_last_le(sortedVec, x)
% last index where sortedVec <= x, or 0 if none
idx = find(sortedVec <= x, 1, 'last');
if isempty(idx), idx = 0; end
end