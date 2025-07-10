%% Supplementary Figure: Combined Analysis

% Load dataset
data = readtable('/Users/ravi/Documents/projects/vo_field/analysis/results/vof_processed_data.csv');

% Mic geometry and positions
refMicXYZ = [0, 0, -0.93];
batXYZ = [data.X, data.Y, data.Z];
vecToMic = batXYZ - refMicXYZ;
distToMic = sqrt(sum(vecToMic.^2, 2));

call_rate = data.Rate;
ipi = data.IPI;
distance = data.Distance;
velocity = data.Velocity;
call_level = data.Level;


% Acoustic constants
ref_dist = 0.1;
freq_kHz = 45;
f = freq_kHz * 1e3;
T_C = 20;
T_K = T_C + 273.15;
h = 50;
P = 101.325;

% Atmospheric attenuation (Bass et al.)
frO = 24 + (4.04e4 * h * (0.02 + h) / ((T_K^2) * (0.391 + h)));
frN = P * (T_K / 293.15)^(-1/2) * (9 + 280 * h * exp(-4.17 * ((T_K/293.15)^(-1/3) - 1)));
alpha = 8.686 * f^2 * ...
    (1.84e-11 * P^(-1) * sqrt(T_K) + ...
    ((0.01275 * exp(-2239.1 / T_K)) / (frO + f^2 / frO)) + ...
    ((0.1068 * exp(-3352 / T_K)) / (frN + f^2 / frN)));

% Transmission loss
TL_spherical = 20 * log10(distToMic ./ ref_dist);
TL_attenuation = alpha * (distToMic - ref_dist);
TL_total = TL_spherical + TL_attenuation;

% Extract and filter
velocity = data.Velocity;
validIdx = velocity <= 30;
fRate = data.Rate(validIdx);
fDurations = data.Duration(validIdx);
fLevels = data.Level(validIdx) + TL_total(validIdx);

% Compute deltas
dCallRate = diff(fRate);
dCallLevel = diff(fLevels);

%% Create combined figure
figure('Position', [100, 100, 1500, 500]);

% --- Subplot 1: Call Rate vs Call Level ---
subplot(2,2,1);
scatter(fRate, fLevels, 10, 'k', 'filled'); hold on;
p = polyfit(fRate(:), fLevels(:), 1);
plot(fRate, polyval(p, fRate), 'r--', 'LineWidth', 1.5);
R2 = corr(fRate(:), fLevels(:))^2;
xlabel('Call Rate (Hz)', 'Interpreter', 'latex');
ylabel('Call Level (dB SPL @10\,cm)', 'Interpreter', 'latex');
title('Call Rate vs.\ Call Level', 'Interpreter', 'latex');
legend('Data', sprintf('Fit: y = %.2fx + %.2f R$^2$ = %.3f', p(1), p(2), R2), ...
       'Location', 'northeast', 'Interpreter', 'latex');
axis tight
pbaspect([1 1 1])
formatLatex(gca)

% --- Subplot 2: Call Duration vs Call Level ---
subplot(2,2,2);
scatter(fDurations, fLevels, 10, 'k', 'filled'); hold on;
p = polyfit(fDurations(:), fLevels(:), 1);
plot(fDurations, polyval(p, fDurations), 'r--', 'LineWidth', 1.5);
R2 = corr(fDurations(:), fLevels(:))^2;
xlabel('Call Duration (s)', 'Interpreter', 'latex');
ylabel('Call Level (dB SPL @10\,cm)', 'Interpreter', 'latex');
title('Call Duration vs.\ Call Level', 'Interpreter', 'latex');
legend('Data', sprintf('Fit: y = %.2fx + %.2f R$^2$ = %.3f', p(1), p(2), R2), ...
       'Location', 'northeast', 'Interpreter', 'latex');
axis square
% pbaspect([1 1 1])
formatLatex(gca)



% --- Subplot 3: Δ Call Rate vs Δ Call Level ---
subplot(2,2,3);
scatter(dCallRate, dCallLevel, 10, 'k', 'filled'); hold on;
p = polyfit(dCallRate(:), dCallLevel(:), 1);
plot(dCallRate, polyval(p, dCallRate), 'r--', 'LineWidth', 1.5);
R2 = corr(dCallRate(:), dCallLevel(:))^2;
xlabel('$\Delta$ Call Rate (Hz)', 'Interpreter', 'latex');
ylabel('$\Delta$ Call Level (dB SPL @10\,cm)', 'Interpreter', 'latex');
title('$\Delta$ Call Rate vs.\ $\Delta$ Call Level', 'Interpreter', 'latex');
legend('Data', sprintf('Fit: y = %.2fx + %.2f $R^2$ = %.3f', p(1), p(2), R2), ...
       'Location', 'best', 'Interpreter', 'latex');

formatLatex(gca)
axis square
% Microphone setup
MicrophonePositions = [ ...
     0,    0,     0;      % Mic 1
    -1.15, 0, -1.10;      % Mic 2
     0,    0, -0.93;      % Mic 3 (reference)
     1.15, 0, -1.10];     % Mic 4
refMicXYZ = MicrophonePositions(3, :);  % Reference mic at [0 0 -0.93]

% Acoustic Parameters
freq_kHz = 45;
f = freq_kHz * 1e3;
T_C = 20;
T_K = T_C + 273.15;
h = 50;  % % RH
P = 101.325;  % kPa
ref_dist = 0.1;  % m

% Atmospheric attenuation (Bass et al. 1995 approximation)
frO = 24 + (4.04e4 * h * (0.02 + h) / ((T_K^2) * (0.391 + h)));
frN = P * (T_K / 293.15)^(-1/2) * (9 + 280 * h * exp(-4.17 * ((T_K/293.15)^(-1/3) - 1)));
alpha = 8.686 * f^2 * ...
        (1.84e-11 * P^(-1) * sqrt(T_K) + ...
        ((0.01275 * exp(-2239.1 / T_K)) / (frO + f^2 / frO)) + ...
        ((0.1068 * exp(-3352 / T_K)) / (frN + f^2 / frN)));  % dB/m

% Transmission loss
vecToMic = batXYZ - refMicXYZ;
distToMic = sqrt(sum(vecToMic.^2, 2));
TL_spherical = 20 * log10(distToMic ./ ref_dist);
TL_attenuation = alpha * (distToMic - ref_dist);
TL_total = TL_spherical + TL_attenuation;

% Filter for valid velocity (<30 m/s)
validIdx = velocity <= 30;
fDists = distance(validIdx);
fLevels = call_level(validIdx) + TL_total(validIdx);  % source level at 10 cm
fRate = call_rate(validIdx);
fVel = velocity(validIdx);
fIPI = fDists ./ fVel;

% Claculate sequencwise diff
dLevels = [];
dDists = [];

% Extract sequence numbers
seqNums = data.SeqNum(validIdx);
seqLevels = fLevels;
seqDists = fDists;

% Get unique sequence IDs
uniqueSeqs = unique(seqNums);

% Loop over sequences and compute local diffs
for i = 1:length(uniqueSeqs)
    thisSeq = uniqueSeqs(i);
    idx = find(seqNums == thisSeq);
    
    if numel(idx) < 2
        continue; % skip sequences with <2 points
    end

    dLevels = [dLevels; diff(seqLevels(idx))];
    dDists  = [dDists; diff(seqDists(idx))];
end
% Hyperbolic fit (t = a / (r + b))
hyperbolic = @(b, x) b(1) ./ (x + b(2));
initial_guess = [1, 0];
fit1 = fitnlm(fRate, fIPI, hyperbolic, initial_guess);
a = fit1.Coefficients.Estimate(1);
b = fit1.Coefficients.Estimate(2);
x_fit = linspace(min(fRate), max(fRate), 500);
y_fit = a ./ (x_fit + b);

% Linear regression (ΔLevel vs. ΔDistance)
p = polyfit(dLevels, dDists, 1);
y_reg = polyval(p, dLevels);
R2 = corr(dLevels, dDists)^2;

% Subplot 2: ΔLevel vs ΔDistance
subplot(2,2,4);
scatter(dLevels, dDists, 10, 'k', 'filled'); hold on;
plot(dLevels, y_reg, 'r--', 'LineWidth', 1.5);
xlabel('$\Delta$ Source Level (dB SPL @10cm)');
ylabel('$\Delta$ Distance Covered (m)');
title('$\Delta$ Call Level vs. $\Delta$ Distance');
legend('Data', sprintf('Fit: y = %.3fx + %.3f R$^2$ = %.3f', p(1), p(2), R2), 'Location', 'best', 'Interpreter', 'latex');

formatLatex(gca)
axis square

%% export graphics
% exportgraphics(gcf, '/Users/ravi/Documents/projects/wingbeat/wingbeat_call_synchrony/fig/level_parameters.pdf', 'Resolution', 300)