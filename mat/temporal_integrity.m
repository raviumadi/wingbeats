%% Load and extract data
data = readtable('/Users/ravi/Documents/projects/vo_field/analysis/results/vof_processed_data.csv');
call_rate = data.Rate;
ipi = data.IPI;
distance = data.Distance;
velocity = data.Velocity;
call_level = data.Level;
batXYZ = [data.X, data.Y, data.Z];

%% Microphone setup
MicrophonePositions = [ ...
     0,    0,     0;      % Mic 1
    -1.15, 0, -1.10;      % Mic 2
     0,    0, -0.93;      % Mic 3 (reference)
     1.15, 0, -1.10];     % Mic 4
refMicXYZ = MicrophonePositions(3, :);  % Reference mic at [0 0 -0.93]

%% Acoustic Parameters
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

%% Transmission loss
vecToMic = batXYZ - refMicXYZ;
distToMic = sqrt(sum(vecToMic.^2, 2));
TL_spherical = 20 * log10(distToMic ./ ref_dist);
TL_attenuation = alpha * (distToMic - ref_dist);
TL_total = TL_spherical + TL_attenuation;

%% Filter for valid velocity (<30 m/s)
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


%% Hyperbolic fit (t = a / (r + b))
hyperbolic = @(b, x) b(1) ./ (x + b(2));
initial_guess = [1, 0];
fit1 = fitnlm(fRate, fIPI, hyperbolic, initial_guess);
a = fit1.Coefficients.Estimate(1);
b = fit1.Coefficients.Estimate(2);
x_fit = linspace(min(fRate), max(fRate), 500);
y_fit = a ./ (x_fit + b);

%% Linear regression (ΔLevel vs. ΔDistance)
p = polyfit(dLevels, dDists, 1);
y_reg = polyval(p, dLevels);
R2 = corr(dLevels, dDists)^2;

%% Plotting: IPI + ΔLevel
figure('Position', [100, 100, 1200, 500]);

% Subplot 1: Call Rate vs IPI
subplot(1,2,1);
plot(x_fit, y_fit, 'r--', 'LineWidth', 2.5); hold on;
scatter(fRate, fIPI, 10, 'k', 'filled');
xlabel('Call Rate (Hz)');
ylabel('$Interval = \frac{Distance}{Velocity}$ (s)');
title('Call Rate vs. Interpulse Interval');
legend(sprintf('Fit: IPI = %.2f / ($C_r$ + %.2f)', a, b), 'Data', 'Location', 'northeast', 'Interpreter', 'latex');
formatLatex(gca)
axis square

% % Subplot 2: ΔLevel vs ΔDistance
% subplot(1,2,2);
% scatter(dLevels, dDists, 10, 'k', 'filled'); hold on;
% plot(dLevels, y_reg, 'r--', 'LineWidth', 1.5);
% xlabel('$\Delta$ Source Level (dB SPL @10cm)');
% ylabel('$\Delta$ Distance Covered (m)');
% title('$\Delta$ Call Level vs. $\Delta$ Distance');
% legend('Data', sprintf('Fit: y = %.3fx + %.3f\nR$^2$ = %.3f', p(1), p(2), R2), 'Location', 'best', 'Interpreter', 'latex');
% formatLatex(gca)



%% -------- Section: Call Rate vs. Distance per Call Fit -------- %%
% Filter again for consistency
call_rate = call_rate(validIdx);
distance_per_call = distance(validIdx);  % d = v / r

% Fit hyperbolic model: d = a / (r + b)
fit2 = fitnlm(call_rate, distance_per_call, hyperbolic, [1, 1]);
a2 = fit2.Coefficients.Estimate(1);
b2 = fit2.Coefficients.Estimate(2);
x2 = linspace(min(call_rate), max(call_rate), 500);
y2 = a2 ./ (x2 + b2);

% Plot
subplot(1,2,2);
scatter(call_rate, distance_per_call, 10, 'k', 'filled'); hold on;
plot(x2, y2, 'r--', 'LineWidth', 1.5);
xlabel('Call Rate (Hz)');
ylabel('Distance per Call (m)');
title('Distance per Call vs. Call Rate');
legend('Data', sprintf('Fit: d = %.2f / ($C_r$ + %.2f)', a2, b2), 'Location', 'northeast', 'Interpreter', 'latex');
formatLatex(gca)
axis square
% grid on;

%% Export Figure
exportgraphics(gcf, '/Users/ravi/Documents/projects/wingbeat/wingbeat_call_synchrony/fig/temporal_integrity.pdf', 'Resolution', 300)