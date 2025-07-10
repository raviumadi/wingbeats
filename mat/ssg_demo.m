% SSG Demo Plot
data = load('ssg_demo.mat');
result_motile = data.result;
data = load('ssg_demo_nomotile.mat');
result_nomotile = data.result;
fs = 192e3;
%% make plot
figure('Position', [100, 100, 1200, 540]);

%--- Panel placement params (all in normalized units) ---
w1 = 0.65;      % Width of waveform plot
w2 = 0.29;      % Width of IPI plot
h  = 0.35;      % Height of each row
gH = 0.1;      % Vertical gap between rows
gW = 0.01;      % Horizontal gap between columns
x0 = 0.08;      % Left margin
y1 = 0.55;      % Bottom of first row
y2 = y1 - h - gH; % Bottom of second row

%--- Row 1: Motile target ---
ax1 = axes('Position', [x0, y1, w1, h]);
plot((1:length(result_motile.audio))./fs, result_motile.audio(:,1), 'k', 'LineWidth', 0.75)
hold on
plot((1:length(result_motile.audio))./fs, 10*result_motile.audio(:,2), 'm', 'LineWidth', 0.75)
% xlabel('Time (s)')
ylabel('Amplitude')
% title('Motile Target: Call Sequence','FontWeight','bold')
xlim([6.5 9.5])
formatLatex(gca)
set(gca,'FontSize',13)
xline([6.75, 7.08], 'Color', 'k', 'Linewidth', 2, 'LineStyle', '--');
xline([7.6, 7.8], 'Color', 'k', 'Linewidth', 2, 'LineStyle', '--');
xline([8.33, 8.42], 'Color', 'k', 'Linewidth', 2, 'LineStyle', '--');

ax2 = axes('Position', [x0+w1+gW, y1, w2, h]);
plot(result_motile.call_rate, result_motile.delta_t(2:end), 'k*', 'MarkerSize', 5)
% xlabel('Call Rate (Hz)')
ylabel('IPI (s)')
title('IPI vs. Call Rate')
axis square
set(gca,'FontSize',13)
formatLatex(gca)
grid on

%--- Row 2: Stationary target ---
ax3 = axes('Position', [x0, y2, w1, h]);
plot((1:length(result_nomotile.audio))./fs, result_nomotile.audio(:,1), 'b', 'LineWidth', 0.75)
hold on
plot((1:length(result_nomotile.audio))./fs, 3*result_nomotile.audio(:,2), 'm', 'LineWidth', 0.75)
xlabel('Time (s)')
ylabel('Amplitude')
% title('Stationary Target: Call Sequence','FontWeight','bold')
xlim([6.5 9.5])
ylim([-0.2 0.2])
set(gca,'FontSize',13)
formatLatex(gca)

ax4 = axes('Position', [x0+w1+gW, y2, w2, h]);
plot(result_nomotile.call_rate, result_nomotile.delta_t(2:end), 'b*', 'MarkerSize', 5)
xlabel('Call Rate (Hz)')
ylabel('IPI (s)')
% title('IPI vs. Call Rate')
axis square
set(gca,'FontSize',13)
formatLatex(gca)
grid on

%--- Figure title ---
annotation('textbox', [0.08, 0.96, 0.84, 0.05], ...
    'String', 'Echolocation Simulation: Motile vs. Stationary Target', ...
    'EdgeColor', 'none', 'FontWeight', 'bold', 'FontSize', 18, ...
    'HorizontalAlignment', 'center', 'Interpreter', 'latex');

%%
exportgraphics(gcf, '/Users/ravi/Documents/projects/wingbeat/wingbeat_call_synchrony/fig/ssg_demo.pdf', 'Resolution', 300)