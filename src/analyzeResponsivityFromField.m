function results = analyzeResponsivityFromField(callTimes, RcMax, varargin)
% analyzeResponsivityFromField - Analyze responsivity and buzz readiness from field call timestamps,
% with optional 3D position data, call duration vector, and velocity vector to estimate Ta, Tb, distance, and contraction onset.
%
% Inputs:
%   callTimes : vector of call timestamps (s)
%   RcMax     : maximum physiological call rate (Hz)
%   varargin  : optional inputs:
%               - batXYZ : Nx3 matrix of 3D bat positions (meters)
%               - t_call : Nx1 vector of call durations (seconds)
%               - vel    : Nx1 vector of measured forward velocities (m/s)
%
% Outputs:
%   results : structure containing
%       - Rc_n : instantaneous call rates
%       - R : responsivity values
%       - readinessTime : time at buzz readiness
%       - readinessIndex : index of buzz readiness
%       - Tb_prime : reaction window at buzz readiness (s)
%       - Ta : echo delay (if batXYZ provided)
%       - Tb : biological reaction time (if batXYZ provided)
%       - distance : Euclidean distance from array (if batXYZ provided)
%       - v_est : estimated or input forward velocity
%       - d_crit : call duration contraction threshold distance (if t_call provided)
%       - contractionIndex : index where t_call drops below Ta + Tb (optional)
%       - callTimes : input call times
%       - t_call : input call durations (if provided)

% Handle optional inputs
batXYZ = [];
t_call = [];
vel = [];
if ~isempty(varargin)
    for i = 1:length(varargin)
        arg = varargin{i};
        if ismatrix(arg) && size(arg,2) == 3
            batXYZ = arg;
        elseif isvector(arg) && isempty(t_call)
            t_call = arg(:); % first vector is t_call
        elseif isvector(arg) && isempty(vel)
            vel = arg(:); % second vector is vel
        end
    end
    if ~isempty(batXYZ) && size(batXYZ,1) ~= length(callTimes)
        error('batXYZ must match callTimes length');
    end
    if ~isempty(t_call)
        if length(t_call) == length(callTimes) - 1
            t_call = [t_call; NaN]; % pad to match callTimes
        elseif length(t_call) > length(callTimes)
            t_call = t_call(1:length(callTimes)); % trim excess
        elseif length(t_call) < length(callTimes)
            error('t_call must match callTimes length or be one less');
        end
    end
    if ~isempty(vel) && length(vel) ~= length(callTimes) - 1
        error('vel must have one less element than callTimes (to match IPI length)');
    end
end

% Compute IPI and call rate
ipi = diff(callTimes);
Rc_n = 1 ./ ipi;

% % Responsivity
% R = abs(1 ./ ipi(2:end) - 1 ./ ipi(1:end-1));
% Rc_n_short = Rc_n(1:end-1);
% readinessTarget = abs(R - (RcMax - Rc_n_short));
% [~, readinessIdx] = min(readinessTarget);
% readinessTime = callTimes(readinessIdx + 1);
% Tb_prime = ipi(readinessIdx);

% Responsivity: inverse of IPI change magnitude
dIpi = diff(ipi);                          % Δt_{n+1} - Δt_n
% R = abs(1 ./ dIpi);                        % Responsivity = |1 / ΔΔt|
R = 1 ./ dIpi;                        % Responsivity = |1 / ΔΔt|
R = [NaN; R];                              % pad to align with ipi index n

% Buzz readiness index: where responsivity is maximal
% [~, readinessIdx] = max(R);
[~, readinessIdx] = min(abs(abs(R) - RcMax));        % Max responsivity
readinessTime = callTimes(readinessIdx);  % +1 because diff shortens vector

% Tb*: IPI change at readiness
if readinessIdx < length(ipi)
    % Tb_prime = abs(ipi(readinessIdx+1) - ipi(readinessIdx));
    Tb_prime = dIpi(readinessIdx+1);
else
    Tb_prime = NaN;
end

% Estimate Ta, Tb, d_crit
Ta = [];
Tb = [];
d_crit = [];
contractionIndex = NaN;
v_est = [];
distance = [];

if ~isempty(batXYZ)
    c = 343;
    d = sqrt(sum(batXYZ(2:end,:).^2, 2));
    Ta = (2 * d) / c;
    kr = 1 / (RcMax * min(Ta));
    Tb = kr .* Ta;

    % Distance from array at each point
    distance = sqrt(sum(batXYZ.^2, 2));
end

% Use input velocity or estimate it
if ~isempty(vel)
    v_est = vel;
elseif ~isempty(batXYZ)
    posDiff = diff(batXYZ);
    dispVec = sqrt(sum(posDiff.^2, 2));
    v_est = dispVec ./ ipi(:);
end

% Calculate d_crit if call durations and velocities are available
if ~isempty(t_call) && ~isempty(v_est) && ~isempty(Tb)
    t_call_ipi = t_call(2:end);
    d_crit = (t_call_ipi * 343 / 2) + Tb .* v_est;

    contractionCheck = t_call_ipi < (Ta + Tb);
    contractionIndex = find(contractionCheck, 1);
end

% Output
results = struct();
results.callTimes = callTimes;
results.Rc_n = Rc_n;
results.R = R;
results.readinessTime = readinessTime;
results.readinessIndex = readinessIdx;
results.ipi = ipi;
results.Tb_prime = Tb_prime;
results.RcMax = RcMax;
results.Ta = Ta;
results.Tb = Tb;
results.distance = distance;
results.v_est = v_est;
results.d_crit = d_crit;
results.contractionIndex = contractionIndex;
results.t_call = t_call;
results.BatXYZ = batXYZ;
end
