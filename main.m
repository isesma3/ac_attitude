%% PROJECT  B -- AE6520 -- SPRING 2025
% author: Ismael Rodriguez Sesma

clc; clear; close all
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultLineLineWidth', 2)

%% INPUTS
% Weight and geometric data (imperial units)
W = 17578;
GRAV = 32.1745;
MASS = W/GRAV;

SURF = 260;
WING_SPAN = 27.5;
CHORD = 10.8;
cogPosit = 0.25*CHORD;

ALT = 35000;
MACH = 0.6;
TRIM_AOA = 8.8;
THETA_INIT = 0;
FLIGHT_PATH = 10;   % [deg]
TURN_RATE = 2;      % [deg/s]

I_XX = 8190;
I_YY = 25900;
I_ZZ = 29100;
I_XZ = -1952;

I_C = I_XX*I_ZZ - I_XZ^2;

% Stability derivatives (all normalized by mass or mass moment of inertia)
X_U = -0.0128;
X_AOA = -20.61;
X_Q = 0;
X_DE = 2.69;
X_DT = 0;

Z_U = -0.114;
Z_AOA = -218.3;
Z_AOA_DOT = 0;
Z_Q = 0;
Z_DE = -23.73;
Z_DT = 0;

M_U = 0.0004;
M_AOA = -5.402;
M_AOA_DOT = -0.16;
M_Q = -0.484;
M_DE = -8.1;
M_DT = 0;

Y_B = -60.38;
Y_P = 0;
Y_R = 0;
Y_DA = -0.478;
Y_DR = 10.46;

L_B = -14.24;
L_P = -0.671;
L_R = 0.464;
L_DA = 7.998;
L_DR = 2.739;

N_B = 7.864;
N_P = -0.004;
N_R = -0.291;
N_DA = -0.139;
N_DR = -3.517;

% Now convert to dimensional values
% Forces are normalized with weight
X_U = MASS*X_U;
X_AOA = MASS*X_AOA;
X_Q = MASS*X_Q;
X_DE = MASS*X_DE;
X_DT = MASS*X_DT;

Z_U = MASS*Z_U;
Z_AOA = MASS*Z_AOA;
Z_AOA_DOT = MASS*Z_AOA_DOT;
Z_Q = MASS*Z_Q;
Z_DE = MASS*Z_DE;
Z_DT = MASS*Z_DT;

Y_B = MASS*Y_B;
Y_P = MASS*Y_P;
Y_R = MASS*Y_R;
Y_DA = MASS*Y_DA;
Y_DR = MASS*Y_DR;

% Moments are normalized with moments of inertia
M_U = I_YY*M_U;
M_AOA = I_YY*M_AOA;
M_AOA_DOT = I_YY*M_AOA_DOT;
M_Q = I_YY*M_Q;
M_DT = I_YY*M_DT;
M_DE = I_YY*M_DE;

L_B = I_XX*L_B;
L_P = I_XX*L_P;
L_R = I_XX*L_R;
L_DA = I_XX*L_DA;
L_DR = I_XX*L_DR;

N_B = I_ZZ*N_B;
N_P = I_ZZ*N_P;
N_R = I_ZZ*N_R;
N_DA = I_ZZ*N_DA;
N_DR = I_ZZ*N_DR;

%% EQUILIBRIUM CONDITION
[temp, soundSpeed] = atmosisa(convlength(ALT, 'ft', 'm'));
equilSpeed = convlength(MACH*soundSpeed, 'm', 'ft');

L_AOA   = - Z_AOA;
L_DE = - Z_DE;
L_Q = - Z_Q;

% Define system of equations
f = @(x) [
    L_AOA * x(1) + L_DE * x(5) + L_Q * deg2rad(TURN_RATE) * sin(x(3)) - MASS * GRAV * (cos(deg2rad(FLIGHT_PATH))/cos(x(3)) - 1);     % Eq 11.1
    M_AOA * x(1) + M_DE * x(5) + M_Q * deg2rad(TURN_RATE) * sin(x(3));                                                              % Eq 11.2
    Y_B  * x(4) - Y_P * deg2rad(TURN_RATE) * x(3) + Y_DA * x(6) + Y_DR * x(7) + Y_R * deg2rad(TURN_RATE) * cos(x(3));               % Eq 11.3
    L_B  * x(4) - L_P * deg2rad(TURN_RATE) * x(2) + L_DA * x(6) + L_DR * x(7) + L_R * deg2rad(TURN_RATE) * cos(x(3));               % Eq 11.4
    N_B  * x(4) - N_P * deg2rad(TURN_RATE) * x(2) + N_DA * x(6) + N_DR * x(7) + N_R * deg2rad(TURN_RATE) * cos(x(3));               % Eq 11.5
    x(3) - atan(equilSpeed * deg2rad(TURN_RATE) / GRAV);                                                                            % Eq 11.6
    x(2) - x(1) - deg2rad(FLIGHT_PATH);                                                                                             % Eq 11.7
];

% x = [delta_alpha, theta, phi, beta, delta_e, delta_a, delta_r]
x0 = deg2rad([0, 0, 0, 0, 0, 0, 0]);  % Initial guess

% Solve
opts = optimoptions('fsolve','FunctionTolerance',1e-10);
[x_sol, fval, exitflag] = fsolve(f, x0, opts);

% Output results
deg = @(r) rad2deg(r);
fprintf('\n--- Trimmed Climbing Turn Solution ---\n');
fprintf('Delta_alpha : %.3f deg\n', deg(x_sol(1)));
fprintf('Theta       : %.3f deg\n', deg(x_sol(2)));
fprintf('Phi         : %.3f deg\n', deg(x_sol(3)));
fprintf('Beta        : %.3f deg\n', deg(x_sol(4)));
fprintf('Delta_e     : %.3f deg\n', deg(x_sol(5)));
fprintf('Delta_a     : %.3f deg\n', deg(x_sol(6)));
fprintf('Delta_r     : %.3f deg\n', deg(x_sol(7)));

%% INERTIAL CROSS-COUPLING
% Convert angle derivatives to speed derivatives (AOA to W | Beta to V)
X_W = X_AOA/equilSpeed;
Z_W = Z_AOA/equilSpeed;
M_W = M_AOA/equilSpeed;
Z_W_DOT = Z_AOA_DOT/equilSpeed;
M_W_DOT = M_AOA_DOT/equilSpeed;

Y_V = Y_B/equilSpeed;
L_V = L_B/equilSpeed;
N_V = N_B/equilSpeed;

% Study coupling
I_X = I_XX;
I_Y = I_YY;
I_Z = I_ZZ;

ROLL_RATE = 0:0.05:4;   % [rad/s]
M_AOA_PRIME = M_AOA + M_AOA_DOT*Z_W/MASS;
M_Q_PRIME = M_Q + M_AOA_DOT;
MU_1 = (I_Z - I_X)/I_Y;
MU_2 = (I_X - I_Y)/I_Z;

nSteps = length(ROLL_RATE);
detCoup = zeros(1, nSteps);

for i = 1:nSteps
  
    coupMat = [Z_W/MASS 1 -ROLL_RATE(i) 0;
        M_AOA_PRIME/I_Y M_Q_PRIME/I_Y -M_AOA_DOT * ROLL_RATE(i)/I_Y ROLL_RATE(i) * MU_1;
        ROLL_RATE(i) 0 Y_V/MASS -1;
        0 ROLL_RATE(i) * MU_2 N_B/I_Z N_R/I_Z];
    
    detCoup(i) = det(coupMat);
end

plot(ROLL_RATE, detCoup)
grid on
xlabel('Roll rate, $P_0$ [rad/s]');
ylabel('Determinant of A matrix');
title('Inertia cross--coupling evaluation')

%% STABILITY & CONTROL MATRICES
% Longitudinal dyamics
longStab = [X_U/MASS X_W/MASS X_Q -GRAV*cos(deg2rad(THETA_INIT));

    Z_U/(MASS - Z_W_DOT) Z_W/(MASS - Z_W_DOT) (Z_Q + MASS*equilSpeed)/(MASS - Z_W_DOT) ...
    -MASS*GRAV*sin(deg2rad(THETA_INIT))/(MASS - Z_W_DOT);

    (M_U + M_W_DOT*Z_U/(MASS - Z_W_DOT))/I_YY (M_W + M_W_DOT*Z_W/(MASS - Z_W_DOT))/I_YY ...
    (M_Q + M_W_DOT*(Z_Q + MASS*equilSpeed)/(MASS - Z_W_DOT))/I_YY ...
    -M_W_DOT*MASS*GRAV*sin(deg2rad(THETA_INIT))/(MASS - Z_W_DOT)/I_YY;
    
    0 0 1 0];

longControl = [X_DE/MASS X_DT/MASS;
    Z_DE/(MASS - Z_W_DOT) Z_DT/(MASS - Z_W_DOT);
    (M_DE + M_W_DOT*Z_DE/(MASS - Z_W_DOT))/I_YY (M_DT + M_W_DOT*Z_DE/(MASS - Z_W_DOT)/I_YY);
    0 0];

% Lateral dyamics
latStab = [Y_V/MASS Y_P/MASS Y_R/MASS-equilSpeed GRAV*cos(deg2rad(THETA_INIT));
    (I_ZZ*L_V + I_XZ*N_V)/I_C (I_ZZ*L_P + I_XZ*N_P)/I_C (I_ZZ*L_R + I_XZ*N_R)/I_C 0;
    (I_XZ*L_V + I_XX*N_V)/I_C (I_XZ*L_P + I_XX*N_P)/I_C (I_XZ*L_R + I_XX*N_R)/I_C 0;
    0 1 tan(deg2rad(THETA_INIT)) 0];    

latControl = [Y_DA/MASS Y_DR/MASS;
    (I_ZZ*L_DA + I_XZ*N_DA)/I_C (I_ZZ*L_DR + I_XZ*N_DR)/I_C;
    (I_XZ*L_DA + I_XX*N_DA)/I_C (I_XZ*L_DR + I_XX*N_DR)/I_C;
    0 0];

% Define transformation matrices to use problem vector selection
transMatLong = [equilSpeed 0 0 0; 0 equilSpeed 0 0; 0 0 1 0; 0 0 0 1];
transMatLat = [equilSpeed 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];

% Compute required matrices
longStabMat = transMatLong\longStab*transMatLong;
longContMat = transMatLong\longControl;

latStabMat = transMatLat\latStab*transMatLat;
latContMat = transMatLat\latControl;

%% DYNAMIC MODES & ANALYSIS
% Obtain eigenvectors and eigenvalues
[longEigenVect, longEigenVals] = eig(longStabMat);
[latEigenVect, latEigenVals] = eig(latStabMat);

% Obtain natural frequencies and damping ratios
[natFreqLong, dampRatioLong] = damp(longStabMat);
[natFreqLat, dampRatioLat] = damp(latStabMat);      % Non-oscillatory modes have a value but should not be considered

% Obtain time periods
timePerLong = 2*pi./abs(imag(diag(longEigenVals)));
timePerLat = 2*pi./abs(imag(diag(latEigenVals)));

% Obtain time-to-half
timeToHalfLong = log(2)./abs(real(diag(longEigenVals)));
timeToHalfLat = log(2)./abs(real(diag(latEigenVals)));

% Obtain time constant (only for non-oscillatory modes)
timeConstLat = 1./abs(real(diag(latEigenVals)));

%% MODE DETERMINATION
% Longitudinal modes
longEig = diag(longEigenVals);
idx = find(imag(longEig) > 0);  

% Initialize state labels
stateLabels = {'\Delta u / V', '\alpha', 'q', '\theta'};

for k = 1:length(idx)
    eigIdx = idx(k);
    modeVec = longEigenVect(:, eigIdx);

    % Normalize eigenvector so that pitch angle (4th state) = 1
    modeVecNorm = modeVec / modeVec(4);

    fprintf('\n')
    fprintf('Mode %d: Eigenvalue = %.4f %+.4fi\n', k, real(longEig(eigIdx)), imag(longEig(eigIdx)));
    fprintf('%-10s | %-10s | %-10s\n', 'State', 'Magnitude', 'Phase (deg)');
    fprintf('----------------------------------------\n');

    for i = 1:4
        mag = abs(modeVecNorm(i));
        phs = rad2deg(angle(modeVecNorm(i)));  % Phase in degrees
        fprintf('%-10s | %-10.4f | %-10.2f\n', stateLabels{i}, mag, phs);
    end

    fprintf('\n');
end

%% Lateral modes

eigLat = diag(latEigenVals);

% Identify modes
spiral_idx = find(imag(eigLat) == 0 & abs(real(eigLat)) < 0.05);       % Small negative real
roll_idx   = find(imag(eigLat) == 0 & abs(real(eigLat)) > 0.1);        % Fast real pole
dutch_idxs = find(imag(eigLat) > 0);                                    % Complex pair

% Define state labels
stateLabels = {'\beta', 'p', 'r', '\phi'};

% Spiral and Roll eigenvectors
fprintf('Spiral Mode: Eigenvalue = %.6f\n', eigLat(spiral_idx));
fprintf('Roll Mode:   Eigenvalue = %.4f\n\n', eigLat(roll_idx));

fprintf('%-10s | %-10s | %-10s\n', 'State', 'Spiral', 'Roll');
fprintf('------------------------------------------\n');
for i = 1:4
    spiral_val = latEigenVect(i, spiral_idx) / latEigenVect(4, spiral_idx);  % Normalize w.r.t phi
    roll_val   = latEigenVect(i, roll_idx)   / latEigenVect(4, roll_idx);
    fprintf('%-10s | %-10.4g | %-10.4g\n', stateLabels{i}, spiral_val, roll_val);
end

% Dutch roll eigenvector (first of the conjugate pair)
fprintf('\nDutch Roll Mode: Eigenvalue = %.5f %+.5fi\n\n', real(eigLat(dutch_idxs(1))), imag(eigLat(dutch_idxs(1))));
fprintf('%-10s | %-10s | %-10s\n', 'State', 'Magnitude', 'Phase (deg)');
fprintf('------------------------------------------\n');
VR = latEigenVect(:, dutch_idxs(1));
VR = VR / VR(4);  % Normalize w.r.t. phi

for i = 1:4
    mag = abs(VR(i));
    phase = angle(VR(i)) * (180/pi);
    fprintf('%-10s | %-10.4g | %-10.1f\n', stateLabels{i}, mag, phase);
end

%% EFFECT OF VARYING PITCH STIFFNESS
k = 0:10;
varPitchStifUp = M_W*(1 + 0.1*k);
varPitchStifDown = M_W*(1 - 0.1*k);
varPitchStif = sort([varPitchStifDown, varPitchStifUp]);

plot(real(diag(longEigenVals)), imag(diag(longEigenVals)), 'kx', 'MarkerSize', 10)
grid on
xlabel('Real part'); ylabel('Imaginary part')
title('Variation of longitudinal eigenvalues with pitch stiffness, $\bar{\mathcal{M}}_\alpha$')
hold on

nValues = length(varPitchStif);
nominalIndex = find(varPitchStif == M_W);

% Plot dummy markers for legend
h_kx = plot(nan, nan, 'kx');
h_rx = plot(nan, nan, 'rx');
h_gx = plot(nan, nan, 'gx');

for i = 1:nValues

    if i == nominalIndex(1) || i == nominalIndex(2)
        marker = 'kx';
    elseif varPitchStif(i) > M_W
        marker = 'rx';
    else
        marker = 'gx';
    end

    longStabMat(3, 2) = (varPitchStif(i) + M_W_DOT*Z_W/(MASS - Z_W_DOT))/I_YY*equilSpeed;
    [~, eigenVal] = eig(longStabMat);

    plot(real(diag(eigenVal)), imag(diag(eigenVal)), marker)
end

% Create the legend using the dummy handles
legend([h_kx, h_rx, h_gx], ...
    'Current value of $\bar{\mathcal{M}}_\alpha$', ...
    'Decrease of $|\bar{\mathcal{M}}_\alpha|$', ...
    'Increase of $|\bar{\mathcal{M}}_\alpha|$', ...
    'Location', 'best')

hold off

%% EFFECT OF VARYING ROLL STIFFNESS
varRollStifUp = L_V*(1 + 0.1*k);
varRollStifDown = L_V*(1 - 0.1*k);
varRollStif = sort([varRollStifDown, varRollStifUp]);

plot(real(diag(latEigenVals)), imag(diag(latEigenVals)), 'kx')
grid on
xlabel('Real part'); ylabel('Imaginary part')
title('Variation of lateral eigenvalues with roll stiffness, $\bar{\mathcal{L}}_\beta$')
hold on

nValues = length(varRollStif);
nominalIndex = find(varRollStif == L_V);

% Plot dummy markers for legend
h_kx = plot(nan, nan, 'kx');
h_rx = plot(nan, nan, 'rx');
h_gx = plot(nan, nan, 'gx');

for i = 1:nValues

    if i == nominalIndex(1) || i == nominalIndex(2)
        marker = 'kx';
    elseif varRollStif(i) > L_V
        marker = 'rx';
    else
        marker = 'gx';
    end

    latStabMat(2, 1) = (I_ZZ*varRollStif(i) + I_XZ*N_V)/I_C*equilSpeed;
    latStabMat(3, 1) = (I_XZ*varRollStif(i) + I_XX*N_V)/I_C*equilSpeed;

    [~, eigenVal] = eig(latStabMat);

    plot(real(diag(eigenVal)), imag(diag(eigenVal)), marker)
end

legend([h_kx, h_rx, h_gx], ...
    'Current value of $\bar{\mathcal{L}}_\beta$', ...
    'Decrease of $|\bar{\mathcal{L}}_\beta|$', ...
    'Increase of $|\bar{\mathcal{L}}_\beta|$', ...
    'Location', 'best')

hold off

%% EFFECT OF VARYING YAW STIFFNESS
varYawStifUp = N_V*(1 + 0.1*k);
varYawStifDown = N_V*(1 - 0.1*k);
varYawStif = sort([varYawStifDown, varYawStifUp])

plot(real(diag(latEigenVals(3,3))), imag(diag(latEigenVals(3,3))), 'kx')
grid on
xlabel('Real part'); ylabel('Imaginary part')
title('Variation of lateral eigenvalues with roll stiffness, $\bar{\mathcal{N}}_\beta$')
hold on

nValues = length(varYawStif);
nominalIndex = find(varYawStif == N_V);

% Plot dummy markers for legend
h_kx = plot(nan, nan, 'kx');
h_rx = plot(nan, nan, 'rx');
h_gx = plot(nan, nan, 'gx');

for i = 1:nValues

    if i == nominalIndex(1) || i == nominalIndex(2)
        marker = 'kx';
    elseif varYawStif(i) < N_V
        marker = 'rx';
    else
        marker = 'gx';
    end

    latStabMat(2, 1) = (I_ZZ*L_V + I_XZ*varYawStif(i))/I_C*equilSpeed;
    latStabMat(3, 1) = (I_XZ*L_V + I_XX*varYawStif(i))/I_C*equilSpeed;

    [~, eigenVal] = eig(latStabMat);

    plot(real(diag(eigenVal)), imag(diag(eigenVal)), marker)
end

legend([h_kx, h_rx, h_gx], ...
    'Current value of $\bar{\mathcal{N}}_\beta$', ...
    'Decrease of $|\bar{\mathcal{N}}_\beta|$', ...
    'Increase of $|\bar{\mathcal{N}}_\beta|$', ...
    'Location', 'best')

hold off

% EOF
