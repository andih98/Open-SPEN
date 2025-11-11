function deriveSPENforwardOperator(seq, sys, simulateSPENfo)
% DERIVESPENFORWARDOPERATOR Calculates the forward operator for a SPEN
% (SPatiotemporal ENcoding) sequence, and optionally simulates the SPEN Forward Operator using a Bloch
% equation model, adapting code from a Stanford University Lecture: "Bloch
% Equation Matrix Simulations" (https://www.stanford.edu/)
%
% Args:
%   seq (mr.Sequence): The sequence object containing the SPEN parameters.
%   sys (mr.SystemSpec): The system specification object.
%   simulateSPENfo (struct/empty): A struct containing simulation
%     parameters (e.g., 'spins_per_voxel') or empty if no simulation is
%     desired.

%% --- Parameter Extraction and Sequence Definition ---

% Extract relevant parameters from the sequence definition
param.fov = seq.getDefinition('FOV');       % The Field of View in the phase encoded dimension [m]
param.Ny = seq.getDefinition('Ny');         % Number of steps in the phase-encoding dimension
param.R = seq.getDefinition('R');           % Parallel imaging acceleration factor
param.rfdur = seq.getDefinition('rfref_dur'); % Duration of the RF pulse [s]
param.sweepBw = seq.getDefinition('sweepBw'); % Bandwidth of the chirped RF pulse [Hz]
param.g = seq.getDefinition('Gspen');       % SPEN gradient amplitude
param.spenNavigator = seq.getDefinition('spenNavigator'); % SPEN navigator flag (unused in this core logic)

% Initialize a new sequence object based on system specs (important for building the sequence)
seq = mr.Sequence(sys);
roDur = param.rfdur; % Readout duration is set equal to the RF pulse duration

%% --- Gradient and RF Pulse Generation ---

% Excitation Gradient (Phase Encoding/SPEN)
% Area = R/FOV * Ny (normalized k-space units, scaled for the simulation context)
gexc = mr.makeTrapezoid('y', sys, 'FlatArea', param.R/param.fov(2)*param.Ny, ...
                        'FlatTime', param.rfdur, 'Delay', sys.rfDeadTime);

% Chirped RF Excitation Pulse (90-degree)
rf = makeChirpedRfPulse('duration', param.rfdur, ...
                         'Delay', sys.rfDeadTime + gexc.riseTime, ...
                         'bandwidth', param.sweepBw, ...
                         'ang', 90, 'n_fac', 40, 'system', sys, 'use', 'excitation');

% Acquisition/Readout Gradient
% The actual gradient applied during the ADC, its area is crucial for the chirp encoding
gacq = mr.makeTrapezoid('y', sys, 'Area', -gexc.flatArea, 'Duration', roDur);

% Analog-to-Digital Converter (ADC) for readout
adc = mr.makeAdc(param.Ny, 'duration', gacq.flatTime, 'delay', gacq.riseTime);

%% --- Sequence Construction ---

% Add excitation block (chirped-RF and Gexc)
seq.addBlock(rf, gexc);
% Add readout block (Gacq and ADC)
seq.addBlock(gacq, adc);

%% --- K-space Calculation and Forward Operator Derivation ---

% Calculate k-space trajectory
[kspace.ktraj_adc, kspace.t_adc, kspace.ktraj, kspace.t_ktraj, kspace.t_excitation, kspace.t_refocusing] = seq.calculateKspacePP();

% Calculate the analytical forward operator kernel
% This function is assumed to compute the kernel based on the k-space trajectory and sequence parameters
kernel = forwardOperatorCalc(kspace, param);

%% --- Calculated Forward Operator Visualization ---

% Calculate the phase of the kernel (used for visualization)
phase = angle(mfft(kernel, 2));

% Plot 1: Quadratic Phase Acquisition Kernel
figure;
imagesc([-param.fov(1)/2, +param.fov(1)/2], linspace(0, 1, size(kernel, 1)), (phase));
colormap("turbo");
axis('tight', 'square');
title('Quadratic Phase Acquisition Kernel');
ylabel('Normalized Acquisition Time');
xlabel('FOV [m]');
c = colorbar;
c.Label.String = 'Phase [rad]';
c.Label.FontSize = 12;

% Plot 2: Calculated Forward Operator (Magnitude)
N = size(kernel, 2); % Number of y-coordinates in the kernel
Kvec = (-N/2:N/2-1) * param.sweepBw / size(kernel, 2); % Frequency vector after FFT shift [Hz]
figure;
imagesc(Kvec, [1, size(kernel, 1)], abs(fliplr(kernel)) / max(abs(fliplr(kernel)), [], "all"));
colormap("parula");
axis('tight', 'square');
title('Calculated Forward-Operator');
ylabel('Acquisitions');
xlabel('Frequency [Hz]');
c = colorbar;
c.Label.String = 'Signal Contribution';
c.Label.FontSize = 12;

% Plot 3: Point Spread Function (PSF) in Image Domain (sum over time points)
psf = real((sum(kernel, 1) / max(real(sum(kernel, 1)))));
figure;
plot(1:length(psf), psf);
title('Image Domain Point Spread Function (Sum of Kernel)');
xlabel('Spatial Position Index');
ylabel('Normalized Amplitude');
grid on;

%% --- Optional Bloch Simulation ---
if ~isempty(simulateSPENfo)
    dt = sys.rfRasterTime; % Time step for simulation (RF raster time)

    % Get gradient waveforms and times from the sequence
    gw_data = seq.waveforms_and_times();
    % Extract the y-gradient waveform
    gy = mr.pts2waveform(gw_data{2}(1, :), gw_data{2}(2, :), sys.rfRasterTime);

    gyWaveform = gy;

    % Align RF waveform with the gradient waveform time
    rfWaveform = [zeros(1, round(gexc.riseTime / dt)) rf.signal];
    rfWaveform = [rfWaveform zeros(1, length(gy) - length(rfWaveform))];

    % Create ADC sampling vector
    samples = zeros(1, seq.blockDurations(2) / dt);
    samples(1:end) = 1; % Assume continuous sampling during the readout block
    adc_samples = [zeros(1, round(seq.blockDurations(1) / dt - 100)) samples];

    % --- Simulation Setup and Parameters ---
    fov1 = param.fov(1) * 1000; % FOV in mm (y-dimension)
    t = 1:length(rfWaveform);   % Time index vector
    rf1 = rfWaveform;           % RF signal
    T1 = 9000000;               % T1 relaxation time [ms] (very long for minimal relaxation)
    T2 = 900000;                % T2 relaxation time [ms] (very long for minimal relaxation)
    dT = dt;                    % Time step
    gamma = 4258;               % Gyromagnetic ratio [Hz/G]

    % Spin Positions
    % A set of positions covering the FOV
    num_spins = param.Ny * simulateSPENfo.spins_per_voxel;
    pos = linspace(-fov1/2, fov1/2 - fov1/num_spins, num_spins);
    pos = pos(:).'; % Row vector [mm]

    % Convert gradient from [T/m] to [G/cm]
    grad = ((gyWaveform / sys.gamma) / 1e-4) / 100; % [G/cm]

    % Convert RF from [rad/s] (implicit in mr library) to rotation angle [rad]
    rfrot = 2 * pi * rf1 * dT; % RF rotation angle over dT [rad]

    % Proton Density (Phantom)
    pd = ones(size(pos));
    % Place a single spin at the center for impulse response testing
    % pd(round(length(pos)/2)) = 1;

    % ADC indices
    adc_idx = find(adc_samples == 1);
    nADC = numel(adc_idx);   % Number of ADC samples
    numPos = numel(pos);     % Number of simulated spin positions
    nRF = numel(rf1);        % Number of time steps

    % Preallocation
    msig2 = complex(zeros(numPos, nADC)); % Matrix to store magnetizations M_xy (time x position)
    msig = complex(zeros(numPos, 1));     % Final transverse magnetization
    m = zeros(3, numPos);                 % Final M_x, M_y, M_z

    % --- Precompute Rotation Matrices and Relaxation Operators ---
    Rrf = cell(nRF, 1);
    A_half = cell(nRF, 1);
    B_half = cell(nRF, 1);
    Rg_half = cell(nRF, 1);

    for k = 1:nRF
        % RF Rotation matrix (around the effective B1 field direction)
        Rrf{k} = throt_fast(abs(rfrot(k)), angle(rfrot(k))); % Assumes 'throt_fast' is available

        % Relaxation operator (T1/T2) for dT/2 (T in ms!)
        [A_half{k}, B_half{k}] = freeprecess(1000 * dT / 2, T1, T2, 0); % Assumes 'freeprecess' is available

        % Gradient Rotation matrix for all positions over dT/2
        % pos/10 converts mm to cm
        angles = 2 * pi * gamma * (pos / 10) * grad(k) * dT / 2; % [1 × numPos] Phase angle [rad]
        Rg_half{k} = reshape(...
            [cos(angles); -sin(angles); zeros(1, numPos);
             sin(angles);  cos(angles); zeros(1, numPos);
             zeros(1, numPos); zeros(1, numPos); ones(1, numPos)], ...
            3, 3, numPos); % 3x3 rotation matrix for each spin
    end

    % --- Main Bloch Simulation Loop ---
    for x = 1:numPos % Loop over all spin positions

        if pd(x) == 0
            continue; % Skip if proton density is zero
        end

        % Progress report
        fprintf('Simulating position %d of %d\n', x, numPos);

        M = [0; 0; 1]; % Initial magnetization: Mz=1
        sig2_loc = complex(zeros(1, nADC)); % Signal for current spin position
        sample = 1;

        for k = 1:nRF % Loop over all time steps
            % 1. Relaxation + Gradient (first dT/2)
            M = A_half{k} * M + B_half{k};
            M = Rg_half{k}(:, :, x) * M;

            % 2. RF Pulse (instantaneous rotation at the center of the step)
            M = Rrf{k} * M;

            % 3. Relaxation + Gradient (second dT/2)
            M = A_half{k} * M + B_half{k};
            M = Rg_half{k}(:, :, x) * M;

            % 4. ADC Acquisition (after the full time step)
            if adc_samples(k)
                % Transverse magnetization is the signal
                sig2_loc(sample) = (M(1) + 1i * M(2)) * pd(x);
                sample = sample + 1;
            end
        end

        % Store results for this spin position
        msig2(x, :) = sig2_loc;
        m(:, x) = M;
        msig(x) = (M(1) + 1i * M(2)); % Final transverse component
    end

    % --- Simulation Results Evaluation & Plotting ---

    % Plot 4: Final Magnetization and Waveforms
    figure;
    subplot(3, 2, 1); plot(pos, abs(msig)); xlabel('Position (mm)'); ylabel('Magnitude'); grid on; title('Magnitude Profile');
    subplot(3, 2, 3); plot(pos, angle(msig)); xlabel('Position (mm)'); ylabel('Phase (rad)'); grid on; title('Phase Profile');
    subplot(3, 2, 5); plot(pos, m(3, :)); xlabel('Position (mm)'); ylabel('$M_z$'); grid on; title('$M_z$ Profile');
    subplot(3, 2, 2); plot(t * dt, rf1); xlabel('Time (s)'); ylabel('RF (a.u.)'); grid on; title('RF vs Time');
    subplot(3, 2, 4); plot(t * dt, grad); xlabel('Time (s)'); ylabel('Grad (G/cm)'); grid on; title('Grad vs Time');
    sgtitle('SPEN Simulation', 'Interpreter', 'latex');
    subplot(3, 2, 6);
    text_string = {
    sprintf('Spins/voxel: %d', simulateSPENfo.spins_per_voxel),
    sprintf('RF-Dauer: %.1e s', param.rfdur),
    sprintf('Sweep BW: %d Hz', param.sweepBw),
    sprintf('R-Faktor: %d', param.R)};
    text(0, 0.5, text_string, 'FontSize', 10, 'Interpreter', 'latex');
    axis off;

    % Simulated Forward Operator (PSF in Frequency Space)
    % Bin the simulated signal (msig2) to match the final image resolution
    % and simulate aliasing/ADC
    psf = binMatrixSum(msig2, simulateSPENfo.spins_per_voxel, 1); % Bin over position (rows)

    % Frequency vector for visualization (after FFT shift)
    N_psf = size(psf, 2);
    Kvec = (-N_psf/2:N_psf/2-1) * param.sweepBw / N_psf; % Frequency vector [Hz]

    % Plot 5: Simulated Forward Operator (Magnitude)
    figure;
    imagesc(Kvec, linspace(0, 1, size(psf, 1)), abs(psf) / max(abs(psf), [], 'all'));
    colormap("parula");
    axis('tight', 'square');
    title('Simulated Forward-Operator');
    ylabel('Normalized Acquisition Time');
    xlabel('Frequency [Hz]');
    c = colorbar;
    c.Label.String = 'Signal Contribution';
    c.Label.FontSize = 12;

    % Scroll Plot 1: PSF in Frequency Space
    scrollPlot(abs(psf) / max(abs(psf), [], 'all'), 1, 'TitlePrefix', 'Simulated PSF Frequency Space: '); % Assumes 'scrollPlot' is available

    % Plot 6: Line plots of the PSF in Frequency Space for different FOV positions
    Kvec_plot = (-N_psf/2:N_psf/2-1) * param.sweepBw / N_psf / 10; % Adjusted K-vector for plotting (e.g., k-space units 1/m)
    figure;
    hold on;
    % Get acquisition time indices corresponding to roughly 1/4, 1/2, 3/4 of the total readout time
    idx1 = round(size(psf, 1) * 0.25);
    idx2 = round(size(psf, 1) * 0.50);
    idx3 = round(size(psf, 1) * 0.75);

    % Plot 1/4 FOV (early time point, wide bandwidth)
    plot(Kvec_plot, abs(psf(idx1, :)) / max(abs(psf(idx1, :))), '-', 'Color', 'k', 'LineWidth', 5);
    % Plot 1/2 FOV (mid time point)
    plot(Kvec_plot, abs(psf(idx2, :)) / max(abs(psf(idx2, :))), '-.', 'Color', [0.4 0.4 0.4], 'LineWidth', 3);
    % Plot 3/4 FOV (late time point, narrow bandwidth)
    plot(Kvec_plot, abs(psf(idx3, :)) / max(abs(psf(idx3, :))), ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 3);

    xlabel('$k_y$ [1/m]', 'Interpreter', 'latex');
    ylabel('Normalized Signal Contribution')
    legend('Position 1/4 Acquisition Time', 'Position 1/2 Acquisition Time', 'Position 3/4 Acquisition Time', 'Location', 'best')
    title('Signal Localization and Convolution in Frequency Domain')
    axis tight
    ylim([min(abs(psf(idx3, :)) / max(abs(psf(idx3, :)))) 1.2])

    % Final PSF in Image Space (by FFT along the frequency/position dimension)
    psf_img = binMatrixSum(msig2, simulateSPENfo.spins_per_voxel, 1); % Bin over position (rows)
    psf_img = binMatrixSum(psf_img, length(psf_img)/param.Ny, 2); % Bin over time (columns) - this reduces the time-domain sampling for the plot
    scrollPlot(abs(mfft(psf_img, 2)) / max(abs(mfft(psf_img, 2)), [], 'all'), 1, 'TitlePrefix', 'Simulated PSF Image Space: '); % FFT along position dimension
end

end
