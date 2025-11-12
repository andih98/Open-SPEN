%% SPEN-SE-EPI Sequence Generation with Pulseq
% This script generates a ss SPEN Spin-Echo EPI sequence with the 
% chirped-RF pulse as excitation pulse, using the Pulseq toolbox for MATLAB.


%% --- Initialization ---
clear all;
close all;
clc;

% Add the path of the Pulseq toolbox here (https://github.com/pulseq/pulseq.git)

% This section defines all user-configurable parameters.

% Imaging Parameters
fov            = [220e-3, 220e-3, 220e-3]; % Field of View [x, y, z] in m
Nx             = 100;                      % Image resolution in x-direction (readout)
Ny             = 100;                      % Image resolution in y-direction (phase encoding)
TE             = 70e-3;                    % Echo time in s
sliceThickness = 5e-3;                     % Slice thickness in m

% SPEN-specific Parameters
R              = 2;     % Phase undersampling due to swept RF pulse. Increases robustness towards off-resonance.
rampSampling   = true;  % `true` for ramp sampling, `false` for regular sampling
ro_os          = 1;     % Oversampling factor in readout direction
spenNavigator  = false; % `true` to enable the SPEN navigator


%% --- System Definition ---
% Definition of the MR system limits.

% Standard system limits
sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
              'MaxSlew', 150, 'SlewUnit', 'T/m/s', ...
              'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6, ...
              'adcDeadTime', 10e-6, 'B0', 3);

% System limits for crusher gradients (allows higher slew rate)
sysCrusher = sys;
sysCrusher.maxSlew = 70 * sysCrusher.gamma; % Slew rate in Hz/m/s

% Initialize the sequence object
seq = mr.Sequence(sys);

deltak    = 1 ./ fov(2); % Pulseq uses k-space units of 1/m
kWidth    = Nx * deltak;

% 1. Excitation Pulse (Chirped Pulse for SPEN)
rfdur   = 20e-3;          % Duration of the RF pulse
sweepBw = R * Ny / rfdur; % Sweep bandwidth for SPEN encoding

gexc = mr.makeTrapezoid('y', sys, 'Amplitude', sweepBw / fov(1), ...
                        'FlatTime', rfdur, 'RiseTime', sys.rfDeadTime);
                        
% First attempt to create the chirped pulse.
% The function `makeChirpedRfPulse2` is assumed to exist.
rf = makeChirpedRfPulse('duration', rfdur, 'delay', sys.rfDeadTime, ...
                         'bandwidth', sweepBw, 'ang', 90, 'n_fac', 40, ...
                         'system', sys, 'use', 'excitation');
                         
% Calculate the bandwidth of the first pulse to use the actual bandwidth
% for a more accurate second calculation. This corrects for deviations
% caused by the discrete pulse shape.
[bw, ~, ~, ~] = mr.calcRfBandwidth(rf);

% Re-create the chirped pulse with the corrected bandwidth
rf = makeChirpedRfPulse('duration', rfdur, 'delay', sys.rfDeadTime, ...
                         'bandwidth', sweepBw * sweepBw / bw, 'ang', 90, 'n_fac', 40, ...
                         'system', sys, 'use', 'excitation');
                         
% Save the final bandwidth for metadata
[bw, ~, ~, ~] = mr.calcRfBandwidth(rf);

% 2. Refocusing Pulse (180° Sinc Pulse) and Crusher
dur_ref = 3e-3;
[rfref, gz] = mr.makeSincPulse(pi, sys, 'Duration', dur_ref, ...
                              'SliceThickness', sliceThickness, ...
                              'apodization', 0.2, 'timeBwProduct', 4, ...
                              'PhaseOffset', 0, 'use', 'refocusing');
% Crusher gradient; the factor of 3 is relatively high, but kept from the original code.
gzCrush = mr.makeTrapezoid('z', sysCrusher, 'Area', -gz.area * 3);

% 3. EPI Readout with Ramp Sampling
% Duration of the blip, rounded up to a multiple of twice the gradient raster time
blip_dur = ceil(2 * sqrt(gexc.flatArea / Ny / sys.maxSlew) / (2 * sys.gradRasterTime)) ...
           * (2 * sys.gradRasterTime);

% Blip gradient in y-direction
gacq = mr.makeTrapezoid('y', sys, 'Area', gexc.flatArea / Ny, 'Duration', blip_dur);

% Readout gradient (gRO) in x-direction with correction for ramp sampling
roDur = 5.2e-4; % Duration of the flat top of the readout gradient

% Additional area for the ramps to center the k-space trajectory
extra_area = (blip_dur / 2)^2 * sys.maxSlew; 
gRO_target_area = kWidth; % Desired k-space width

gRO = mr.makeTrapezoid('x', sys, 'Area', gRO_target_area + extra_area, ...
                       'duration', roDur + blip_dur);

% Correction of the gradient amplitude to cover exactly `kWidth`, since the
% ADC is not active during the entire gradient duration.
actual_area = gRO.area - gRO.amplitude / gRO.riseTime * (blip_dur/2)^2 / 2 ...
                       - gRO.amplitude / gRO.fallTime * (blip_dur/2)^2 / 2;
gRO.amplitude = gRO.amplitude / actual_area * gRO_target_area;

% Update gradient areas after amplitude adjustment
gRO.area = gRO.amplitude * (gRO.flatTime + gRO.riseTime / 2 + gRO.fallTime / 2);
gRO.flatArea = gRO.amplitude * gRO.flatTime;

% 4. ADC Definition for Ramp Sampling
adcDwellNyquist = deltak / gRO.amplitude / ro_os;
adcDwell = floor(adcDwellNyquist * 1e7) * 1e-7; % Round down dwell time to 100 ns
adcSamples = floor(roDur / adcDwell / 4) * 4;  % Number of samples, divisible by 4 (Siemens requirement)

adc = mr.makeAdc(adcSamples, 'Dwell', adcDwell, 'Delay', blip_dur / 2);

% Align ADC in time to the center of the gradient plateau
time_to_center = adc.dwell * ((adcSamples - 1) / 2 + 0.5); % Time to the center of the middle sample
adc.delay = round((gRO.riseTime + gRO.flatTime / 2 - time_to_center) * 1e6) * 1e-6; % Round delay to 1 µs

% The actual matrix size in the x-direction is determined by the number of samples
Nx = adcSamples;

% 5. Prepare Blip Gradients for EPI Readout
gacq_parts = mr.splitGradientAt(gacq, blip_dur / 2, sys);
[gacq_blipup, gacq_blipdown, ~] = mr.align('right', gacq_parts(1), 'left', gacq_parts(2), gRO);
gacq_blipdownup = mr.addGradients({gacq_blipdown, gacq_blipup}, sys);
durAcq = Ny * mr.calcDuration(gRO); % Total duration of the EPI readout

% 6. Other Gradients (SPEN Rephaser, Navigators, Crushers)
gre = mr.makeTrapezoid('y', sys, 'Area', -gexc.area / 2); % Rephaser after excitation
gre2 = mr.makeTrapezoid('y', sys, 'Area', -gexc.flatArea / 2); % Rephaser after refocusing

gxCrushArea = 0; % Placeholder for additional dephasing
gNavPre  = mr.makeTrapezoid('x', sysCrusher, 'Area', -gRO.area / 2);
gxCrush1 = mr.makeTrapezoid('x', sysCrusher, 'Area', gRO.area / 4 + gxCrushArea);
gxCrush2 = mr.makeTrapezoid('x', sysCrusher, 'Area', -gRO.area / 4 + gxCrushArea);

% SPEN Navigator (optional)
quadPhaseNav = 100;
quadPhaseNavDur = quadPhaseNav * (round(R * roDur / quadPhaseNav / sys.adcRasterTime) * sys.adcRasterTime);
spenRO = mr.makeTrapezoid('y', sys, 'flatArea', gexc.flatArea, 'FlatTime', quadPhaseNavDur);
spenADC = mr.makeAdc(quadPhaseNav, sys, 'Duration', spenRO.flatTime, 'Delay', spenRO.riseTime);
spenNav = mr.makeTrapezoid('y', sys, 'Area', -spenRO.area / 2 - deltak * R / 2);


%% --- Timing Calculation (Delays) ---
% Calculation of the necessary delays to meet the exact echo time (TE).

% Delay 1: Between excitation and refocusing pulse
duration_pre_refocus = mr.calcDuration(gexc) / 2 ...  % Half of the excitation block
                     + mr.calcDuration(gre) ...
                     + 3 * mr.calcDuration(gRO) ...     % k0-Navigator
                     + 2 * mr.calcDuration(gNavPre) ...
                     + mr.calcDuration(gxCrush1, gzCrush) ...
                     + mr.calcDuration(gz) / 2;         % Half of the refocusing block
if spenNavigator
    duration_pre_refocus = duration_pre_refocus ...
                         + 5 * mr.calcDuration(spenRO) ... % SPEN-Navigator
                         + 2 * mr.calcDuration(spenNav);
end
delay1 = TE / 2 - duration_pre_refocus;
delay1 = round(max(delay1, 0) / sys.gradRasterTime) * sys.gradRasterTime;

% Delay 2: Between refocusing pulse and EPI readout
duration_post_refocus = mr.calcDuration(gz) / 2 ...       % Second half of the refocusing block
                      + mr.calcDuration(gxCrush2, gzCrush, gre2) ...
                      + (1 + Ny / 2) * mr.calcDuration(gRO); % Assumption: Echo is in the center of the readout
delay2 = TE / 2 - duration_post_refocus;
delay2 = round(max(delay2, 0) / sys.gradRasterTime) * sys.gradRasterTime;


%% --- Sequence Assembly ---
% The individual events are added to the sequence object in the correct order.

% 1. SPEN Excitation and Rephasing
seq.addBlock(rf, gexc);
seq.addBlock(gre);

% 2. k0-Navigator for EPI phase corrections (3 echoes)
seq.addBlock(gNavPre);
seq.addBlock(gRO, adc);
gRO = mr.scaleGrad(gRO, -1);
seq.addBlock(gRO, adc);
gRO = mr.scaleGrad(gRO, -1);
seq.addBlock(gRO, adc);
seq.addBlock(gNavPre);

% 3. SPEN Navigator (optional)
if spenNavigator
    seq.addBlock(spenNav);
    for k = 1:5 % Five readouts with alternating polarity
        seq.addBlock(spenRO, spenADC);
        spenRO = mr.scaleGrad(spenRO, -1);
    end
    seq.addBlock(spenNav);
end

% 4. First Half of the Spin Echo
seq.addBlock(mr.makeDelay(delay1));
seq.addBlock(gxCrush1, gzCrush);

% 5. Slice-selective Refocusing Pulse
seq.addBlock(rfref, gz);

% 6. Second Half of the Spin Echo
seq.addBlock(gxCrush2, gzCrush, gre2);
seq.addBlock(mr.makeDelay(delay2));

% 7. Continuous EPI Readout
for i = 1:Ny
    if i == 1
        seq.addBlock(gRO, gacq_blipup, adc);   % First line
    elseif i == Ny
        seq.addBlock(gRO, gacq_blipdown, adc); % Last line
    else
        seq.addBlock(gRO, gacq_blipdownup, adc); % Middle lines
    end
    gRO = mr.scaleGrad(gRO, -1); % Invert polarity for the next line
end


%% --- Finalization and Export ---

% Check the sequence for compliance with system limits
fprintf('Checking the sequence...\n');
[ok, error_report] = seq.checkTiming;
if ok
    fprintf('Timing check successful.\n');
    fprintf('Sequence duration: %.3f ms, TE: %.3f ms\n', seq.duration() * 1e3, TE * 1e3);
else
    fprintf('Timing check failed!\n');
    disp(error_report);
end


% Definitions for the .seq file (metadata)
seq_name = 'SPEN_SE_EPI';
seq.setDefinition('Name', seq_name);
seq.setDefinition('FOV', fov);
seq.setDefinition('Nx', Nx);
seq.setDefinition('Ny', Ny);
seq.setDefinition('TE', TE);
seq.setDefinition('sliceThickness', sliceThickness);
seq.setDefinition('ReceiverGainHigh', 1);
seq.setDefinition('DurAcq', durAcq);
seq.setDefinition('R', R);
seq.setDefinition('sweepBw', bw);
seq.setDefinition('rfref_dur', rfdur); 
seq.setDefinition('Gspen', gexc.amplitude);
seq.setDefinition('SeqDur', seq.duration());
seq.setDefinition('RampSampling', rampSampling);
seq.setDefinition('gROriseFlatAmp', [gRO.riseTime, gRO.flatTime, gRO.amplitude]);
seq.setDefinition('adcDelayDwellNsamples', [adc.delay, adc.dwell, adcSamples]);
seq.setDefinition('ro_os', ro_os);
seq.setDefinition('quadPhaseNavDur', quadPhaseNavDur);
seq.setDefinition('quadPhaseNav', quadPhaseNav);
seq.setDefinition('spenNavigator', spenNavigator);

% Optional: Write sequence file and install on scanner
% seq.write([seq_name '.seq']);
% seq.install('');

% Optional: PNS Calculation (Peripheral Nerve Stimulation)
% [pns_ok, pns_n, pns_c, tpns] = seq.calcPNS('.asc');

%% Check
rep = check(seq, TE, [], []);
fprintf([rep{:}]);

%%
% Optional: Plot the RF pulse profile
[rf_bw, ~, rf_spectrum, rf_w] = mr.calcRfBandwidth(rf);
figure; plot(rf_w, abs(rf_spectrum));
title('Excitation pulse profile (small-angle approximation)');
xlabel('Frequency [Hz]');
xlim(2 * [-rf_bw, rf_bw]);
grid on

%%
simulateSPENfo.spins_per_voxel = 10;

% simulateSPENfo=[];

deriveSPENforwardOperator('exc',seq,sys,simulateSPENfo)