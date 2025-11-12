%% SPEN-SE-EPI Sequence Generation with Pulseq
% This script generates a Multi-Slice SPEN Spin-Echo EPI sequence 
% with the chirped-RF pulse as refocusing pulse, full-refocusing, 
% and optional diffusion-weighting, using the Pulseq toolbox for MATLAB.

clear all
close all
clc
%% 1. SYSTEM LIMITS AND SETUP
% Define MR system limits.
sys = mr.opts('MaxGrad', 60, 'GradUnit', 'mT/m', ...
    'MaxSlew', 180, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6, ...
    'adcDeadTime', 10e-6, 'B0', 3);

% Create separate system objects for crushers and diffusion gradients
% This allows for different slew rate limits if necessary.
sysCrusher = sys;
sysCrusher.maxSlew = 70 * sysCrusher.gamma; % Custom slew for crushers
sysDiff = sys;
sysDiff.maxSlew = 90 * sysDiff.gamma;      % Custom slew for diffusion

% Initialize the sequence object
seq = mr.Sequence(sys);

%% 2. SEQUENCE PARAMETERS 
fov = [220e-3, 220e-3, 220e-3]; % Field of View [x, y, z] (m)
Nx = 100;                      % Matrix size in x (readout)
Ny = 100;                      % Matrix size in y (SPEN/EPI)
TR = 1000e-3;                  % Repetition Time (s)
sliceThickness = 5e-3;         % Slice thickness (m)
numSlices = 10;                % Number of slices
avg = 1;                       % Number of averages

% Calculate slice spacing to cover FOV(z)
spacing = round((fov(3) - numSlices * sliceThickness) / (numSlices - 1), 6);

% k-space parameters
deltak = 1 / fov(2);           % k-space step (m^-1)
kWidth = Nx * deltak(1);       % Total k-space width in x

%% 3. SPEN/EPI PARAMETERS
R = 2;                         % SPEN undersampling factor due to swept RF pulse
rampSampling = 1;              % Flag for using trapezoidal readout gradients (true)
ro_os = 2;                     % Readout oversampling factor
randomSampling = 0;            % Flag for random blip amplitudes
rfref_dur = 31e-3;             % Duration of the SPEN refocusing RF pulse (s)

% Diffusion gradient amplitudes (in T/m). Set to 0 for b=0 acquisition.
x_gDiff = 60e-3;
y_gDiff = 0e-3;
z_gDiff = 0e-3;
distributeDiffArea = 0.9;      % Factor for diffusion gradient distribution (before ref. - after ref.)
addDiffDur = 0e-3;             % Additional duration for diffusion gradient after ref.

%% 4. RF PULSE DEFINITIONS
% --- Excitation pulse (90-degree)
[rf, gz, gzRe] = mr.makeSincPulse(pi/2, sysCrusher, 'Duration', 2e-3, ...
    'SliceThickness', sliceThickness, 'apodization', 0.42, ...
    'timeBwProduct', 4, 'PhaseOffset', 0, 'use', 'excitation');

%% 5. Hybrid/SPEN -SPACE SAMPLING STRATEGY
% Total k-space area to be covered in the y-direction (SPEN encoding)
area = 2 * Ny * deltak * R;

% Generate the vector of blip areas ('rs')
if randomSampling
    % Generate random variations around the mean blip area
    rs = (rand(1, Ny) * 2 - 1);
    rs = (rs - mean(rs)) * area / (Ny) / 3;
    rs = [ones(1, Ny) * area / (Ny) + rs, 0];
else
    % Uniform blip areas
    rs = [0, ones(1, Ny) * area / Ny];
end

%% 6. SPEN REFOCUSING PULSE AND GRADIENT
% This is the core of SPEN: a chirped (swept) RF pulse played simultaneously
% with a strong gradient ('gref') to encode spatial information.

% Calculate sweep bandwidth (SPEN-cond.)
sweepBw = R * Ny / rfref_dur;

% Define the SPEN gradient ('gref')
gref = mr.makeTrapezoid('y', sys, 'Amplitude', sweepBw / fov(1), ...
    'FlatTime', rfref_dur, 'Delay', sys.rfDeadTime);

% --- SPEN chirped RF refocusing pulse (180-degree)
% This uses a custom function 'makeChirpedRfPulse' (assumed to be on path)
rfref = makeChirpedRfPulse('duration', rfref_dur, ...
    'Delay', sys.rfDeadTime + gref.riseTime, 'bandwidth', sweepBw, ...
    'ang', 180, 'n_fac', 40, 'system', sys, 'use', 'refocusing');

% Bandwidth correction:
% Calculate the actual bandwidth of the pulse
[bw, f0, M_xy_sta, F1] = mr.calcRfBandwidth(rfref);
% Re-create the pulse using the corrected bandwidth for better accuracy
rfref = makeChirpedRfPulse('duration', rfref_dur, ...
    'delay', sys.rfDeadTime + gref.riseTime, 'bandwidth', sweepBw * sweepBw / bw, ...
    'ang', 180, 'n_fac', 40, 'system', sys, 'use', 'refocusing');

% Calculate pulse power properties
[bw, f0, M_xy_sta, F1] = mr.calcRfBandwidth(rfref);
[total_energy, peak_pwr, rf_rms] = mr.calcRfPower(rfref, sys.rfRasterTime);

%% 7. EPI READOUT DEFINITION
% Calculate blip duration, rounding to the gradient raster time
blip_dur = ceil(2 * sqrt(max(rs) / sys.maxSlew) / 10e-6 / 2) * 10e-6 * 2;

% Calculate readout duration (SPEN-cond. / full-refocusing)
roDur = 2 * rfref_dur / Ny - blip_dur

% --- Readout gradient (gRO) with ramp sampling
% Calculate extra area needed to compensate for ramp times
extra_area = blip_dur / 2 * blip_dur / 2 * sys.maxSlew;
gRO = mr.makeTrapezoid('x', sys, 'Area', kWidth + extra_area, 'duration', roDur + blip_dur);

% Rescale gradient to achieve the exact target area 'kWidth'
actual_area = gRO.area - gRO.amplitude / gRO.riseTime * blip_dur / 2 * blip_dur / 2 / 2 ...
                       - gRO.amplitude / gRO.fallTime * blip_dur / 2 * blip_dur / 2 / 2;
gRO.amplitude = gRO.amplitude / actual_area * kWidth;
gRO.area = gRO.amplitude * (gRO.flatTime + gRO.riseTime / 2 + gRO.fallTime / 2);
gRO.flatArea = gRO.amplitude * gRO.flatTime;

% --- ADC definition
% Calculate dwell time and number of samples for ramp sampling
adcDwellNyquist = deltak / gRO.amplitude / ro_os;
adcDwell = floor(adcDwellNyquist * 1e7) * 1e-7; % Round down to 100 ns
adcSamples = floor(roDur / adcDwell / 4) * 4;   % Samples must be divisible by 4

% Create ADC object
adc = mr.makeAdc(adcSamples, 'Dwell', adcDwell, 'Delay', blip_dur / 2);

% Re-align ADC to the center of the gradient flat-top
time_to_center = adc.dwell * ((adcSamples - 1) / 2 + 0.5);
adc.delay = round((gRO.riseTime + gRO.flatTime / 2 - time_to_center) * 1e6) * 1e-6;
Nx = adcSamples; % Update Nx to match the actual number of samples

% Create copies for the first readout line (which has different timing)
gRO1 = gRO;
gRO1.delay = blip_dur / 2;
adc_first = adc;
adc_first.delay = adc_first.delay + blip_dur / 2;
durAcq = Ny * mr.calcDuration(gRO); % Total acquisition duration

%% 8. OTHER PULSE/GRADIENT DEFINITIONS
% --- Fast recovery pulse (moves magnetization to -z)
restoreRF = mr.makeBlockPulse(-pi, sys, 'duration', 600e-6, 'use', 'inversion');

% SPEN rephasing gradient
gre = mr.makeTrapezoid('y', sysCrusher, 'Area', -gref.flatArea, 'delay', mr.calcDuration(gref));
gy = mr.addGradients({gref, gre}, sys); % Combine SPEN gradient and its rephaser

% Readout pre-/re-phasers for k-space navigators
gx1 = mr.makeTrapezoid('x', sysCrusher, 'Area', gRO.area / 4);
gx2 = mr.makeTrapezoid('x', sysCrusher, 'Area', -gRO.area / 4);
gNavPre = mr.makeTrapezoid('x', sysCrusher, 'Area', -gRO.area / 2);

% --- Crusher gradients (placed around the 180-deg pulse)
gzCrush1 = mr.makeTrapezoid('z', sysCrusher, 'area', -gz.area * 4);
gzCrush2 = mr.makeTrapezoid('z', sysCrusher, 'area', -gz.area * 4, 'Delay', mr.calcDuration(gref) - gref.fallTime);
dummyz = mr.makeTrapezoid('z', sysCrusher, 'Area', 0, 'Duration', mr.calcDuration(gref) - sys.rfDeadTime);
gzCrush = mr.addGradients({dummyz, gzCrush2}, sys); % Combined crusher block

gxCrushArea = 0; % Additional x-crushing (if needed)
gxCrush1 = mr.makeTrapezoid('x', sysCrusher, 'Area', gx1.area + gxCrushArea);
gxCrush2 = mr.makeTrapezoid('x', sysCrusher, 'Area', gx2.area + gxCrushArea, 'Delay', mr.calcDuration(gref) - gref.fallTime);
dummyx = mr.makeTrapezoid('x', sysCrusher, 'Area', 0, 'Duration', mr.calcDuration(gref) - sys.rfDeadTime);
gxCrush = mr.addGradients({dummyx, gxCrush2}, sys);
gxCrush1.delay = mr.calcDuration(gzCrush1) - mr.calcDuration(gxCrush1);

% --- Gradient moment nulling (at the end of TR)
null_y = mr.makeTrapezoid('y', sysCrusher, 'Area', -gref.area + 2 * R * deltak);
null_x = mr.makeTrapezoid('x', sysCrusher, 'Area', gRO.area / 2);

% --- Gradient spoiling (at the end of TR)
spoilDur = 4e-3;
gROSpoil = mr.makeTrapezoid('x', 'Area', 6 * Nx * deltak(1), 'system', sysCrusher, 'Duration', spoilDur);
gacqSpoil = mr.makeTrapezoid('y', 'Area', 6 * Nx * deltak(1), 'system', sysCrusher, 'Duration', spoilDur);
gzSpoil = mr.makeTrapezoid('z', 'Area', 6 * Nx * deltak(1), 'system', sysCrusher, 'Duration', spoilDur);
gzS = gzSpoil.amplitude;
gxS = gROSpoil.amplitude;
gyS = gacqSpoil.amplitude;

%% 9. TIMING CALCULATIONS (TE AND DELAYS)
% Calculate timing delays to achieve the target Echo Time (TE)
rf2RO = (blip_dur / 2 + mr.calcDuration(gzCrush) - gref.delay - gref.riseTime);
minus = gz.flatTime / 2 + gz.fallTime + mr.calcDuration(gzRe) + 2 * mr.calcDuration(gNavPre) ...
      + 3 * mr.calcDuration(gRO) + mr.calcDuration(gzCrush1) + gref.delay + gref.riseTime;
delTau0 = round((rf2RO - minus) / sys.gradRasterTime) * sys.gradRasterTime;

% Calculate the final TE
TE = 2 * (gz.flatTime / 2 + gz.fallTime + mr.calcDuration(gzRe) + delTau0 ...
        + mr.calcDuration(gzCrush1, gxCrush1) + gref.riseTime + rfref_dur / 2);

% Apply the calculated delay to the crusher gradients
gxCrush1.delay = gxCrush1.delay + delTau0;
gzCrush1.delay = delTau0;

%% 10. FAT SATURATION PULSE
sat_ppm = -3.45; % Chemical shift of fat
rf_fs = mr.makeGaussPulse(110 * pi / 180, 'system', sys, 'Duration', 8e-3, ...
    'bandwidth', abs(sat_ppm * 1e-6 * sys.B0 * sys.gamma), ...
    'freqPPM', sat_ppm, 'use', 'saturation');
% Compensate for frequency-offset induced phase
rf_fs.phasePPM = -2 * pi * rf_fs.freqPPM * rf_fs.center;
gz_fs = mr.makeTrapezoid('z', sysCrusher, 'delay', mr.calcDuration(rf_fs), 'Area', 1 / 1e-4); % Spoiler

%% 11. SLICE ORDERING
% Calculate slice positions and reorder for interleaved acquisition
if numSlices > 1
    % Calculate linearly spaced positions centered at 0
    pos = (sliceThickness + spacing) * ((1:numSlices) - 1 - (numSlices - 1) / 2);
    newPos = zeros(size(pos));
    
    % Check if the number of slices is odd or even
    if mod(numSlices, 2) == 1
        % --- ODD NUMBER LOGIC ---
        % 1. Find the center slice index
        centerIndex = ceil(numSlices / 2);
        % 2. Put the center slice (position 0) first
        newPos(1) = pos(centerIndex);
        % 3. Create left and right vectors WITHOUT the center slice
        l = pos(1 : centerIndex - 1);
        r = pos(centerIndex + 1 : end);
        % 4. Apply interleaving to the rest
        newPos(2:2:end) = r; % Fills even indices (2, 4, ...)
        newPos(3:2:end) = l; % Fills odd indices (3, 5, ...)
    else
        % --- EVEN NUMBER LOGIC (Standard Interleaving) ---
        l = pos(1 : length(pos) / 2);
        r = pos(length(pos) / 2 + 1 : end);
        newPos(1:2:end) = r;
        newPos(2:2:end) = l;
    end
else
    % If only one slice, position is 0
    newPos = 0;
end

%% 12. DIFFUSION GRADIENT DEFINITIONS
% These gradients combine crusher and diffusion weighting.
% If gDiff variables (x_gDiff, etc.) are 0, they act as simple crushers.

% --- X-axis diffusion gradients ---
xDiff1 = mr.makeTrapezoid('x', sysDiff, 'Area', (mr.calcDuration(gxCrush1) / 2) * x_gDiff * sys.gamma);
xDiff2 = mr.makeTrapezoid('x', sysDiff, 'Area', -distributeDiffArea * xDiff1.area + gxCrush1.area);
xDiff2.delay = mr.calcDuration(xDiff1) + mr.calcDuration(gxCrush1) - (mr.calcDuration(xDiff1) + mr.calcDuration(xDiff2)) + addDiffDur;
xDiff_first = mr.addGradients({xDiff1, xDiff2}, sysDiff); % First diffusion block (before 180)

xDiff3 = mr.makeTrapezoid('x', sysDiff, 'Duration', mr.calcDuration(gzCrush1) - gzCrush1.delay + addDiffDur, 'Area', (xDiff1.area + xDiff2.area + 2 * gxCrush2.area), 'Delay', gxCrush2.delay);
xDiff_second = mr.addGradients({dummyx, xDiff3}, sysDiff); % Second diffusion block (after 180)

% --- Z-axis diffusion gradients ---
zDiff1 = mr.makeTrapezoid('z', sysDiff, 'Area', (mr.calcDuration(gxCrush1) / 2) * z_gDiff * sys.gamma);
zDiff2 = mr.makeTrapezoid('z', sysDiff, 'Area', -distributeDiffArea * zDiff1.area + gzCrush1.area);
zDiff2.delay = mr.calcDuration(zDiff1) + mr.calcDuration(gzCrush1) - (mr.calcDuration(zDiff1) + mr.calcDuration(zDiff2)) + addDiffDur;
zDiff_first = mr.addGradients({zDiff1, zDiff2}, sysDiff);

zDiff3 = mr.makeTrapezoid('z', sysDiff, 'Duration', mr.calcDuration(gzCrush1) - gzCrush1.delay + addDiffDur, 'Area', (zDiff1.area + zDiff2.area), 'Delay', gxCrush2.delay);
zDiff_second = mr.addGradients({dummyz, zDiff3}, sysDiff);

% --- Y-axis diffusion gradients ---
yDiff1 = mr.makeTrapezoid('y', sysDiff, 'Area', (mr.calcDuration(gxCrush1) / 2) * y_gDiff * sys.gamma);
yDiff2 = mr.makeTrapezoid('y', sysDiff, 'Area', -yDiff1.area);
yDiff2.delay = mr.calcDuration(yDiff1) + mr.calcDuration(gzCrush1) - (mr.calcDuration(yDiff1) + mr.calcDuration(yDiff2)) + addDiffDur;
yDiff_first = mr.addGradients({yDiff1, yDiff2}, sysDiff);

%% 13. SEQUENCE LOOP (BUILDING BLOCKS)
for av = 1:avg
    for slice = 1:numSlices
        % --- Start of TR ---
        seq.addBlock(rf_fs, gz_fs); % 1. Fat saturation
        
        % Set slice-specific RF and ADC properties
        adc.phaseOffset = mod(slice, 2) * pi;    % Alternate ADC phase
        adc_first.phaseOffset = mod(slice, 2) * pi;
        rf.freqOffset = gz.amplitude * newPos(slice); % Set slice frequency offset
        rf.phaseOffset = pi / 2 - 2 * pi * rf.freqOffset * mr.calcRfCenter(rf); % Phase correction
        
        % --- Excitation ---
        seq.addBlock(rf, gz);        % 2. Excitation pulse
        seq.addBlock(gzRe);          % 3. Slice rephasing
        
        % --- Navigator Echoes (pre-180 pulse) ---
        seq.addBlock(gNavPre);
        seq.addBlock(gRO, adc);
        gRO = mr.scaleGrad(gRO, -1);
        seq.addBlock(gRO, adc);
        gRO = mr.scaleGrad(gRO, -1);
        seq.addBlock(gRO, adc);
        seq.addBlock(gNavPre);
        
        % --- First Diffusion/Crusher Block ---
        seq.addBlock(xDiff_first, zDiff_first, yDiff_first);
        
        % --- Refocusing ---
        seq.addBlock(rfref, gy, xDiff_second, zDiff_second);
        
        % --- b-value waveform calculation (run once) ---
        if slice == 1 && av == 1
            [bx,by,bz]=calc_bval(seq,sys);
        end
        
        % --- EPI Readout Loop ---
        for i = 1:Ny
            % Define current and next y-blips
            gacq_current = mr.makeTrapezoid('y', sys, 'Area', rs(i), 'Duration', blip_dur);
            gacq_next = mr.makeTrapezoid('y', sys, 'Area', rs(i + 1), 'Duration', blip_dur);
            
            % Split blips for alignment
            gacq_parts_current = mr.splitGradientAt(gacq_current, blip_dur / 2, sys);
            gacq_parts_next = mr.splitGradientAt(gacq_next, blip_dur / 2, sys);
            
            if i == 1
                % First line: Align end of gRO1 with start of first blip-half
                [gacq_full, gacq_blipup, ~] = mr.align('left', gacq_current, 'right', gacq_parts_next(1), gRO1);
                gacq_first = mr.addGradients({gacq_full, gacq_blipup}, sys);
            else
                % Middle lines: Align end of last blip-half with start of next blip-half
                [gacq_blipdown, gacq_blipup, ~] = mr.align('left', gacq_parts_current(2), 'right', gacq_parts_next(1), gRO);
                gacq_blipdownup = mr.addGradients({gacq_blipdown, gacq_blipup}, sys);
            end
            
            % Add readout blocks
            if i == 1
                seq.addBlock(gRO1, gacq_first, adc_first); % First line
            elseif i == Ny
                seq.addBlock(gRO, gacq_blipdown, adc); % Last line
            else
                seq.addBlock(gRO, gacq_blipdownup, adc); % Middle lines
            end
            
            gRO = mr.scaleGrad(gRO, -1); % Flip readout polarity
        end
        
        % --- End of Readout ---
        seq.addBlock(restoreRF); % Fast recovery pulse
        seq.addBlock(null_x, null_y); % Gradient moment nulling
        
        % --- Spoiling ---
        sp = 0.9; % Spoiler amplitude scaling
        gROSpoil.amplitude = gROSpoil.amplitude * sp;
        gacqSpoil.amplitude = gacqSpoil.amplitude * sp;
        gzSpoil.amplitude = gzSpoil.amplitude * sp;
        seq.addBlock(gROSpoil, gacqSpoil, gzSpoil);
        
        % --- TR Delay ---
        if slice == 1 && av == 1
            % Calculate TR delay once
            delayTR = TR - (seq.duration);
        end
        seq.addBlock(mr.makeDelay(delayTR)); % Add delay to fill TR
    end
end

%% 15. SEQUENCE CHECKS
% Check if SPEN conditions are met (full refocusing)
% 1. 2 * RF refocusing duration == total acquisition duration
assert(abs(round(Ny * mr.calcDuration(gRO) / 2, 6) - round(rfref_dur, 6)) < 1e6)
% 2. 2 * SPEN gradient area == total y-blip area
assert(abs(gref.flatArea * 2) - abs(sum(rs)) < 1e-3)

% Run Pulseq checker
name = 'MS_SPEN_SE_EPI';
rep = check(seq, TE, TR, delayTR);
fprintf([rep{:}]);

%%
% Optional: Plot the RF pulse profile
[rf_bw, ~, rf_spectrum, rf_w] = mr.calcRfBandwidth(rfref);
figure; plot(rf_w, abs(rf_spectrum));
title('Refocusing pulse profile (small-angle approximation)');
xlabel('Frequency [Hz]');
xlim(2 * [-rf_bw, rf_bw]);
grid on

%% 16. SET DEFINITIONS (Metadata for .seq file)
seq.setDefinition('Name', name);
seq.setDefinition('FOV', fov);
seq.setDefinition('NumSlices', numSlices);
seq.setDefinition('Ny', Ny);
seq.setDefinition('Nx', Nx);
seq.setDefinition('ReceiverGainHigh', 1);
seq.setDefinition('TE', TE);
seq.setDefinition('TR', TR);
seq.setDefinition('DurAcq', durAcq);
seq.setDefinition('RO_dur', roDur);
seq.setDefinition('R', R);
seq.setDefinition('sweepBw', bw);
seq.setDefinition('rfref_dur', rfref_dur);
seq.setDefinition('Gspen', gref.amplitude);
seq.setDefinition('SeqDur', seq.duration());
seq.setDefinition('rampSampling', rampSampling);
seq.setDefinition('ro_os', ro_os);
seq.setDefinition('randomSampling', randomSampling);
seq.setDefinition('bval', [bx, by, bz]);
seq.setDefinition('avg', avg)

%% 17. PNS CALCULATION AND EXPORT
% Calculate Peripheral Nerve Stimulation (PNS)
% Replace with your system's gradient file path
% [pns_ok, pns_n, pns_c, tpns] = seq.calcPNS('*.asc');

% --- Logic for dynamic filename ---
b_str = ''; % Initialize diffusion string
if x_gDiff ~= 0
    b_str = strcat('_bx',int2str(round(bx)));
end
if y_gDiff ~= 0
    b_str = strcat('_by',int2str(round(by)));
end
if z_gDiff ~= 0
    b_str = strcat('_bz',int2str(round(bz)));
end

% If no diffusion gradient was active, mark it as b0
if isempty(b_str)
    b_str = '_b0';
end

filename = sprintf('SPEN_R%d_avg%d%s.seq', R, avg, b_str);

% path='..\';

% --- Write the sequence file ---
% seq.write(strcat(path,filename));


%%
simulateSPENfo.spins_per_voxel = 10;

% simulateSPENfo=[];

deriveSPENforwardOperator('ref',seq,sys,simulateSPENfo)





function [bx,by,bz]=calc_bval(seq,sys)
    
    % Calculate k-space and get refocusing time
    [~, ~, ~, ~, kspace.t_excitation, kspace.t_refocusing] = seq.calculateKspacePP();
    % Get all waveforms
    [gw, rfw, adcw] = seq.waveforms_and_times();
    
    % Process X-gradient waveform
    xwf(1, :) = [gw{1, 1}(1, 1) - sys.gradRasterTime, gw{1, 1}(1, 1):sys.gradRasterTime:gw{1, 1}(1, end)];
    xwf(2, :) = [0, mr.pts2waveform(gw{1, 1}(1, :), gw{1, 1}(2, :), sys.gradRasterTime), 0];
    xwf(2, find(abs(xwf(1, :) - kspace.t_refocusing) < 1e-6):end) = -1 * xwf(2, find(abs(xwf(1, :) - kspace.t_refocusing) < 1e-6):end); % Flip post-180
    xwf(2, find(xwf(1, :) <= kspace.t_excitation)) = 0; % Zero before excitation
    
    % Process Y-gradient waveform
    ywf(1, :) = [gw{1, 2}(1, 1) - sys.gradRasterTime, gw{1, 2}(1, 1):sys.gradRasterTime:gw{1, 2}(1, end)];
    ywf(2, :) = [0, mr.pts2waveform(gw{1, 2}(1, :), gw{1, 2}(2, :), sys.gradRasterTime), 0];
    ywf(2, find(abs(ywf(1, :) - kspace.t_refocusing) < 1e-6):end) = -1 * ywf(2, find(abs(ywf(1, :) - kspace.t_refocusing) < 1e-6):end);
    ywf(2, find(ywf(1, :) <= kspace.t_excitation)) = 0;
    
    % Process Z-gradient waveform
    zwf(1, :) = [gw{1, 3}(1, 1) - sys.gradRasterTime, gw{1, 3}(1, 1):sys.gradRasterTime:gw{1, 3}(1, end)];
    zwf(2, :) = [0, mr.pts2waveform(gw{1, 3}(1, :), gw{1, 3}(2, :), sys.gradRasterTime), 0];
    zwf(2, find(abs(zwf(1, :) - kspace.t_refocusing) < 1e-6):end) = -1 * zwf(2, find(abs(zwf(1, :) - kspace.t_refocusing) < 1e-6):end);
    zwf(2, find(zwf(1, :) <= kspace.t_excitation)) = 0;

    % Calculate b-values from the processed gradient waveforms
    bx = (sum((cumsum(xwf(2, :)) * sys.gradRasterTime).^2) * sys.gradRasterTime * (2 * pi)^2) / 1e6;
    by = (sum((cumsum(ywf(2, :)) * sys.gradRasterTime).^2) * sys.gradRasterTime * (2 * pi)^2) / 1e6;
    bz = (sum((cumsum(zwf(2, :)) * sys.gradRasterTime).^2) * sys.gradRasterTime * (2 * pi)^2) / 1e6;

end
