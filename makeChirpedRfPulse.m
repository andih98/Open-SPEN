%% Open SPEN in Pulseq - Chirp RF pulse function 
% Andreas Holl
% Division of Medical Physics, Department of Diagnostic and Interventional Radiology,
% University Medical Center Freiburg, Faculty of Medicine, University of Freiburg, Freiburg, Germany
% Email: andreas.holl@uniklinik-freiburg.de
% July. 23, 2025
% A part of the following code was taken from 
% https://github.com/mikgroup/sigpy/tree/main/sigpy/mri.

function rf = makeChirpedRfPulse(varargin)
    % --- Parse input parameters ---
    validPulseUses = mr.getSupportedRfUse();
    persistent parser
    if isempty(parser)
        parser = inputParser;
        parser.FunctionName = 'makeChirpedRfPulse2';
        
        % RF Parameters
        addOptional(parser, 'system', mr.opts(), @isstruct);
        addParamValue(parser, 'duration', 0, @isnumeric);
        addParamValue(parser, 'ang', 0, @isnumeric); % Flip angle in degrees
        addParamValue(parser, 'freqOffset', 0, @isnumeric);
        addParamValue(parser, 'freqPPM', 0, @isnumeric);
        addParamValue(parser, 'phasePPM', 0, @isnumeric);
        addParamValue(parser, 'phaseOffset', 0, @isnumeric);
        addParamValue(parser, 'n_fac', 40, @isnumeric); % Shape factor for amplitude modulation
        addParamValue(parser, 'bandwidth', 0, @isnumeric);
        addParamValue(parser, 'adiabaticity', 0.4, @isnumeric);
        
        addParamValue(parser, 'delay', 0, @isnumeric);
        addParamValue(parser, 'dwell', 0, @isnumeric);
        addOptional(parser, 'use', '', @(x) any(validatestring(x,validPulseUses)));
    end
    parse(parser, varargin{:});
    opt = parser.Results;

    if opt.dwell == 0
        opt.dwell = opt.system.rfRasterTime;
    end

    Nraw = round(opt.duration / opt.dwell + eps);
    % Number of points must be divisible by 4 (library/hardware requirement)
    N = floor(Nraw / 4) * 4; 
    
    % Time vector from 0 to pulse duration
    t = linspace(0, opt.duration, N);
    
    % --- Pulse design: Amplitude and Phase ---
    
    % 1. Amplitude modulation (AM)
    am = 1 - abs(cos(pi * t / opt.duration)).^opt.n_fac;

    % 2. Frequency modulation (FM) - linear chirp
    fm = linspace(-opt.bandwidth / 2, opt.bandwidth / 2, N) * 2 * pi;

    % 3. Phase evolution (PM) by integrating the frequency
    % This is the correct way to obtain the phase for a complex signal.
    pm = cumsum(fm) * opt.dwell;

    % We find the pulse center to center the phase and to determine the 
    % rate of frequency change (roc_fm0) for the amplitude calculation.
    [~, ifm] = min(abs(fm)); % Find index closest to frequency=0
    
    if fm(ifm) == 0
        pm0 = pm(ifm);
        am0 = am(ifm);
        roc_fm0 = abs(fm(ifm + 1) - fm(ifm - 1)) / (2 * opt.dwell);
    else % If 0 is not hit exactly, interpolate
        if fm(ifm) * fm(ifm+1) < 0, b=1; else, b=-1; end
        pm0 = (pm(ifm)*fm(ifm+b) - pm(ifm+b)*fm(ifm)) / (fm(ifm+b)-fm(ifm));
        am0 = (am(ifm)*fm(ifm+b) - am(ifm+b)*fm(ifm)) / (fm(ifm+b)-fm(ifm));
        roc_fm0 = abs(fm(ifm) - fm(ifm+b)) / opt.dwell;
    end
    
    % Center the phase correctly
    pm = pm - pm0;

    % 5. Calculate amplitude factor 'a' based on adiabaticity
    % This is the crucial step that makes the pulse robust.
    a = ((roc_fm0 * opt.adiabaticity)^0.5) / (2 * pi * am0);

    % 6. Generate complex RF signal
    signal = a * am .* exp(1i * pm);

    % --- Scaling for the desired flip angle ---
    % The scaling factors are empirical to compensate for the reduced efficiency
    % of FM pulses compared to pure AM pulses.
    % *1.9 for pi, 1.04 for pi/2 from D.Kunz1986 
    % https://doi.org/10.1002/mrm.1910030303
    % --> The integrated area of an FM RF pulse must be greater
    % than that of a constant-frequency amplitude-modulated
    % pulse to obtain the same tip angle. This is because the FM
    % pulse is on-resonance for any given spin for only a frac-
    % tion of the duration of that pulse.
    if opt.ang == 180
        signal = signal * ((opt.ang / 360) / abs(sum(signal) * opt.dwell)) * 2.1;
    elseif opt.ang == 90
        signal = signal * ((opt.ang / 360) / abs(sum(signal) * opt.dwell)) * 1.07;
    else
        % General scaling for other angles (without empirical factor)
        if sum(signal) ~= 0
           signal = signal * (opt.ang * pi/180) / (2 * pi * sum(signal) * opt.dwell);
        end
    end
    
    % If N was rounded down, pad with zeros
    if (N ~= Nraw)
        Npad = Nraw - N;
        signal = [zeros(1, floor(Npad/2)) signal zeros(1, Npad-floor(Npad/2))];
        N = Nraw;
        t = (0:N-1) * opt.dwell;
    end
    
    % --- Create Pulseq RF structure ---
    rf.type = 'rf';
    rf.signal = signal;
    rf.t = t;
    rf.shape_dur = N * opt.dwell;
    rf.freqOffset = opt.freqOffset;
    rf.phaseOffset = opt.phaseOffset;
    rf.deadTime = opt.system.rfDeadTime;
    rf.ringdownTime = opt.system.rfRingdownTime;
    rf.delay = opt.delay;
    rf.center = length(signal)/2*opt.dwell;
    rf.freqPPM=opt.freqPPM;
    rf.phasePPM=opt.phasePPM;
    if ~isempty(opt.use)
        rf.use = opt.use;
    else
        rf.use = 'excitation';
    end
    
    if rf.deadTime > rf.delay
        rf.delay = rf.deadTime;
    end
end