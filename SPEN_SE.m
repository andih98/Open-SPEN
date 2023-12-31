clear all
close all
clc

gyrom=42.577e6;
% Define FOV and resolution
fov = 256e-3;
sliceThickness = 5e-3;
Nx = 100;
Ny = Nx;
tacq=2.56e-3;
sweepBw=500000;
gexc=(sweepBw/(gyrom*fov))*1e3;


% Define sequence parameters
TE = 9e-3;
TR = 600e-3;
alphaRef=180;

% set system limits
sys = mr.opts('MaxGrad',40,'GradUnit','mT/m',...
    'MaxSlew',200,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6);

% Create a new sequence object
seq=mr.Sequence(sys);

% Create a chirped 90°RF pulse
duration=2.5e-3; 
rf= makeWurstPulse('wurst','duration',duration,'bandwidth',sweepBw,'n_fac',40,'system',sys); 

% Create a slice selective sinc 180° RF pulse 
[rfref, gz] = mr.makeSincPulse(pi,sys,'Duration',2e-3,...
    'SliceThickness',5e-3,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',0,'use','refocusing');
gzRe = mr.makeTrapezoid('z', 'Amplitude', -gz.amplitude*3, 'Duration', (4e-3)/2);

% Define other gradients and ADC events
deltak = 1/fov; % Pulseq toolbox defaults to k-space units of m^-1
%gx = mr.makeTrapezoid('x', 'FlatArea',Nx*deltak, 'FlatTime', duration);
gx = mr.makeTrapezoid('x', 'Amplitude',gexc, 'Duration', duration);
adc = mr.makeAdc(Nx, 'Duration', gx.flatTime, 'Delay', gx.riseTime);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;
dur_gyPre=2e-3;

% Calculate timing
delayTE1 = round((TE/2-mr.calcDuration(rf)/2-mr.calcDuration(rfref)/2-dur_gyPre)/seq.gradRasterTime)*seq.gradRasterTime;
delayTE2 = round((TE/2-mr.calcDuration(rfref)/2-mr.calcDuration(gzRe)-mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR = round((TR - mr.calcDuration(rf)-dur_gyPre-delayTE1-mr.calcDuration(rfref)-mr.calcDuration(gzRe)-delayTE2-mr.calcDuration(adc))/seq.gradRasterTime)*seq.gradRasterTime;


% Loop over phase encodes and define sequence blocks
for i=1:Ny
    
    seq.addBlock(rf,gx);
    seq.addBlock(mr.makeDelay(delayTE1));
    gyPre = mr.makeTrapezoid('y', 'Area', -phaseAreas(i)/2, 'Duration', dur_gyPre);
    seq.addBlock(gyPre, gzRe);
    seq.addBlock(rfref,gz);
    gyPre = mr.makeTrapezoid('y', 'Area', phaseAreas(i)/2, 'Duration', dur_gyPre);
    seq.addBlock(gyPre, gzRe);
    seq.addBlock(mr.makeDelay(delayTE2));
    seq.addBlock(gx, adc);
    seq.addBlock(mr.makeDelay(delayTR));

end

% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (~ok)
    fprintf('Timing check failed! Sequence probably will not run on the scanner.\n'); 
end

% export definitions
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'SPEN_SE_bw500e3'); % if submitting a sequence please write your name to the Name field of the definition section

% seq.write('SPEN_SE_test.seq')       % Write to pulseq file

%seq.plot('timeRange', [0 2*TR])

seq.install('siemens');

rep = seq.testReport;
fprintf([rep{:}]);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              