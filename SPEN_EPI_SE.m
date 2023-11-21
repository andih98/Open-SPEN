clear all
close all
clc


% Define FOV and resolution
fov = 256e-3;
Nx = 100;
Ny = Nx;

% Define other gradients and ADC events
deltak = 1/fov; % Pulseq toolbox defaults to k-space units of m^-1
kWidth = Nx*deltak;

% Define sequence parameters
TE = 40e-3;

% set system limits
sys = mr.opts('MaxGrad',40,'GradUnit','mT/m',...
    'MaxSlew',200,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6,'B0',3);

% Create a new sequence object
seq=mr.Sequence(sys);

% Create Gexc
rf_dur=5.85e-3; 
gexc=mr.makeTrapezoid('y',sys,'Amplitude',35e-3*sys.gamma, 'FlatTime', rf_dur);
gexcRe=mr.makeTrapezoid('y',sys,'Area',(gexc.area-deltak*Ny)/2);

% Calculate RF-bw
sweepBw=gexc.amplitude*fov*20;

% Create a chirped 90°RF pulse
rf=makeChirpedRfPulse('wurst','duration',rf_dur,'delay',gexc.riseTime,'bandwidth',sweepBw, ...
    'ang',90,'n_fac',40,'system',sys); 

% Create a slice selective sinc 180° RF pulse, SS and crusher
dur_ref=6e-3;
sliceThickness=5e-3;
[rfref, gz] = mr.makeSincPulse(pi,sys,'Duration',dur_ref,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',0,'use','refocusing');

gzCrush = mr.makeTrapezoid('z', 'Amplitude', -gz.amplitude, 'Duration', (dur_ref)/2);

% Create readout gradient, crusher and adc
dwell=6e-5;
gRO = mr.makeTrapezoid('x',sys,'FlatArea',kWidth,'FlatTime',30e-5);%rf_dur

adc = mr.makeAdc(Nx,sys,'Duration',gRO.flatTime,'Delay',gRO.riseTime);

gxCrush1 = mr.makeTrapezoid('x',sys, 'Area', gRO.area/4);
gxCrush2 = mr.makeTrapezoid('x',sys, 'Area', -gRO.area/4);

% Create acquisition gradient
acqDur = ceil(2*sqrt(deltak/sys.maxSlew)/10e-6)*10e-6;
gacq = mr.makeTrapezoid('y',sys,'Area',deltak,'duration',acqDur);

% Delay
delay1=round((TE/2-mr.calcDuration(gexc)/2-mr.calcDuration(gzCrush)-mr.calcDuration(rfref)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delay2=round((TE/2-mr.calcDuration(rfref)/2-mr.calcDuration(gzCrush)-mr.calcDuration(gexcRe)-mr.calcDuration(gRO)/2)/seq.gradRasterTime)*seq.gradRasterTime;

% Create Blocks and loop over acqusition
seq.addBlock(rf,gexc)

seq.addBlock(mr.makeDelay(delay1))
seq.addBlock(gzCrush,gxCrush1)
seq.addBlock(rfref,gz)
seq.addBlock(gzCrush,gxCrush2)
seq.addBlock(gexcRe)

seq.addBlock(mr.makeDelay(delay2))
for i=1:Ny
    seq.addBlock(gRO,adc);            % Read one line of k-space
    seq.addBlock(gacq);               % Phase blip
    gRO.amplitude = -gRO.amplitude;   % Reverse polarity of read gradient
end
seq.addBlock(mr.makeDelay(1));
seq.install('siemens')
%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% export and visualization
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'SPEN_EPI_14_11_23');
%seq.write('SPEN_SE_EPI.seq');   % Output sequence for scanner
seq.plot();             % Plot sequence waveforms

%% calculate trajectory 
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
%[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace();

%% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis

figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display

%% sanity checks
TE_check=(t_refocusing(1)-t_excitation(1))*2;
fprintf('intended TE=%.03f ms, actual spin echo TE=%.03fms\n', TE*1e3, TE_check*1e3); 

rep = seq.testReport;
fprintf([rep{:}]);
return
%% EPI reconstruction
% read in rawdata
[scan_ID,path] = uigetfile('*.dat','select the rawdata*.dat file','C:/Users/andre/Documents/MATLAB/MR_data/161123');
s = convertStringsToChars(strcat(path,scan_ID));
raw=mapVBVD(s);
rawdata=permute(raw.image{''},[3 1 2]);

% eliminate linear phase shift
rawdata(:,1:end-1,:)=rawdata(:,2:end,:);

% flip even k-space lines because of EPI
for k=1:size(rawdata,3)
for i=2:2:length(rawdata)
    rawdata(i,:,k)=fliplr(rawdata(i,:,k));
end
end
figure; imab(rawdata)

% eliminate ghosts
rawdata1=zeros(size(rawdata));
rawdata2=rawdata1;
rawdata1(1:2:end,:,:)=rawdata(1:2:end,:,:);
rawdata2(2:2:end,:,:)=rawdata(2:2:end,:,:);

imdata1=fftshift(fft(fftshift(rawdata1,1),[],1),1);
imdata2=fftshift(fft(fftshift(rawdata2,1),[],1),1);
for i=1:size(imdata1,3)
    reco1(:,:,i)=new_fresnel_conv(squeeze(imdata1(:,:,i)),[0,1],[2.3,0],[20/128,0]);
end
for i=1:size(imdata2,3)
    reco2(:,:,i)=new_fresnel_conv(squeeze(imdata2(:,:,i)),[0,1],[2.3,0],[20/128,0]);
end

const_phase_diff_cmplx=sum(reco2(1,:).*conj(reco1(1,:)));
const_phase_diff_cmplx=const_phase_diff_cmplx/abs(const_phase_diff_cmplx);

for i=1:size(reco1,3)
    image(:,:,i)=abs(reco1(:,:,i)-conj(const_phase_diff_cmplx)*reco2(:,:,i));
end

image=sqrt((sum(abs(image(:,:,1:size(image,3))).^2,3)));
figure; imab(imrotate(image,-90))
colormap('gray')
return
%% 2D FFT reconstruction
% clear all
% close all
% clc

[scan_ID,path] = uigetfile('*.dat','select the rawdata*.dat file','C:/Users/andre/Documents/MATLAB/MR_data/141123');
s = convertStringsToChars(strcat(path,scan_ID)) ;
raw=mapVBVD(s);
rawdata=permute(raw.image{''},[3 1 2]);
imdata=conj(fftshift(fft(fftshift(rawdata,1),[],1),1));
for i=1:size(imdata,3)
    reco(:,:,i)=new_fresnel_conv(squeeze(conj(imdata(:,:,i))),[0,1],[2.3,0],[20/128,0]);
end
image=sqrt((sum(abs(reco(:,:,1:size(reco,3))).^2,3)));
figure; imab(imrotate(image,-90))    
colormap('gray')
figure; imab(imrotate(angle(reco(:,:,1)),-90))