% GNU_Radio_Postprocess_PDPs_Doppler
%
% Authors: Professor Chris Anderson
%          CDR Thomas W. Tedesso
%
% Program inputs a data file from software defined radio receiver and a
% PN-sequence used in a DSS signal used for channel measurments.  The
% received signal is correlated with the pn-sequence to determine the delay
% spread and doppler characteristics of the channel.  The program processes
% segments of the input file (ten thousand pn-sequences) and produces as an
% output a 3D image of the output magnitude vs time vs delay spread.  It
% also produces a 3D image of magnitude vs delay spread vs doppler spread.
% 
%

clear all
close all
set(0, 'DefaultAxesFontSize', 14)
set(0, 'DefaultAxesFontWeight','Bold')
fs = input('enter the sampling frequency in MHz >>');
fs = fs*1e6;
Ts = 1/fs;

pnseq = csvread('pntest.txt');
%
pnseq = 2*pnseq - 1;  % Make pn seq an antipodal sequence
%
pnlength = length(pnseq);  % Length of the PN Sequence in Chips
%
pntime = 0:Ts:Ts*(pnlength-1);

% Plot PN Sequence 
figure(1)
plot(pntime*1e6,pnseq)
xlabel('Time (us)')
ylabel('Amplitude (V)')

% Load the captured data file 
%
%filename = 'pn_doppler_test_1.dat';
filename = '2015.09.25.14.56.16.dat';
datafile = fopen(filename, 'rb');
%
% Read and process data file in blocks
numpnseq = 1e4;             % number of pn-sequence blocks to process
segsize = pnlength*numpnseq;            % segment size to process
%
% Process each data segment.  Check for end of file at each iteration.
counter = 0;
while ~feof(datafile)
    counter = counter+1                         % Tracks number of iterations
    z = fread (datafile, [2,segsize], 'float');
    data = z(1,:) + z(2,:)*1i; 
    [r, c] = size (data); 
    data = reshape (data, c, r);
    
    % Note: Originally we had two samples/chip, probably only need 1 for the
    % final implementation
    ds_data = downsample(data,2);
    % Now, we need to handle correlating the Inphase and Quadrature in orger to
    % generate the Magnitude and Phase PDP's
    i_data = real(ds_data);
    q_data = imag(ds_data);
    % Run the cross-correlation, generate the time vector, and plot the results
    % Note that the time vector starts at an offset of the PN length
    pdp_i = xcorr(i_data,pnseq');
    pdp_q = xcorr(q_data,pnseq'); 
    L = length(pdp_i);
    % 
    % Since the first L/2 -1 correlation results = 0, these values are
    % discarded to preserve memory and minimize calculations
    %
    pdp_i = pdp_i(floor(L/2):L);
    pdp_q = pdp_q(floor(L/2):L);
    pdp_mag = abs(pdp_i+1i*pdp_q);
    pdp_mag_db = mag2db(pdp_mag);
    pdp_phase = atan2d(pdp_q,pdp_i); 
    pdp_time = 0:Ts:(L/2+1)*Ts;
    %
    % Display magnitude and phase plots
    %
    figure(2)
    subplot(2,1,1)
    plot(pdp_time*1e6,pdp_mag_db,'k-')
    xlabel('Time (us)')
    ylabel('Relative Power (dB)')
    axis([1000 4000 0 80])
    subplot(2,1,2)
    plot(pdp_time*1e6,pdp_phase,'r-')
    xlabel('Time (us)')
    ylabel('Phase (degrees)')
    axis([1000 4000 -180 180])
    %
    % This part grabs a sample index on the PDP and generates the mag/phase
    % data for Doppler Spectrum generation
    %
    M = 3;                                     % Multiplier for Threshold Detector
    min_peak_height = M*mean(pdp_mag);          % Threshold 
    [pdp_abs_pk, pdp_start_index] = findpeaks(pdp_mag,'MinPeakHeight', min_peak_height, 'NPeaks', 50); % determine the peaks
    start = pdp_start_index<=pdp_start_index(1)+pnlength-1;         % find the start index of the peaks (must be within pnseq)
    pdp_start_index = pdp_start_index(start);
    sample_increment = pnlength;  % assuming 1 sample per chip
    pdp_stop_index = pdp_start_index + sample_increment.*floor((length(pdp_mag)-pdp_start_index)/sample_increment);
    %
    % Preallocate matrix to store individual doppler data
    %
    CC = ceil((length(pdp_mag)-pdp_start_index(1))/sample_increment);
    RR = length(pdp_start_index);
    pdp_doppler_i = zeros(RR, CC);             
    pdp_doppler_q = zeros(RR,CC);
    for d = 1:RR
        pdp_doppler_i(d,:) = pdp_i(pdp_start_index(d):sample_increment:pdp_stop_index(d));
        pdp_doppler_q(d,:) = pdp_q(pdp_start_index(d):sample_increment:pdp_stop_index(d));
    end
    fs_doppler = fs/(2*sample_increment);
    time_doppler = 0:1/fs_doppler:(1/fs_doppler)*(length(pdp_doppler_i)-1);
    %
    % Plot the time domain signal for the doppler shifted signals
    %
    for d=1:RR
        figure
        plot(time_doppler*1e3,pdp_doppler_i(d,:),'r-');
        hold on
        plot(time_doppler*1e3,pdp_doppler_q(d,:),'b-');
        xlabel('Time (ms)')
        ylabel('Amplitude (V)')
    end
    %
    % Graph time delay spread
    %
    % Reshape matrix of magnitude of correlation vector into a matrix to allow
    % graphing magnitude vs delay spread vs time.  The intial pnlength-1 index
    % points are disregarded due to them indicating negative lags.     
    %
    R = length(pnseq);
    C = floor(length(pdp_mag_db)/R);
    pdp_mag_db_matrix = reshape(pdp_mag_db(1:R*C), R, C); 
    figure
    delay = (0:2/fs:(R-1)*2/fs)*1e6;
    time = (0:R*2/fs:(C-1)*(R*2)/fs) *1000;
    surf(time,delay, pdp_mag_db_matrix);
    ylabel('\mu sec');
    xlabel('msec');
    zlabel('dB');
    shading interp
    %
    % Graph Doppler spectrum for each time delayed signal
    %
    K = 1024;                % size of fft;
    PDP_DOPPLER = fft(pdp_doppler_i.',K);  
    f = ((-K/2:K/2-1)*fs_doppler/K)*0.001;
    for m = 1:RR
        figure
        pwelch(pdp_doppler_i(m,:),512,256,512,fs_doppler);
        %plot(f,fftshift(mag2db(abs(PDP_DOPPLER(:,m)))));
        %xlabel('kHz');
        %ylabel('db');
    end

    R = length(pnseq);
    C = floor(length(pdp_i)/R);
    pdp_i_matrix = reshape(pdp_i(1:R*C), R, C); 

    PDP_I_MATRIX = fft(pdp_i_matrix',K);
    figure
    surf(delay,f,mag2db(fftshift(abs(PDP_I_MATRIX))));
    shading interp
    xlabel('\mu sec')
    ylabel('kHz');
end
fclose(datafile);