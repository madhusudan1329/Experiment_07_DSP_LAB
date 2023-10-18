a = 3;
% Last 3 digits of Register number is 284
% alpha is 3
%1

% Define the parameters
N = 140;  % Window length
n = 0:N-1;  % Discrete time index

% Calculate the window functions
rectangular_window = ones(1, N);  % Rectangular window
hamming_window = 0.54 - 0.46 * cos(2 * pi * n / (N - 1));  % Hamming window
hanning_window = 0.5 * (1 - cos(2 * pi * n / (N - 1)));  % Hanning window

% Plot the window functions
subplot(3, 1, 1);
plot(n, rectangular_window);
title('Rectangular Window');
xlabel('n');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(n, hamming_window);
title('Hamming Window');
xlabel('n');
ylabel('Amplitude');

subplot(3, 1, 3);
plot(n, hanning_window);
title('Hanning Window');
xlabel('n');
ylabel('Amplitude');

sgtitle('Window Functions in MATLAB');

%2

% Define the DFT length
N = 1024;  % Length of DFT

% Define a vector of different window lengths
window_lengths = [32, 64, 128, 256];

% Initialize a matrix to store the frequency domain representations
spectrum = zeros(length(window_lengths), N);

% Calculate and plot the spectrum for each window length
figure;

for i = 1:length(window_lengths)
    L = window_lengths(i);
    
    % Generate the Blackman window of length L
    w = hamming(L);
    
    % Zero-pad the window to match the DFT length
    w_padded = [w; zeros(N - L, 1)];
    
    % Compute the DFT of the window
    W = fft(w_padded, N);
    
    % Normalize the magnitude by the actual length
    normalized_spectrum = abs(W) / L;
    
    % Store the normalized spectrum for plotting
    spectrum(i, :) = normalized_spectrum;
    
    % Plot the spectrum
    subplot(length(window_lengths), 1, i);
    plot((0:N-1) / N, normalized_spectrum);
    title(['hanning Window Spectrum (Length = ' num2str(L) ')']);
    xlabel('Normalized Frequency');
    ylabel('Magnitude');
end

sgtitle('Spectrum of Hanning Window for Different Lengths')



a = 3;
Fc = 1/3;
N = 21;
Hf = fdesign.lowpass('N,Fc',N,Fc);
Hd1 = design(Hf,'window','window',@blackman,'systemobject',true);
Hd2 = design(Hf,'window','window',@rectwin,'systemobject',true);
hfvt = fvtool(Hd1,Hd2,'Color','White');
legend(hfvt,'Hamming window design', ...
       'rectangular window design');
figure;

%2 

[num, den] = tf(Hd1);
bode(num, den);

%Time-Domain Windowing
% Create a time-domain signal
t = 0:0.001:1;
signal = sin(2*pi*50*t) + sin(2*pi*100*t);

% Apply a Hanning window to the signal
windowed_signal = signal .* hann(length(signal));

% Plot the original signal and the windowed signal
subplot(2, 1, 1);
plot(t, signal);
title('Original Signal');
xlabel('Time');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t, windowed_signal);
title('Windowed Signal (Hanning Window)');
xlabel('Time');
ylabel('Amplitude');
figure;

%Window-Based FIR filter

% Create an example signal
t = 0:0.001:1;
signal = sin(2*pi*50*t) + sin(2*pi*100*t) + 0.2*randn(size(t));

% Design a low-pass FIR filter using window-based design
filter_order = 51;
cutoff_frequency = 80; % Hz
nyquist = 0.5 * 1000; % Nyquist frequency for a 1 kHz sampling rate
normalized_cutoff = cutoff_frequency / nyquist;
filter_coeffs = fir1(filter_order, normalized_cutoff);

% Apply the FIR filter to the signal
filtered_signal = filter(filter_coeffs, 1, signal);

% Plot the original signal and the filtered signal
subplot(2, 1, 1);
plot(t, signal);
title('Original Signal');
xlabel('Time');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t, filtered_signal);
title('Filtered Signal (Low-Pass FIR Filter)');
xlabel('Time');
ylabel('Amplitude');