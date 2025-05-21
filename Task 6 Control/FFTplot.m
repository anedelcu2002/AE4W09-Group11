% --- FFT and Plotting Script for Fore-Aft Acceleration and Pitch ---

% Inputs (replace with your actual data)
% t       - time vector (s)
% acc_FA  - fore-aft acceleration (m/s^2)
% pitch   - pitch angle or rate (deg or deg/s)

% Sampling information
Fs = 1 / mean(diff(Time));     % Sampling frequency
N = length(Time);              % Number of samples
f = (0:N-1)*(Fs/N);         % Frequency vector

% Perform FFT
acc_FA_fft = abs(fft(NcIMUTAxs)) / N;
pitch_fft  = abs(fft(BlPitch1)) / N;

% Keep only first half of spectrum (real signals are symmetric)
half_idx = 1:floor(N/2);
f = f(half_idx);
acc_FA_fft = 2 * acc_FA_fft(half_idx);
pitch_fft  = 2 * pitch_fft(half_idx);

% Plotting
figure;
subplot(2,1,1);
plot(f, acc_FA_fft, 'b', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('FFT of Fore-Aft Acceleration for 0.0625 gain');
grid on;

subplot(2,1,2);
plot(f, pitch_fft, 'r', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('FFT of Pitch for 0.0625 gain');
grid on;

% Compute area under the amplitude spectra
acc_area = trapz(f, acc_FA_fft);
pitch_area = trapz(f, pitch_fft);

% Display the results
fprintf('Spectral Area (Fore-Aft Acceleration) for 0.0625 gain: %.4f\n', acc_area);
fprintf('Spectral Area (Pitch) for 0.0625 gain: %.4f\n', pitch_area);
