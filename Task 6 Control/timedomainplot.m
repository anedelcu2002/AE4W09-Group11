% --- FFT and Plotting Script for Fore-Aft Acceleration and Pitch ---

% Inputs (replace with your actual data)
% t       - time vector (s)
% acc_FA  - fore-aft acceleration (m/s^2)
% pitch   - pitch angle or rate (deg or deg/s)

% Plotting
figure;
subplot(2,1,1);
plot(Time, NcIMUTAxs, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Fore-Aft Acceleration for 0 gain');
grid on;

subplot(2,1,2);
plot(Time, BlPitch1, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Pitch for 0 gain');
grid on;
