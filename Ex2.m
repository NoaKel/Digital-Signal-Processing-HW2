clear all; close all; clc;

% Question 1 %
% Signal DFT
load('signal.mat');
fs = 2e3;
N = length(x);
Xd = fftshift(fft(x));
f = (fs/N)*(1:(N/2)); 
Xd_DB = mag2db(abs(Xd));
Xd_DB = Xd_DB((N/2)+1:N);
% plot 1.1
figure(1); subplot(2,1,1);
plot(f,Xd_DB); 
xlabel('f [Hz]');
ylabel('Xd [dB]'); 
title('Signal DFT');
% plot 4.1
figure(4); subplot(3,1,1);
plot(f,Xd_DB); 
xlabel('f [Hz]');
ylabel('X_f [dB]'); 
title('Signal DFT');

% Question 4 %
% Divide signal to 256 and DFT each sig
Ts = 1/fs;
N_sig = 256;
M_DFT = fft_windows(x,N_sig);
M_DFT_DB = mag2db(abs(M_DFT));
M_DFT_DB = M_DFT_DB(:,(N_sig/2+1):N_sig);
m = size(M_DFT,1);
f_star = (fs/N_sig)*(1:(N_sig/2)); 
t = linspace(0,N*Ts,m); 
[f_matrix, t_matrix] = meshgrid(f_star,t);
% plot 2.1
figure(1); subplot(2,1,2);
mesh(f_matrix,t_matrix,M_DFT_DB);
xlabel('f [Hz]');
ylabel('t [sec]');
zlabel('Xd [dB]');
title('Divided Signal DFT - 3D');
view(10,70);

% Question 5 %
% Different manipulations on the sig of time - 0.55 sec
t0 = 0.55;
selected_row = ceil(t0/(N_sig*Ts));
% Plot 2.1
figure(2); subplot(2,3,1);
plot(f_star,M_DFT_DB(selected_row,:));
xlabel('f [Hz]');
ylabel('Xd [dB]'); 
title('row DFT at t=0.55 sec , row of matrix, Rec');

% Question 6 %
% Different manipulations on the sig of time - 0.55 sec
x_samples = x(((selected_row-1)*N_sig+1+N_sig/4):(selected_row*N_sig-N_sig/4));
Xd_selected_row = fftshift(fft(x_samples));
Xd_selected_row_DB = mag2db(abs(Xd_selected_row));
Xd_selected_row_DB = Xd_selected_row_DB((N_sig/4)+1:(N_sig/2));
f_selected_row = linspace(0,fs/2,N_sig/4); 
% plot 2.2
figure(2); subplot(2,3,2);
plot(f_selected_row,Xd_selected_row_DB); 
xlabel('f [Hz]');
ylabel('Xd [dB]'); 
title('row DFT at t=0.55 sec, 128 samples, Rec');

% Question 7 %
% Different manipulations on the sig of time - 0.55 sec
x_padded = [x_samples, zeros(1,128)];
Xd_padded = fftshift(fft(x_padded));
Xd_padded_DB = mag2db(abs(Xd_padded));
Xd_padded_DB = Xd_padded_DB((N_sig/2+1):N_sig);
figure(2); subplot(2,3,3);
% plot 2.3
plot(f_star,Xd_padded_DB);
xlabel('f [Hz]');
ylabel('Xd [dB]'); 
title('row DFT at t=0.55 sec, 128 samples padded, Rec');

% Question 8 %
% Different manipulations on the sig of time - 0.55 sec with Blackman
% window

Xd_Blackman = fftshift((fft(x(((selected_row-1)*N_sig+1):((selected_row)*N_sig)).*blackman(N_sig)')));
Xd_Blackman_DB = mag2db(abs(Xd_Blackman));
Xd_Blackman_DB = Xd_Blackman_DB((N_sig/2+1):N_sig);
% plot 2.4
figure(2); subplot(2,3,4);
plot(f_star,Xd_Blackman_DB);
xlabel('f [Hz]');
ylabel('Xd [dB]'); 
title('row DFT at t=0.55 sec, row of matrix, Blackman');

x_samples_blackman = x_samples.*(blackman(N_sig/2)');
Xd_selected_row_Blackman = fftshift(fft(x_samples_blackman));
Xd_selected_row_Blackman_DB = mag2db(abs(Xd_selected_row_Blackman));
Xd_selected_row_Blackman_DB = Xd_selected_row_Blackman_DB((N_sig/4+1):N_sig/2);
% plot 2.5
figure(2); subplot(2,3,5);
plot(f_selected_row,Xd_selected_row_Blackman_DB);
xlabel('f [Hz]');
ylabel('Xd [dB]'); 
title('row DFT at t=0.55 sec, 128 samples, Blackman');

x_padded_Blackman = [x_samples_blackman, zeros(1,128)];
Xd_padded_Blackman = fftshift(fft(x_padded_Blackman));
Xd_padded_Blackman_DB = mag2db(abs(Xd_padded_Blackman));
Xd_padded_Blackman_DB = Xd_padded_Blackman_DB((N_sig/2+1):N_sig);
% plot 2.6
figure(2); subplot(2,3,6);
plot(f_star,Xd_padded_Blackman_DB);
xlabel('f [Hz]');
ylabel('Xd [dB]'); 
title('row DFT at t=0.55 sec, 128 samples padded, Blackman');

% Question 12 %
% Zeros and poles map of filter
% First filter
h1 = [1,1];
zeros1 = [1,1];
poles1 = [1,0];
H1 = tf(zeros1, poles1); 
% plot 3.1
figure(3); subplot(1,2,1);
pzmap(H1); 
axis ([-1 1 -1 1]); 
set(findall(gcf,'type','line'),'linewidth', 3, 'markersize', 12); grid on;
title('H1 Poles and Zeros Map');

% Second filter
h2 = [1,0,0,-1];
zeros2 = [1,0,0,-1];
poles2 = [1,0,0,0];
% plot 3.2
H2 = tf(zeros2, poles2);
figure(3); subplot(1,2,2);
pzmap(H2); 
axis ([-1 1 -1 1]); 
set(findall(gcf,'type','line'),'linewidth', 1, 'markersize', 12); grid on;
title('H2 Poles and Zeros Map');

% Question 13 %
% Singal convolution with filters
% First filter
Xf1 = fftshift(fft(conv(h1,x)));
Xf1_DB = mag2db(abs(Xf1));
Xf1_DB = Xf1_DB((N/2)+1:N);
% plot 4.2
figure(4); subplot(3,1,2);
plot(f,Xf1_DB); 
xlabel('f [Hz]');
ylabel('Xf1 [dB]'); 
title('DFT of x filtered with H1');

% Second filter
Xf2 = fftshift(fft(conv(h2,x)));
Xf2_DB = mag2db(abs(Xf2));
Xf2_DB = Xf2_DB((N/2)+1:N);
% plot 4.3
figure(4); subplot(3,1,3);
plot(f,Xf2_DB); 
xlabel('f [Hz]');
ylabel('Xf1 [dB]'); 
title('DFT of x filtered with H2');