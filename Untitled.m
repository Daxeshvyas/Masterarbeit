clc
clear all
close all

d = daqread('corona_20201124T101707');

plot(d)
fs = 256000;
t=linspace(0,length(d)/fs,length(d));

 figure(1)
plot(t,d)
 xlabel('Time (S)')
 ylabel('Amplitude')
title('Original Signal')

 

m = length(d);       % original sample length
y = fft(d,m);
 xdft1 = y;
 
   xdft1 = xdft1(1:m/2+1);
  psdx = (1/(fs*m)) * abs(xdft1).^2;
  
   psdx(2:end-1) = 2*psdx(2:end-1);
    freq18 = 0:fs/m:fs/2;
    y18 = 10*log10(psdx);
    figure(2)
  plot(freq18/1000,y18)
  grid on
 
xlabel('Frequency')
ylabel('Intensity')
title('FFT of Signal')

wl = 512;
window = hamming(wl);
 novarlap = wl / 2;
nfft = 2^nextpow2(wl);
figure(11)
spectrogram(d, window, novarlap , 20,fs, 'yaxis'); 
yt = get(gca, 'YTick');
set(gca, 'YTick',yt, 'YTickLabel',yt/1E+3)
ylabel('Frequency (khz)')
xlabel(' Time(S)')
title('Spectrogram')
colorbar

% figure(15)
% bandpass(d,[100 200],fs)
% n = 50;  
% w = 1200/(fs/2);  
% B = fir1(n,w,'high');  
% freqz(B,1,128,8000);  
% figure(15)  
% [h,w]=freqz(B,1,128,8000);  
% plot(w,abs(h)); % Normalized Magnitude Plot  

f = 0.78125;
n = 50;
a = fir1(n , f, 'low');
% b = fir1(n, f, 'low');
% o = filter(a , 1, d);
p = filter(a, 1, d);
% fvtool(1, 1);
figure(8)
subplot(2,1,1);
plot(t,d)
ylabel('Amplitude')
xlabel(' Time(S)')
title('Original ohne Filter')

subplot(2,1,2);
plot(t,p)
ylabel('Amplitude')
xlabel(' Time(S)')
title('Filtered Signal')



z = hilbert(p);
analytic_signal = hilbert(z);
amplitude_envelope = abs(analytic_signal);
instantaneous_phase = unwrap(angle(analytic_signal));
instantaneous_frequency = (diff(instantaneous_phase) /(2.0*pi) * fs);

figure(9)
subplot(2,1,1)
plot(t, real(analytic_signal),t, imag(analytic_signal))
ylabel('Amplitude')
xlabel(' Time(S)')
title('Hilbert Signal')
legend('real','imaginary')
subplot(2,1,2)
plot(t, amplitude_envelope,'r')
ylabel('Amplitude')
xlabel(' Time(S)')
title('Envelope')

% subplot(3,1,3)
% plot(t, fft(analytic_signal))



 xdft2 = p;
 
   xdft2 = xdft2(1:m/2+1);
  psdx1 = (1/(fs*m)) * abs(xdft2).^2;
  
   psdx1(2:end-1) = 2*psdx1(2:end-1);
    freq = 0:fs/m:fs/2;
    y = 10*log10(psdx1);
    figure(10)
  plot(freq/1000,y)
  grid on
xlabel('Frequency')
ylabel('Intensity')
title('FFT of Hilbert Transform Signal')




Fn = fs/2;                                             % Nyquist Frequency (Hz)
E = normalizefreq(d);
Wp = [1   20]/Fn;                                         % Passband Frequency (Normalised)
Ws = [0.5   21]/Fn;                                         % Stopband Frequency (Normalised)
Rp =   1;                                                   % Passband Ripple (dB)
Rs = 150;                                                   % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                             % Filter Order
[Z,P,k] = cheby2(n,Rs,Ws);                                  % Filter Design
[sosbp,gbp] = zp2sos(Z,P,k);                                % Convert To Second-Order-Section For Stability
figure(3)
freqz(sosbp, 2^16, fs)                                      % Filter Bode Plot
filtered_signal = filtfilt(sosbp, gbp, E);    % Filter Signal