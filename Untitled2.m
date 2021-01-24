d = daqread('corona_20201124T101206');
plot(d)
fs = 256000;
t=linspace(0,length(d)/fs,length(d));
 F = sin(2*pi*50*t);


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

f = 0.6;
n = 6;
a = fir1(n , f, 'high');
b = fir1(n, f, 'low');
o = filter(a , 1, d);
p = filter(b, 1, o);

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
envelope = abs (hilbert(z));
analytic_signal = hilbert(z);
amplitude_envelope = abs(analytic_signal);
instantaneous_phase = unwrap(angle(analytic_signal));
instantaneous_frequency = (diff(instantaneous_phase) /(2.0*pi) * fs);

figure(9)
subplot(3,1,1)
plot(t, real(analytic_signal),t, imag(analytic_signal))
ylabel('Amplitude')
xlabel(' Time(S)')
title('Hilbert Signal')
legend('real','imaginary')
subplot(3,1,2)
plot(t, amplitude_envelope,'r')
ylabel('Amplitude')
xlabel(' Time(S)')
title('Envelope')

 subplot(3,1,3)
% plot(t, fft(analytic_signal))
plot(t,p,'b-',t,imag(hilbert(p)),'k--',t,envelope,'r:')


 xdft2 = p;
 
   xdft2 = xdft2(1:m/2+1);
  psdx1 = (1/(fs*m)) * abs(xdft2).^2;
  
   psdx1(2:end-1) = 2*psdx1(2:end-1);
    freq = 0:fs/m:fs/2;
    y = 10*log10(psdx1);
    figure(10)
  plot(freq/1000,y)
xlabel('Frequency')
ylabel('Intensity')
title('FFT of Hilbert Transform Signal')

figure(18)
[pxx,f]= pwelch(d, window,novarlap, nfft,fs);    % W/Hz power spectral density
PdB_Hz= 10*log10(pxx); 
plot(f/1000,PdB_Hz)  
xlabel('Frequency (khz)')
ylabel('PSD (dB/Hz)')

