f1=10;f2=100;%the frequencies of sines signal that needs filtered
fs=256000;%sampling frequency
m=(0.3*f1)/(fs/2);%define tansition bandwidth
M=round(8/m);%define the window length
N=M-1;%define the order of filter
b=fir1(N,0.5*f2/(fs/2));%use the firl function to design a filter
%Input parameters are respectively the order number and the cutoff
%frequency of filter
figure(1)
[h,f]=freqz(b,1,512);%amplitude-frequency characteristic graph
plot(f*fs/(2*pi),20*log10(abs(h)))%parameters are respectivelyfrequency and amplitude
xlabel('frequency/Hz');ylabel('gain/dB');title('The gain response oflowpass filter');
figure(2)
subplot(211)
d = daqread('corona_20201124T101707');
t=linspace(0,length(d)/fs,length(d));
plot(t,d)
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram beforefiltering');
axis([0 0.1 -2 2]);
subplot(212)
Fs=fft(d,512);%transform the signal to frequency domain
AFs=abs(d);%take the amplitude
f=(0:255)*fs/512;%frequency sampling
plot(f,AFs(1:256));%plot the frequency domain diagram before filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domaindiagram before filtering');
figure(3)
sf=filter(b,1,d);%use filter function to filter
subplot(211)
plot(t,sf)%plot the signal graph after filtering35
xlabel('time/s');ylabel('amplitude');title('Time-domain diagram afterfiltering');
axis([0.1 0.2 -2 2]);
subplot(212)
Fsf=fft(sf,512);%frequency-domain diagram after filtering
AFsf=abs(Fsf);%the amplitude
f=(0:255)*fs/512;%frequency sampling
plot(f,AFsf(1:256))%plot the frequency domain diagram after filtering
xlabel('frequency/Hz');ylabel('amplitude');title('Frequency-domaindiagram after filtering');
