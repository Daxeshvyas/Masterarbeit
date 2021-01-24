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


 xdft2 = p;
 