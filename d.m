clc
clear all
close all
daqinfo = daqread('corona_20201124T101707', 'info');
daqinfo.ObjInfo.Channel