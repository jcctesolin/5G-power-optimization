clc; clear; close all;
L = 4;
v = [1 0 0 0];
G = [0.5; 0.8; 1.0; 0.6];
P = [5; 5; 6; 5];
AP = [1:4]';

disp(G.^2.*P)