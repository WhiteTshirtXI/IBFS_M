clear all, close all, clc



cmap = zeros(64,3);
cmap(32,:) = [1 1 1];
cmap(1:31,1) = linspace(0,1,31);
cmap(1:31,2) = linspace(0,1,31);
cmap(1:31,3) = 1;
cmap(33:end,1) = 1;
cmap(33:end,2) = linspace(1,0,32);
cmap(33:end,3) = linspace(1,0,32);


save('cmap.mat','cmap');