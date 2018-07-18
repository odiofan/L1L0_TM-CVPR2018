clear;clc;
%% Parameters
lambda1 = 0.3;  
lambda2 = lambda1*0.01;
lambda3 = 0.1;

%% read files
hdr = hdrimread('1.hdr');
[hei,wid,channel] = size(hdr);

tic;
%% transformation
hdr_h = rgb2hsv(hdr);
hdr_l = hdr_h(:,:,3);
hdr_l = log(hdr_l+0.0001);
hdr_l = nor(hdr_l);

%%  decomposition
[D1,D2,B2] = Layer_decomp(hdr_l,lambda1,lambda2,lambda3);

%% Scaling
sigma_D1 = max(D1(:));
D1s = R_func(D1,0,sigma_D1,0.8,1);
% sigma_D2 = max(D2(:));
% D2s = R_func(D2,0,sigma_D2,0.9,1);
B2_n= compress(B2,2.2,1);
hdr_lnn = 0.8*B2_n + D2 + 1.2*D1s;

%% postprocessing
hdr_lnn = nor(clampp(hdr_lnn,0.005,0.995));
out_rgb = hsv2rgb((cat(3,hdr_h(:,:,1),hdr_h(:,:,2)*0.6,hdr_lnn)));
toc;

figure,imshow(out_rgb)




