function [D1,D2,B2] = Layer_decomp(img,lambda1,lambda2,lambda3)

%% first scale decomposition
[hei,wid,~] = size(img);
[D1,B1] = L1L0(img,lambda1,lambda2);


%% second scale decomposition
scale = 4;
B1_d = imresize(B1,[round(hei/scale),round(wid/scale)],'bilinear');
[~,B2_d] = L1(B1_d,lambda3);
B2_r = imresize(B2_d,[hei,wid],'bilinear');
B2 = bilateralFilter(B2_r,nor(B1),0,1,min(wid,hei)/100,0.05);


D2 = B1 - B2; 

end