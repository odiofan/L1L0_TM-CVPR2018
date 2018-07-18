%% This function performs the TV like decomposition

function [D,B] = L1(S,lambda)

iter = 10;
[hei,wid] = size(S);

fx = [1, -1];
fy = [1; -1];
otfFx = psf2otf(fx,[hei,wid]);
otfFy = psf2otf(fy,[hei,wid]);
DxDy = abs(otfFy).^2 + abs(otfFx).^2;
B = S;
C = zeros(hei,wid);
D = zeros(hei,wid);
T1 = zeros(hei,wid);
T2 = zeros(hei,wid);
ro = 1;
for i = 1:iter
    
    %% for B
    Nomi = fft2(S) + ro * conj(otfFx) .* fft2(C + T1./ro) + ro * conj(otfFy) .* fft2(D + T2./ro);
    Denomi = 1 + ro * DxDy;
    B_new = real(ifft2(Nomi./Denomi));
    
    GradxB = -imfilter(B_new,fx,'circular');
    GradyB = -imfilter(B_new,fy,'circular'); 
    %% for C, D
    BB1 = GradxB - T1./ro;
    BB2 = GradyB - T2./ro;
    C_new = sign(BB1) .* max(abs(BB1) - lambda.*1./ro ,0);
    D_new = sign(BB2) .* max(abs(BB2) - lambda.*1./ro ,0);
    
    %% for T1, T2
    T1_new = T1 + ro * (C_new - GradxB);    
    T2_new = T2 + ro * (D_new - GradyB);
    
    %% for ro
    ro = ro * 2;
    
    %% update
    B = B_new;
    C = C_new;
    D = D_new;
    T1 = T1_new;
    T2 = T2_new;
end

D = S - B;

end