%% This function performs the hybrid L1-L0 decomposition

function [D,B]= L1L0(S,lambda1,lambda2)

iter = 15;
[hei,wid] = size(S);

fx = [1, -1];
fy = [1; -1];
otfFx = psf2otf(fx,[hei,wid]);   
otfFy = psf2otf(fy,[hei,wid]);
DxDy = abs(otfFy).^2 + abs(otfFx).^2;

B = S;
C = zeros(hei,wid*2);
E = zeros(hei,wid*2);
L1 = zeros(hei,wid*2);
L2 = zeros(hei,wid*2);
ro1 = 1;
ro2 = 1;
DiffS = [-imfilter(S,fx,'circular'),-imfilter(S,fy,'circular')];
for i = 1:iter
    
    CL = C + L1./ro1;
    EL = DiffS - E - L2./ro2;
    %% for B
    C1L1 = CL(:,1:wid);
    C2L2 = CL(:,1+wid:end);
    E1L3 = EL(:,1:wid);
    E2L4 = EL(:,1+wid:end);
    
    Nomi = fft2(S) + ro1.*conj(otfFx).*fft2(C1L1) + ro1.*conj(otfFy).*fft2(C2L2) ...
        + ro2.*conj(otfFx).*fft2(E1L3) + ro2.*conj(otfFy).*fft2(E2L4);
    Denomi = 1 + (ro1 + ro2) .* DxDy;
    B_new = real(ifft2(Nomi./Denomi));
    DiffB = [-imfilter(B_new,fx,'circular'),-imfilter(B_new,fy,'circular')];
    
    %% for C1, C2, shrinkage
    BL = DiffB - L1./ro1;
    C_new = sign(BL) .* max(abs(BL) - lambda1./ro1 ,0);
    
    %% for E1, E2
    BL = DiffS - DiffB - L2./ro2;
    E_new = BL;
    temp = BL.^2;
    t = temp < 2.*lambda2./ro2;
    E_new(t) = 0;
    
    %% for Li,i=1,2,3,4
    L1_new = L1 + ro1 * (C_new - DiffB);
    L2_new = L2 + ro2 * (E_new - DiffS + DiffB);
    
    %% for ro
    ro1 = ro1 *4;
    ro2 = ro2 *4;
    
    %% update
    B = B_new;
    C = C_new;
    E = E_new;
    L1 = L1_new;
    L2 = L2_new;

end

D = S - B;


end