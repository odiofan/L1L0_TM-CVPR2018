function y = compress(x,gamma,W)

if nargin<3
    W = 1;
end
    
y = W * ((x./W) .^ (1./gamma));

end