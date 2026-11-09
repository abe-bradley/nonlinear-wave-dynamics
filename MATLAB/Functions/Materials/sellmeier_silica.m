function [beta0,beta1,beta2] = sellmeier_silica(lambda0)
    % Fused silica (Malitson Sellmeier). lam [um]; outputs: beta0 [rad/um], beta1 [ns/um], beta2 [ns^2/um]
    
    c = 299792458; % um/us
    c = 1e-3*c; % um/ns
    
    B = [0.6961663 0.4079426 0.8974794];
    C = [0.0684043 0.1162414 9.896161].^2; % um^2
    
    sz = size(lambda0); lambda0 = lambda0(:);
    l2  = lambda0.^2;   den = l2 - C;
    
    S   = 1 + sum(B.*l2./den,2);                                       % n^2
    Sp  = sum((-2*B.*C).*lambda0./den.^2,2);                               % d(n^2)/dlam
    Spp = sum((-2*B.*C)./den.^2 + (8*B.*C).*l2./den.^3,2);             % d^2(n^2)/dlam^2
    
    n   = sqrt(S);
    n1  = Sp./(2*n);
    n2  = Spp./(2*n) - (Sp.^2)./(4*n.^3);
    
    beta0   = 2*pi*n./lambda0;
    beta_l  = 2*pi*(n1./lambda0 - n./lambda0.^2);
    beta_ll = 2*pi*(n2./lambda0 - 2*n1./lambda0.^2 + 2*n./lambda0.^3);
    
    A   = -(lambda0.^2)./(2*pi*c);                                         % d/dw = A d/dlam,  w=2*pi*c/lam
    A_l = -(lambda0)./(pi*c);
    
    beta1 = A.*beta_l;
    beta2 = (A.*A_l).*beta_l + (A.^2).*beta_ll;
    
    beta0 = reshape(beta0,sz); beta1 = reshape(beta1,sz); beta2 = reshape(beta2,sz);
end
