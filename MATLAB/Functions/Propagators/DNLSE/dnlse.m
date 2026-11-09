function dadT = dnlse(T,a,H)
    I = abs(a).^2;
    NL = I.*a;
    dadT = -1i*(H*a + NL);
end
