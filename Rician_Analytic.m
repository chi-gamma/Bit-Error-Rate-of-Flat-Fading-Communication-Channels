function BER_T = Rice_TH(k,MOD,EbN0dB)
EbN0 = 10.^(EbN0dB/10);
for l = 1 : length(MOD)
    M = MOD(l);
    A = (4/log2(M) * (1 - 1/sqrt(M)));
    B = exp(-1*k);
    for i = 1 : length(EbN0dB)
        C = (k + 1)/EbN0(i);
        summ = 0;
        for j = 1 : sqrt(M)/2;
            b2g = ( ((2*j-1)^2) * (3*log2(M)) )/(M-1);
            syms n summ2
            D = ( (k^n) * (C^(n+1)) )/factorial(n);
            E = D * (2 / ( b2g^(n+1) ));
            F =  1/(24*(( C/b2g+(1/2) )^(n+1))) + 1/(8*(( C/b2g+(2/3) )^(n+1)));
            G = E * F;
            summ2 = symsum(G, n, 0, 100);
            H = B * summ2;
            summ = H + summ;
        end
        BER_T(l,i) = A * summ;
    end
end    