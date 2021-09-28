function BER_T = Nakagami_TH(m,MOD,EbN0dB)
EbN0 = 10.^(EbN0dB/10);
for l = 1 : length(MOD)
    M = MOD(l);
    A = (4/log2(M) * (1 - 1/sqrt(M)));
    for i = 1 : length(EbN0dB)
        summ = 0;
        for h = 1 : sqrt(M)/2
            B = ((3 * ((2*h-1)^2)* log2(M))/(M - 1));
            C = m/(EbN0(i) * B);
            D = ((C^m)/12)*(C+0.5)^-m;
            E = ((C^m)/4)*(C + (2/3))^-m;
            F = D + E;
            summ = summ + F;
        end
        BER_T(l,i) = A * summ;
    end
end