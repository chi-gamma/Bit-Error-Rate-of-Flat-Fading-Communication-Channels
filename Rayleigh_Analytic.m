function BER_T = Rayleigh_TH(MOD,EbN0dB)
EbN0 = 10.^(EbN0dB/10);
for l = 1 : length(MOD)
    M = MOD(l);
    A = (2/log2(M) * (1 - 1/sqrt(M)));
    for i = 1 : length(EbN0dB)
        summ = 0;
        for x = 1 : sqrt(M)/2
            B = ((3 * ((2*x-1)^2)* log2(M))/(M - 1)) * EbN0(i);
            C = 1 - sqrt(B/(2+B));
            summ = summ + C;
        end
        BER_T(l,i) = A * summ;
    end
end