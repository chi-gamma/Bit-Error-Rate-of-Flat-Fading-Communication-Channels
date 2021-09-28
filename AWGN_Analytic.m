function BER_T = AWGN_TH(MOD,EbN0dB)
EbN0 = 10.^(EbN0dB/10);
for l = 1 : length(MOD)
    M = MOD(l);
    A = (4/log2(M) * (1 - 1/sqrt(M)));
    for i = 1 : length(EbN0dB)
        summ = 0;
        for h = 1 : sqrt(M)/2
            B = (2*h - 1);
            C = 0.5*erfc(B*sqrt((1.5*log2(M)*EbN0(i))/(M-1)));
            summ = summ + C;
        end
        BER_T(l,i) = A * summ;
    end
end