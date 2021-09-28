clc;
clear all;
N = 3*(2^19);
data = randi([0 1],1,N); %Generating a uniformly distributed random 1s and 0s
totPower = 1;
K = 15;%[0 1 2 3 4 5 6 7 8 9 10 20] ; % The k-factors for the Rician simulation
EbN0dB = 0 : 30;
EbN0 = 10.^(EbN0dB/10);
sigm = 1;
for j = 1 : length(K)
    k_f = K(j);
    s = sqrt(k_f/(k_f+1)*totPower); %Non-Centrality Parameter
    sigma = totPower/sqrt(2*(k_f+1));
    MOD = [4 16 64];
    for k = 1 : length(MOD)
        Mod = MOD(k);
        ns = log2(Mod); % Number of bits per symbol for M-QAM
        b = length(data)/ns;
        serial_data = reshape(data,ns,b); % Serial to Parallel Conversion for QAM
        mod = modem.qammod('M', Mod, 'PHASEOFFSET', 0, 'SYMBOLORDER','BINARY','INPUTTYPE', 'BIT'); % Create modulation object
        mod_data = modulate(mod,serial_data);
        for i = 1 : length(EbN0dB)
            h = ((sigma*randn(1,length(mod_data))+s)+1i*(randn(1,length(mod_data))*sigma+0)); %Rician Fading - single tap            
            %%  AWGN CHANNEL
            SNR = EbN0dB + 10*log10(ns);
            a = h.*mod_data;
            y = awgn(a,SNR(i),'measured');
            y_rec = y./h; % Channel Equalization || Fading Phase Compensation
            %%  QAM Demodulation
            data1 = y_rec(:);
            qam_data = data1.';
            demod = modem.qamdemod('M', Mod, 'PHASEOFFSET', 0, 'SYMBOLORDER', 'BINARY', 'OUTPUTTYPE', 'BIT'); % Create demodulation object
            demod_data = demodulate(demod,qam_data); % Modulating data to QAM
            %%  Symbols to Msg Conversion
            data1 = demod_data(:);
            received_msg = data1.';
            received_msg = received_msg(1,1:length(data));
            [ErrBits BER(k,i)] = biterr(received_msg,data);  % Calculating BER by comparing the transmitted and received
                        
        end
        BER_T = Rice_TH(k_f,MOD,EbN0dB);
    end
   
    figure(j)
    semilogy(EbN0dB,BER(1,:),'r--','Linewidth',1);
    hold on
    semilogy(EbN0dB,BER_T(1,:),'kd','Linewidth',1);
    hold on
    semilogy(EbN0dB,BER(2,:),'b--','Linewidth',1);
    hold on
    semilogy(EbN0dB,BER_T(2,:),'kv','Linewidth',1);
    hold on
    semilogy(EbN0dB,BER(3,:),'m--','Linewidth',1);
    hold on
    semilogy(EbN0dB,BER_T(3,:),'ko','Linewidth',1);
    grid on
    xlabel('E_{b}/N_{0} (dB)');
    ylabel('Bit Error Rate');
    legend('4-QAM(MC)','4-QAM(TH)','16-QAM(MC)','16-QAM(TH)','64-QAM(MC)','64-QAM(TH)');
    title(['The BER curve for the Rician (k = ' num2str(K(j)) ') Fading Channel'])
    axis([0 30 10^-5 1])
end