clc;
clear all;
N = 3*2^19;
data = randi([0 1],1,N);
M = [4 16 64]; % M-ary
EbN0dB = 0:30;
EbN0 = 10.^(EbN0dB/10);
for k = 1 : length(M)
     A = (2/log2(M(k)) * (1 - 1/sqrt(M(k)))); % The first factor in the theoretical M-QAM BER
    %%  M-QAM Modulation
    ns = log2(M(k)); % Number of bits per symbol for M-QAM   
    b = length(data)/ns;
    serial_data = reshape(data,ns,b); % Serial to Parallel Conversion for QAM
    mod = modem.qammod('M', M(k), 'PHASEOFFSET', 0, 'SYMBOLORDER','BINARY','INPUTTYPE', 'BIT'); % Create modulation object
    mod_data = modulate(mod,serial_data);
    %%  AWGN CHANNEL
    SNR = EbN0dB + 10*log10(ns);
    for i = 1 : length(EbN0dB)       
        h = 1/sqrt(2)*(randn(1,length(mod_data))+1i*randn(1,length(mod_data))); %Rayleigh Flat Fading factor- single tap
        a = h.*mod_data;        
        received_data = awgn(a,SNR(i),'measured');  % Adding AWGN
        estimated_data = received_data./h; % Channel Equalization
        %%  QAM Demodulation
        data1 = estimated_data(:);
        qam_data = data1.';
        demod = modem.qamdemod('M', M(k), 'PHASEOFFSET', 0, 'SYMBOLORDER', 'BINARY', 'OUTPUTTYPE', 'BIT'); % Create demodulation object
        demod_data = demodulate(demod,qam_data); % Modulating data to QAM
        %%  Symbols to Msg Conversion
        data1 = demod_data(:);
        received_msg = data1.';
        received_msg = received_msg(1,1:length(data));
        [ErrBits BER(i,k)] = biterr(received_msg,data);  % Calculating BER by comparing the transmitted and received message      
        %%  Theoretical BER Calculation
        BER_T = Rayleigh_TH(M,EbN0dB);

    end
end

semilogy(EbN0dB,BER(:,1),'r--','Linewidth',1);
hold on
semilogy(EbN0dB, BER_T(1,:),'kd');
hold on
semilogy(EbN0dB,BER(:,2),'b--');
hold on
semilogy(EbN0dB, BER_T(2,:),'kv');
hold on
semilogy(EbN0dB,BER(:,3),'m--');
hold on
semilogy(EbN0dB, BER_T(3,:),'ko');
grid on
legend('4-QAM(MC)','4-QAM(TH)','16-QAM(MC)','16-QAM(TH)','64-QAM(MC)','64-QAM(TH)');
title('The BER curve for the Rayleigh Fading channel');
ylabel('Bit Error Rate');
xlabel('E_{b}/N_{0} (dB)');
axis([0 30 10^-5 1])