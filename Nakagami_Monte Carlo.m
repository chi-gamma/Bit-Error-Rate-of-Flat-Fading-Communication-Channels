clc;
clear all;
N = 3*(2^19);
data = randi([0 1],1,N); %Generating a uniformly distributed random 1s and 0s
M = [0.5 1 2 3 4 5 6 7 8 9 10 20]; % The m parameters for Nakagami simulation
EbN0dB = 0 : 30;
EbN0 = 10.^(EbN0dB/10);
sigm = 1;
for j = 1 : length(M)
    m = M(j);
    MOD = [4 16 64];
    for k = 1 : length(MOD)
        Mod = MOD(k);
        ns = log2(Mod); % Number of bits per symbol for M-QAM
        b = length(data)/ns;
        serial_data = reshape(data,ns,b); % Serial to Parallel Conversion for QAM
        mod = modem.qammod('M', Mod, 'PHASEOFFSET', 0, 'SYMBOLORDER','BINARY','INPUTTYPE', 'BIT'); % Create modulation object
        mod_data = modulate(mod,serial_data);
        for i = 1 : length(EbN0dB)
            ph_angle = unifrnd(-pi,pi,[1,length(mod_data)]);
            gain = gamrnd(m, sigm/m, 1, length(mod_data));
            h = sqrt(gain).*exp(1i*ph_angle);            
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
            %%  Theoretical BER Calculation
            BER_T = Nakagami_TH(m,MOD,EbN0dB);            
        end
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
    title(['The BER curve for the Nakagami-' num2str(M(j)) ' Fading Channel'])
    axis([0 30 10^-5 1])
end