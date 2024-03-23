warning off
close all
clear

%% 0-1 Switching


num_ant=4:34;
rndmUsers=5;
SNR=5;
avgBER=zeros(rndmUsers,length(num_ant));
BF="SB";

for num_Tx_ant=1:length(num_ant)
    num_ant(num_Tx_ant)
    for usrIdx=1:rndmUsers
        H=load("64sc_50rndm_4ue_6GHz/"+num2str(num_ant(num_Tx_ant))+"x4/"+num2str(usrIdx)+".mat");
        H=H.channel;
        for snrIdx=1:length(SNR)
            BER= OFDM_MIMO(H,BF);
            avgBER(usrIdx,num_Tx_ant)=mean(BER);
        end
    end
end
    

figure
semilogy(num_ant,mean(avgBER))

%% Phase Shifters

num_ant=4:34;
rndmUsers=5;
SNR=5;
avgBER=zeros(rndmUsers,length(num_ant));
BF="PB";

for num_Tx_ant=1:length(num_ant)
    num_ant(num_Tx_ant)
    for usrIdx=1:rndmUsers
        H=load("64sc_50rndm_4ue_6GHz/"+num2str(num_ant(num_Tx_ant))+"x4/"+num2str(usrIdx)+".mat");
        H=H.channel;
        for snrIdx=1:length(SNR)
            BER= OFDM_MIMO(H,BF);
            avgBER(usrIdx,num_Tx_ant)=mean(BER);
        end
    end
end

hold on
semilogy(num_ant,mean(avgBER))
xlabel("Number of Basestation Antennas")
ylabel("Average BER")
legend("RF Switching","Phase Shifter")