%   

warning off
close all
clear
%% 0-1 Switching

num_ant=4:34;
SINR_per_num_ant=zeros(1,length(num_ant));

num_sc=64;
num_iter=1e2;
rndmUsers=100;
SINR=zeros(rndmUsers,num_sc,num_iter);

for num_Tx_ant=1:length(SINR_per_num_ant)
    %H=load("64sc_4ue_channel_4.2GHz/channel_"+num2str(num_ant(num_Tx_ant))+"x4_64sc.mat");
    %H=load("64sc_4ue_channel_28GHz/channel_"+num2str(num_ant(num_Tx_ant))+"x4_64sc.mat");
    for usrIdx=1:rndmUsers
        H=load("64sc_100rndm_4ue_6GHz/"+num2str(num_ant(num_Tx_ant))+"x4/"+num2str(usrIdx)+".mat");

        H=H.channel;
        channel_noise_metric=0.001;
        
        for i=1:num_iter
            B_Switching= generate_full_rank_binary_matrix(num_ant(num_Tx_ant));
            for j=1:num_sc
                channel_noise=complex(channel_noise_metric.*randn([4,num_ant(num_Tx_ant)]),channel_noise_metric.*randn([4,num_ant(num_Tx_ant)]));
                H_hat=H(:,:,j)+channel_noise;
    
                
                % Rx mat
                R= (H(:,:,j)*B_Switching)/ (H_hat*B_Switching);
                
                SNR=20;
                %S=diag(abs(R)).^2;
                S=sum(abs(R).^2,2);
                noise_pow=S .* 10^(-SNR/10);
                
                %AWGN_noise=complex(sqrt(noise_pow/2).*randn([4,1]),sqrt(noise_pow/2).*randn([4,1]));
                AWGN_noise=complex(sqrt(noise_pow/2).*(2*randi([0,1],4,1)-1),sqrt(noise_pow/2).*(2*randi([0,1],4,1)-1));
            
                % SINR calc 
                SINR(usrIdx,j,i) = SINR_func(R,AWGN_noise);
            end
        end
        
        SINR_per_num_ant(num_Tx_ant)=10*log10(sum(sum(sum(SINR)))/num_sc/num_iter/rndmUsers);
    end
end

figure
plot(SINR_per_num_ant,num_ant)


%% Phase Shifters


num_ant=4:34;
SINR_per_num_ant=zeros(1,length(num_ant));

num_sc=64;
num_iter=1e2;
rndmUsers=100;
SINR=zeros(rndmUsers,num_sc,num_iter);

for num_Tx_ant=1:length(SINR_per_num_ant)
    %H=load("64sc_4ue_channel_4.2GHZ/channel_"+num2str(num_ant(num_Tx_ant))+"x4_64sc.mat");
    %H=load("64sc_4ue_channel_28GHz/channel_"+num2str(num_ant(num_Tx_ant))+"x4_64sc.mat");
    
    for usrIdx=1:rndmUsers
        H=load("64sc_100rndm_4ue_6GHz/"+num2str(num_ant(num_Tx_ant))+"x4/"+num2str(usrIdx)+".mat");
    
        H=H.channel;
        channel_noise_metric=0.001;
            
        for i=1:num_iter
            channel_noise=complex(channel_noise_metric.*randn([4,num_ant(num_Tx_ant)]),channel_noise_metric.*randn([4,num_ant(num_Tx_ant)]));
            phi=exp(-1i*(angle(H(:,:,32)+channel_noise))).';
            % phi=(H(:,:,32)+channel_noise)';
            for j=1:num_sc
                H_hat=H(:,:,j)+channel_noise;
                
                %phi=H_hat';
                % Rx mat
                R= (H(:,:,j)*phi)/ (H_hat*phi);
                
                SNR=20;
                %S=diag(abs(R)).^2;
                S=sum(abs(R).^2,2);
                noise_pow=S .* 10^(-SNR/10);
                
                %AWGN_noise=complex(sqrt(noise_pow/2).*randn([4,1]),sqrt(noise_pow/2).*randn([4,1]));
                AWGN_noise=complex(sqrt(noise_pow/2).*(2*randi([0,1],4,1)-1),sqrt(noise_pow/2).*(2*randi([0,1],4,1)-1));
                
                % SINR calc 
                SINR(usrIdx,j,i) = SINR_func(R,AWGN_noise);
            end
        end
        SINR_per_num_ant(num_Tx_ant)=10*log10(sum(sum(sum(SINR)))/num_sc/num_iter/rndmUsers);
    end
end

hold on
plot(SINR_per_num_ant,num_ant)
xlabel("SINR(dB)")
ylabel("Number of Transmitter Antennas")
legend("RF Switching","Phase Shifter")
%% Functions
function SINR= SINR_func(R, AWGN_noise)
%    [num_row,num_col]=size(R);
    S_values=diag(abs(R)).^2;
    I_values=sum(abs(R).^2,2)-S_values;
    %SINR=mean(10*log10(S_values./(I_values+abs(AWGN_noise).^2)));
    SINR=mean(S_values./(I_values+abs(AWGN_noise).^2));
end

function matrix = generate_full_rank_binary_matrix(m)
    if(m==4)
        matrix=eye(4);
        return
    end
    while true
        matrix = randi([0 1], m, 4);  % Generate random binary matrix
        if rank(matrix) == min(m, 4)  % Check if matrix has full rank
            break;
        end
    end
end