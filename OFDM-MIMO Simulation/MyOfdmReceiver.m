
 function [decoded_data]= MyOfdmReceiver(data)
 %% run transmitter code to load sts and lts and other parameters 
 OFDM_TX; 
 
%% Rx processing params

rx_data = data;          % run OFDM tx code to get raw_rx_dec


LTS_CORR_THRESH = 0.8;         % Normalized threshold for LTS correlation
% Usage: Find all peaks whose magnitude is greater than 0.8 times
% the maximum magnitude after cross correlation (Packet Detection)

% Repeat the following code for each packet

%% Packet Detection

% ideas: Cross correlation of received signal with LTS or use STS to detect the packet?

% length_samples= length(rx_data) - 200; %cross corelation based on STS
% sample=16;
% 
% while( sample < length_samples) %cross corelation based on STS
% 
%     output(sample)= rx_data(sample-length(sts_t) + (1:length(sts_t))) * rx_data(sample + (1:length(sts_t)))' ./norm(rx_data(sample+(1:length(sts_t))))^2; 
% 
%     sample= sample+1;
% 
% end
% 
% output= output./max(abs(output));
% 

% cross correlation based on STS or LTS


num_packet=1;
i=1;
E=abs(rx_data).^2;
E=E/max(E);
while(true)
    if(E(i) >= 0.001) % 0.001 for all real world and 
        rx_data=rx_data(i:end);

        know_seq=sts_t;
        len=length(know_seq);
        R=zeros(length(rx_data)-len+1,1); 
        for m=1:length(R)
            R(m)=sum(conj(know_seq).*rx_data(m:m+len-1))/sum(conj(rx_data(m:m+len-1)).*rx_data(m:m+len-1));  
        end
        R=abs(R)/max(abs(R));


        [~, locations] = findpeaks(R,'MinPeakHeight',0);
        indexes=length(locations);

        k=1;
        while(k<=length(locations)) 
            flag=0;
            m=29;
            for j=1:29
                flag=1;
                if(any(locations==locations(k)+16*m) ~= 1)
                    flag=0;
                    break
                end
                m=m-1;
            end
            if(flag==1)
                good_locs=locations(k);
                packet(num_packet,:)=rx_data(1,good_locs:30*length(sts_t)+2.5*length(lts_t)+(N_SC+CP_LEN)*N_OFDM_SYMS+good_locs-1);
                rx_data=rx_data(1,30*length(sts_t)+2.5*length(lts_t)+(N_SC+CP_LEN)*N_OFDM_SYMS+good_locs+1:end);
                E=abs(rx_data).^2;
                E=E/max(E);
                break
            end
            k=k+1;
        end


        if(length(rx_data)<30*length(sts_t)+2.5*length(lts_t)+(N_SC+CP_LEN)*N_OFDM_SYMS)
            break
        else
            num_packet=num_packet+1;
            i=1;
        end

    else
        i=i+1;
        continue
    end
end

% good_locs=101;
% packet(1,:)=rx_data(1,good_locs:30*length(sts_t)+2.5*length(lts_t)+(N_SC+CP_LEN)*N_OFDM_SYMS+good_locs-1);

% 
% know_seq=sts_t;
% len=length(know_seq);
% R=zeros(length(rx_data)-len+1,1);  
% for m=1:length(R)
%     R(m)=sum(conj(know_seq).*rx_data(m:m+len-1))/sum(conj(rx_data(m:m+len-1)).*rx_data(m:m+len-1));  
% end
% R=R/max(abs(R));
% % [~, locations] = findpeaks(abs(R), 'MinPeakHeight', 0.6);
% [~, locations] = findpeaks(abs(R),'MinPeakHeight', 0.2);

% 
% ind=1;
% indexes=length(locations);
% i=1;
% while(i<=length(locations)) 
%     flag=0;
%     m=29;
%     for j=1:29
%         flag=1;
%         if(any(locations==locations(i)+16*m) ~= 1)
%             flag=0;
%             break
%         end
%         m=m-1;
%     end
%     if(flag==1)
%         good_locs(ind)=locations(i);
%         ind=ind+1;
%         i=find(locations==locations(i)+length(sts_t)*29);
%     end
%     i=i+1;
% end


% sliding window based on STS or LTS
% len=length(sts_t);
% R=zeros(length(rx_data)-2*len+1,1);
% for m=1:length(rx_data)-2*len+1
%     R(m)=sum(abs(rx_data(m+len:m+2*len-1)).*abs(rx_data(m+len:m+2*len-1)))./sum(abs(rx_data(m:m+len-1)).*abs(rx_data(m:m+len-1)));  
% end

% [~, locations] = findpeaks(abs(R), 'MinPeakHeight', 0.7);
% start_time=locations(1)+len;

%rx_data=tx_vec;

% Output: Single packet extracted from rx_data
% with knowledge of preamble (LTS) indices and payload vector indices
%% CFO estimation and correction
% Use two copies of LTS for cross-correlation (Reference: Thesis)

%rx_data=rx_data(2:end);
[num_of_packet,packet_size]=size(packet);
for i=1:num_of_packet
    sts=packet(i,1:16*30);
    
    lts1(i,:)=packet(i,16*30+32+1:16*30+32+64);
    lts2(i,:)=packet(i,16*30+32+64+1:16*30+32+2*64);
    
    CFO=phase(lts1(i,:)*lts2(i,:)')/2/pi/N_SC;
    packet(i,:)=packet(i,:) .* exp(1i*2*pi*CFO*[0:packet_size-1]/SAMP_FREQ);
    lts1(i,:)=packet(i,16*30+32+1:16*30+32+64);
    lts2(i,:)=packet(i,16*30+32+64+1:16*30+32+2*64);
end

% delta_f_vec=phase(lts2./lts1)/2/64/pi; % corrosponding CFO 
% CFO=sum(delta_f_vec)/length(delta_f_vec); % average CFO
% rx_data=rx_data .* exp(1i*2*pi*CFO*[0:length(rx_data)-1]);


% res_CFO1=(phase(lts1(end))-phase(lts1(1)))/64;
% res_CFO2=(phase(lts2(end))-phase(lts2(1)))/64;
% res_CFO=(res_CFO1+res_CFO2)/2;
% rx_data=rx_data .* exp(1i*2*pi*res_CFO*[0:length(rx_data)-1]);




% Output: Packet with each value multiplied by CFO correction factor

%% CP Removal
% Refer to the process used to add CP at TX
% Converting vector back to matrix form will help

packet=packet(:,30*length(sts_t)+2.5*length(lts_t)+1:30*length(sts_t)+2.5*length(lts_t)+(N_SC+CP_LEN)*N_OFDM_SYMS);
[num_of_packet,packet_size]=size(packet);
% OFDM_syms=zeros(num_of_packet,packet_size/N_OFDM_SYMS,N_OFDM_SYMS);

for i=1:num_of_packet
    OFDM_syms(i,:,:)=reshape(packet(i,:),[packet_size/N_OFDM_SYMS,N_OFDM_SYMS]);
    OFDM_syms_t(i,:,:)=OFDM_syms(i,CP_LEN+1:end,:);
end
% Output: CP free payload matrix of size (N_SC * N_OFDM_SYMS)

%% FFT
% Refer to IFFT perfomed at TX
for i=1:num_of_packet
    OFDM_syms_f(i,:,:)=fft(squeeze(OFDM_syms_t(i,:,:)),N_SC, 1);
end


% Output: Symbol matrix in frequency domain of same size


%% Channel estimation and correction
% Use the two copies of LTS and find channel estimate (Reference: Thesis)
% Convert channel estimate to matrix form and equlaize the above matrix

for i=1:num_of_packet
    H1=fft(lts1(i,:))./lts_f;
    H2=fft(lts2(i,:))./lts_f;
    H=(H1+H2)/2;
    OFDM_syms_eq(i,:,:)=squeeze(OFDM_syms_f(i,:,:))./repmat(H.',1,N_OFDM_SYMS);
end

% Output : Symbol equalized matrix in frequency domain of same size

%% Advanced topics: 
%% SFO estimation and correction using pilots
% SFO manifests as a frequency-dependent phase whose slope increases
% over time as the Tx and Rx sample streams drift apart from one
% another. To correct for this effect, we calculate this phase slope at
% each OFDM symbol using the pilot tones and use this slope to
% interpolate a phase correction for each data-bearing subcarrier.

for i=1:num_of_packet

    received_pilots = squeeze(OFDM_syms_eq(i,SC_IND_PILOTS,:))./repmat(pilots,1,N_OFDM_SYMS);
    received_pilots_new=zeros(4,N_OFDM_SYMS);
    received_pilots_new(1,:)=received_pilots(3,:);
    received_pilots_new(2,:)=received_pilots(4,:);
    received_pilots_new(3,:)=received_pilots(1,:);
    received_pilots_new(4,:)=received_pilots(2,:);
    
    SFO=zeros(1,N_OFDM_SYMS);
    res_CFO=zeros(1,N_OFDM_SYMS);
    f=-32:31;
    Delta_f(i,:,:)=zeros(N_SC,N_OFDM_SYMS);
    for j=1:N_OFDM_SYMS
        coeff=polyfit([44-64 58-64 8 22],flip(phase(received_pilots_new(:,j))),1);
        SFO(j)=coeff(1);
        res_CFO(j)=coeff(2);
        Delta_f(i,:,j)=coeff(1)*f+coeff(2);
    end
end

% Output: Symbol equalized matrix with pilot phase correction applied

%% Phase Error Correction using pilots
% Extract the pilots and calculate per-symbol phase error

for i=1:num_of_packet
    OFDM_syms_eq(i,1:32,:)=OFDM_syms_eq(i,1:32,:)./exp(1i*Delta_f(i,end-31:end,:));
    OFDM_syms_eq(i,33:64,:)=OFDM_syms_eq(i,33:64,:)./exp(1i*Delta_f(i,1:32,:));
    % 
    rx_syms(i,:,:)=OFDM_syms_eq(i,SC_IND_DATA,:); % data extraction 
    final_rx_syms(i,:)=rx_syms(i,:);
end
% Output: Symbol equalized matrix with pilot phase correction applied
% Remove pilots and flatten the matrix to a vector rx_syms


%% Demodulation


% FEC decoder   
for i=1:num_of_packet
    figure
    scatter(real(final_rx_syms(i,:)), imag(final_rx_syms(i,:)),'filled');
    title(' Signal Space of received bits');
    xlabel('I'); ylabel('Q');
    Demap_out(i,:) = demapper(final_rx_syms(i,:),MOD_ORDER,1); %scale for simulation, and 1 for real world data
    
    % viterbi decoder
    decoded_data(i,:) = vitdec(Demap_out(i,:),trel,7,'trunc','hard');
end
% decoded_data is the final output corresponding to tx_data, which can be used
% to calculate BER