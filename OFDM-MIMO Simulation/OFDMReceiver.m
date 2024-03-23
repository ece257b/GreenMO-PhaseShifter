function decodedBits=OFDMReceiver(prm,RxAntStreams)
    
    decodedBits=zeros(prm.numUsers,prm.number_of_bits+prm.trellis_end_length);
    % for i=1:prm.numUsers
    %     OFDMPacket_tmp= Rx_Ant_Signal(i,i:prm.numUsers:end);
    % 
    %     % removing preamble
    %     %OFDMPacket= OFDMPacket_tmp(prm.numUsers*1.5*64+1:end);
    % 
    %     OFDMSymbols_with_CP= reshape(OFDMPacket, prm.N_SC+ prm.CP_LEN,prm.N_OFDM_SYMS);
    %     OFDMSymbols(i,:,:)=OFDMSymbols_with_CP(prm.CP_LEN+1:end, :);
    % end
    RxAntStreams=RxAntStreams(:,prm.numTxAntenna*1.5*length(prm.lts_t)+1:end);
    for i=1:prm.numUsers
        OFDMPacket_tmp= RxAntStreams(i,1:(prm.CP_LEN+prm.N_SC)*prm.N_OFDM_SYMS);
        
        % removing preamble
        %OFDMPacket= OFDMPacket_tmp(prm.numUsers*1.5*64+1:end);

        OFDMSymbols_with_CP= reshape(OFDMPacket_tmp, prm.N_SC+ prm.CP_LEN,prm.N_OFDM_SYMS);
        OFDMSymbols=OFDMSymbols_with_CP(prm.CP_LEN+1:end, :);
        % FFT
        % Refer to IFFT perfomed at TX
        OFDM_syms_f=fft(OFDMSymbols,prm.N_SC, 1);
        

        rx_syms=OFDM_syms_f(prm.SC_IND_DATA,:); % data extraction 
        final_rx_syms=rx_syms(:);
        Demap_out = demapper(final_rx_syms,prm.MOD_ORDER,1); %scale for simulation, and 1 for real world data
        % viterbi decoder
        decodedBits(i,:) = vitdec(Demap_out,prm.trel,7,'trunc','hard');
        
    end  
    
end
    
