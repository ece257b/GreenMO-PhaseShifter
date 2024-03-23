function [tx_data, tx_payload_vec, tx_syms]=generate_ofdm_data(prm,channel_f,phi)
    

    ifft_in_mat=zeros(prm.numUsers,prm.N_SC,prm.N_OFDM_SYMS);
    % Repeat the pilots across all OFDM symbols
    pilots_mat = repmat(prm.pilots, 1, prm.N_OFDM_SYMS);
    for i=1:prm.numUsers
        tx_bits = randi(2, 1, prm.number_of_bits) - 1;
        
        % Forward Error Correction
        tx_data(i,:) = double([tx_bits zeros(1,prm.trellis_end_length) ]);    % 8 bits padding
        trel = poly2trellis(7, [171 133]);              % Define trellis
        tx_code = convenc(tx_data(i,:),trel);            % convultional encoder
        
        % bits to signal space mapping these are you are x_k from the class
        tx_syms(i,:) = mapping(tx_code', prm.MOD_ORDER, 1);
        
        % figure(1);
        % scatter(real(tx_syms), imag(tx_syms),'filled');
        % title(' Signal Space of transmitted bits');
        % xlabel('I'); ylabel('Q');
        
        % Reshape the symbol vector to a matrix with one column per OFDM symbol,
        tx_syms_mat = reshape(tx_syms(i,:), length(prm.SC_IND_DATA), prm.N_OFDM_SYMS);
        
        


        % Insert the data and pilot values; other subcarriers will remain at 0
        ifft_in_mat(i,prm.SC_IND_DATA, :)   = tx_syms_mat;
        ifft_in_mat(i,prm.SC_IND_PILOTS, :) = pilots_mat;

    end

    % Perform precoding
    % 
    % ifft_in_mat_ant=zeros(prm.numTxAntenna,prm.N_SC,prm.N_OFDM_SYMS);
    % for i=1:prm.N_OFDM_SYMS
    %     for j=1:prm.N_SC
    %         ifft_in_mat_ant(:,j,i)=channel_f(:,:,j)\ifft_in_mat(:,j,i);
    %         %ifft_in_mat(:,j,i)=inv(channel_f(:,:,j)*phi)*ifft_in_mat(:,j,i);
    %     end
    % end
    % tx_payload_vec=zeros(prm.numTxAntenna,(prm.N_SC+prm.CP_LEN)*prm.N_OFDM_SYMS);
    % for i=1:prm.numTxAntenna
    %     % Perform the IFFT --> frequency to time translation
    %     tx_payload = ifft(squeeze(ifft_in_mat_ant(i,:,:)), prm.N_SC, 1);
    %     % Insert the cyclic prefix
    %     if(prm.CP_LEN > 0)
    %         tx_cp = tx_payload((end-prm.CP_LEN+1 : end), :);
    %         % Reshape to a vector
    %         tx_payload_vec(i,:) = reshape([tx_cp; tx_payload], 1, numel([tx_cp; tx_payload]));
    %     end
    % 
    % end  

    switch prm.BF
        case("DB")
            ifft_in_mat_ant=zeros(prm.numTxAntenna,prm.N_SC,prm.N_OFDM_SYMS);
            for i=1:prm.N_OFDM_SYMS
                for j=[prm.SC_IND_DATA,prm.SC_IND_PILOTS]
                    ifft_in_mat_ant(:,j,i)=channel_f(:,:,j)\ifft_in_mat(:,j,i);
                    %ifft_in_mat(:,j,i)=inv(channel_f(:,:,j)*phi)*ifft_in_mat(:,j,i);
                end
            end
            tx_payload_vec=zeros(prm.numTxAntenna,(prm.N_SC+prm.CP_LEN)*prm.N_OFDM_SYMS);
            for i=1:prm.numTxAntenna
                % Perform the IFFT --> frequency to time translation
                tx_payload = ifft(squeeze(ifft_in_mat_ant(i,:,:)), prm.N_SC, 1);
                % Insert the cyclic prefix
                if(prm.CP_LEN > 0)
                    tx_cp = tx_payload((end-prm.CP_LEN+1 : end), :);
                    % Reshape to a vector
                    tx_payload_vec(i,:) = reshape([tx_cp; tx_payload], 1, numel([tx_cp; tx_payload]));
                end

            end  

         case("SB")
            ifft_in_mat_ant=zeros(prm.numUsers,prm.N_SC,prm.N_OFDM_SYMS);
            for i=1:prm.N_OFDM_SYMS
                for j=1:prm.N_SC
                    ifft_in_mat_ant(:,j,i)=(channel_f(:,:,j)*phi)\ifft_in_mat(:,j,i);
                    %ifft_in_mat(:,j,i)=inv(channel_f(:,:,j)*phi)*ifft_in_mat(:,j,i);
                end
            end
            tx_payload_vec=zeros(prm.numUsers,(prm.N_SC+prm.CP_LEN)*prm.N_OFDM_SYMS);
            for i=1:prm.numUsers
                % Perform the IFFT --> frequency to time translation
                tx_payload = ifft(squeeze(ifft_in_mat_ant(i,:,:)), prm.N_SC, 1);
                % Insert the cyclic prefix
                if(prm.CP_LEN > 0)
                    tx_cp = tx_payload((end-prm.CP_LEN+1 : end), :);
                    % Reshape to a vector
                    tx_payload_vec(i,:) = reshape([tx_cp; tx_payload], 1, numel([tx_cp; tx_payload]));
                end

            end  
        case("PB")
            ifft_in_mat_ant=zeros(prm.numUsers,prm.N_SC,prm.N_OFDM_SYMS);
            for i=1:prm.N_OFDM_SYMS
                for j=1:prm.N_SC
                    ifft_in_mat_ant(:,j,i)=(channel_f(:,:,j)*phi)\ifft_in_mat(:,j,i);
                    %ifft_in_mat(:,j,i)=inv(channel_f(:,:,j)*phi)*ifft_in_mat(:,j,i);
                end
            end
            tx_payload_vec=zeros(prm.numUsers,(prm.N_SC+prm.CP_LEN)*prm.N_OFDM_SYMS);
            for i=1:prm.numUsers
                % Perform the IFFT --> frequency to time translation
                tx_payload = ifft(squeeze(ifft_in_mat_ant(i,:,:)), prm.N_SC, 1);
                % Insert the cyclic prefix
                if(prm.CP_LEN > 0)
                    tx_cp = tx_payload((end-prm.CP_LEN+1 : end), :);
                    % Reshape to a vector
                    tx_payload_vec(i,:) = reshape([tx_cp; tx_payload], 1, numel([tx_cp; tx_payload]));
                end

            end  
    end
end