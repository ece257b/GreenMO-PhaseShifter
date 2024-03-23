function [preamble]=generate_preamble(prm,phi)
    % 
    % preamble_legacy=zeros(prm.numUsers,length(sts_t_rep)+2.5*length(lts_t));
    % 
    % for i=1:prm.numUsers
    %     preamble_legacy(1,:)
    % end
    % preamble_legacy(1,:) = [sts_t_rep, lts_t(33:64), lts_t, lts_t];
    % preamble_legacy(2,:) = [circshift(sts_t_rep, [0, prm.TX_SPATIAL_STREAM_SHIFT]), zeros(1, 160)];

    preamble_mimo= zeros(prm.numTxAntenna,prm.numTxAntenna*1.5*length(prm.lts_t));
    for i=1:prm.numTxAntenna
        preamble_mimo(i,(i-1)*1.5*length(prm.lts_t)+1:i*1.5*length(prm.lts_t))=[prm.lts_t(33:64), prm.lts_t];
    end

    if(prm.BF=="SB" || prm.BF=="PB")
        preamble_mimo=pinv(phi)*preamble_mimo;
    end
    preamble=[preamble_mimo];
end