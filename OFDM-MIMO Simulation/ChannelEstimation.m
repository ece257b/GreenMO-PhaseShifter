function [channel_f]=ChannelEstimation(prm,RxAntStreams)
    RxPreamble=RxAntStreams(:,1:1.5*length(prm.lts_t)*prm.numTxAntenna);
    channel_f=zeros(prm.numUsers,prm.numTxAntenna,prm.N_SC);
    for i=1:prm.numTxAntenna
        for j=1:prm.numUsers
            Rx_lts_t= RxPreamble(j,(i-1)*1.5*length(prm.lts_t)+1:i*1.5*length(prm.lts_t));
            Rx_lts_t=Rx_lts_t(length(prm.lts_t)/2+1:end);
            channel_f(j,i,:)=fft(Rx_lts_t)./prm.lts_f;
        end
    end
end