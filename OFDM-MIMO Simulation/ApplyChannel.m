function [RxAntStreams]=ApplyChannel(prm,channel_t,TxAntStreams)
    RxAntStreams= zeros(prm.numUsers,length(TxAntStreams)+prm.N_SC-1);
    for i=1:prm.numUsers
        for j=1:prm.numTxAntenna
            RxAntStreams(i,:)= RxAntStreams(i,:)+ conv(TxAntStreams(j,:),squeeze(channel_t(i,j,:)));
        end
    end
end