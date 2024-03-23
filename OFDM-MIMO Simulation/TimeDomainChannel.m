function [channel_t]=TimeDomainChannel(prm,channel_f)
channel_t=zeros(size(channel_f));
    for i=1:prm.numTxAntenna
        for j=1:prm.numUsers
            channel_t(j,i,:)=(ifft(channel_f(j,i,:)));
        end
    end

end