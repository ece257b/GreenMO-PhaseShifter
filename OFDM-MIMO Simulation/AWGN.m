function RxAntStreamsNoisy=AWGN(prm, RxAntStreams)
    RxAntStreamsNoisy=zeros(size(RxAntStreams));
    for i=1:prm.numUsers
        noise_pow=sum(abs(RxAntStreams(1,:).^2))/length(RxAntStreams)*10^(-prm.snr/10);
        UserNoise=complex(sqrt(noise_pow/2)*randn(1,length(RxAntStreams)),sqrt(noise_pow/2)*randn(1,length(RxAntStreams)));
        RxAntStreamsNoisy(i,:)=RxAntStreams(i,:)+UserNoise;
    end
end