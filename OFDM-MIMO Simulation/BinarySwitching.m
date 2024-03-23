function    TxAntStreams=BinarySwitching(prm,mimo_signals,phi)
    
    % preSwitched=(pinv(phi))*mimo_signals;
    % size_samps=length(preSwitched);

    size_samps=length(mimo_signals);

    N=50;
    w = blackmanharris(N)';
    % figure
    % plot(w);
    % title("BlackmanHarris Window(N=50)")
    % xlabel("Sample Index")
    % ylabel("Amplitude")
    w=w./(sum(w.^2)/N);

    interleaved_t=zeros(1,prm.numUsers*length(mimo_signals));
    for i=1:prm.numUsers
        interleaved_t(i:prm.numUsers:end)=delayseq(mimo_signals(i,:),-(i-1)/prm.numUsers);
        % interleaved_t(i:prm.numUsers:end)=delayseq(mimo_signals(i,:),-(i-1)*N);
    end
    
    DAC_input_t=upsample(interleaved_t,N);
    DAC_output_t=conv(DAC_input_t,w);
    DAC_output_t=DAC_output_t(1:length(DAC_input_t));
    
    Switching=zeros(prm.numTxAntenna,length(DAC_output_t));

    for i=1:prm.numTxAntenna
        Switching(i,:)=repmat(repelem(phi(i,:),N),1,length(DAC_output_t)/prm.numUsers/N);
    end    

    Switch_output_t=Switching.*repmat(DAC_output_t,prm.numTxAntenna,1);


    Switch_output_f=fftshift(fft(Switch_output_t,length(Switch_output_t),2),2);
    % Switch_output_f=fftshift(fft(Switch_output_t(1,:)));

    % Compensating the LPF phase
    phase_cmpnst_switch_output_f= Switch_output_f.*repmat(exp(-1i*phase( fftshift( fft( w,length( Switch_output_t))))),prm.numTxAntenna,1);
    
    % Compensating the delay phase
    % nfft=length(Switch_output_t);
    % fbins = 2 * pi * ((0:1:nfft-1) - floor(nfft/2)) / nfft;
    % 
    % % 
    % phase_cmpnst_switch_output_f(2,:)= phase_cmpnst_switch_output_f(2,:).* exp(-1j * -N * fbins);
    Switch_output_f=phase_cmpnst_switch_output_f;
    
    % LPF to remove side lobes
    % TxAntStreams_f=Switch_output_f(:,length(Switch_output_f)/2-size_samps/2+1:length(Switch_output_f)/2+size_samps/2);
    % TxAntStreams=ifft(ifftshift(TxAntStreams_f,2), size_samps, 2);

    TxAntStreams_f=Switch_output_f(:,length(Switch_output_f)/2-size_samps/2+1:length(Switch_output_f)/2+size_samps/2);
    TxAntStreams=ifft(ifftshift(TxAntStreams_f,2), size_samps, 2)*size_samps;



    % 
    % for i=1:prm.numTxAntenna
    %     TxAntStreams(i,:)=(TxAntStreams(i,:));
    % end
    %TxAntStreams=TxAntStreams./sqrt(sum(w.^2));


    % figure
    % subplot(2,1,1)
    % plot((abs(fftshift(fft(DAC_output_t)))));
    % ylabel("Amplitude")
    % xlabel("Frequency (MHz)")
    % title("Analog signal freq domain")
    % 
    % subplot(2,1,2)
    % plot((abs(Switch_output_f(1,:))));
    % ylabel("Amplitude")
    % xlabel("Frequency (MHz)")
    % title("Switching output freq domain")



    % figure
    % subplot(2,2,1)
    % plot(abs(fftshift(fft(TxAntStreams(1,:)))))
    % hold on
    % plot(abs((fftshift(fft(mimo_signals(1,:))))))
    % legend('Radiated','Original')
    % 
    % subplot(2,2,2)
    % plot(unwrap(angle((fftshift(fft(TxAntStreams(1,:)))))))
    % hold on
    % plot(unwrap(angle((fftshift(fft(mimo_signals(1,:)))))))
    % legend('Radiated','Original')
    % 
    % subplot(2,2,3)
    % plot(abs(fftshift(fft(TxAntStreams(2,:)))))
    % hold on
    % plot(abs((fftshift(fft(mimo_signals(2,:))))))
    % legend('Radiated','Original')
    % 
    % subplot(2,2,4)
    % plot(unwrap(angle((fftshift(fft(TxAntStreams(2,:)))))))
    % hold on
    % plot(unwrap(angle((fftshift(fft(mimo_signals(2,:)))))))
    % legend('Radiated','Original')
end