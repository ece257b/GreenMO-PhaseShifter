function delayed_sig = delayseq(curr_sig, frac_del)
    nfft = numel(curr_sig);
    fbins = 2 * pi * ((0:1:nfft-1) - floor(nfft/2)) / nfft;
    X = fftshift(fft(curr_sig, nfft));
    delayed_sig = ifft(ifftshift(X .* exp(-1j * frac_del * fbins)));
end