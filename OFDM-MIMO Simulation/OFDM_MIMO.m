function BER=OFDM_MIMO(channel_f,BF)
%% test simulation 
rng(1)

% Beamforming Type 
prm.BF=BF; % DB: Digital Beamforming, HB: Hybrid Beamforming, SB: Switching Beamforming, PB= Phase-Shifter Beamforming
% Waveform params
prm.N_OFDM_SYMS             = 500;         % Number of OFDM symbols
prm.MOD_ORDER               =  1;          % Modulation order in power of 2 (1/2/4/6 = BSPK/QPSK/16-QAM/64-QAM)
% TX_SCALE                = 1.0;         % Scale for Tx waveform ([0:1])
prm.TX_SPATIAL_STREAM_SHIFT= 3; 

%Noise parameters
prm.noise_flag=0;
prm.snr=10;

% OFDM params
prm.SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
prm.SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
prm.N_SC                    = 64;                                     % Number of subcarriers
prm.CP_LEN                  = 16;                                     % Cyclic prefix length
prm.N_DATA_SYMS             = prm.N_OFDM_SYMS * length(prm.SC_IND_DATA);      % Number of data symbols (one per data-bearing subcarrier per OFDM symbol)

prm.channel_coding = .5; % coding rate
prm.trellis_end_length = 8; % bits for trellis to end
prm.trel=poly2trellis(7,[171, 133]);

prm.number_of_bits= (prm.N_DATA_SYMS * prm.MOD_ORDER - 2*prm.trellis_end_length) * prm.channel_coding;

prm.pilots = [1 1 -1 1].';

% LTS
prm.lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
prm.lts_t = ifft(prm.lts_f, 64);

% % STS
% sts_f = zeros(1,64);
% sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
% sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
% sts_t = ifft(sqrt(13/6).*sts_f, 64);
% sts_t = sts_t(1:16);
% sts_t_rep = repmat(sts_t, 1, 30);

prm.channel_estimation=0; % 0: adding noise, 1: estimation using lts
prm.channel_noise_metric=0.0001;
% User and Antenna parameters

[prm.numUsers, prm.numTxAntenna, ~]=size(channel_f);
channel_f=channel_f(1:prm.numUsers,1:prm.numTxAntenna,:);

% Converting freq. domain channel to time domain
[channel_t]=TimeDomainChannel(prm,channel_f);

% Generating switching matrix
if(prm.BF=="DB")
    phi=eye(prm.numTxAntenna);
elseif(prm.BF=="SB")
    phi = generate_full_rank_binary_matrix(prm);
elseif(prm.BF=="PB")
    phi = exp(1i*angle(channel_f(:,:,33)'));
end

% Generating prerambles for different precodings
[preamble]=generate_preamble(prm,phi);

% Sending preambles
mimo_signals= preamble;

% Generating antena streams for different precodings
if(prm.BF=="SB" || prm.BF=="PB")
    TxAntStreams=BinarySwitching(prm,mimo_signals,phi);
elseif(prm.BF=="DB")
    TxAntStreams=mimo_signals;
end

% Passing preambles through the channel for channel estimation
[RxAntStreams]=ApplyChannel(prm, channel_t, TxAntStreams);


% Adding noise to the received signals
if(prm.noise_flag==1)
    RxAntStreamsNoisy=AWGN(prm, RxAntStreams);
else
    RxAntStreamsNoisy=RxAntStreams;
end

% Channel estimation using lts sequences
if(prm.channel_estimation==0) % Generating estimated channel by adding noise to real channel
    noise_metric=prm.channel_noise_metric;
    channel_noise=complex(sqrt(noise_metric/2)*randn(size(channel_f)),sqrt(noise_metric/2)*randn(size(channel_f)));
    est_channel_f=channel_noise+channel_f;
elseif(prm.channel_estimation==1)
    [est_channel_f]=ChannelEstimation(prm,RxAntStreamsNoisy);
end

% % Ploting the estimated channel and physical channel
% figure
% subplot(8,1,1)
% plot(abs(squeeze(channel_f(1,1,:))))
% hold on
% plot(abs(squeeze(est_channel_f(1,1,:))))
% hold off
% title("User 1- Tx 1")
% legend("Real Channel","Estimated Channel")
% 
% subplot(8,1,2)
% plot(abs(squeeze(channel_f(2,1,:))))
% hold on
% plot(abs(squeeze(est_channel_f(2,1,:))))
% hold off
% title("User 1- Tx 2")
% legend("Real Channel","Estimated Channel")
% 
% subplot(8,1,3)
% plot(abs(squeeze(channel_f(3,1,:))))
% hold on
% plot(abs(squeeze(est_channel_f(3,1,:))))
% hold off
% title("User 1- Tx 3")
% legend("Real Channel","Estimated Channel")
% 
% subplot(8,1,4)
% plot(abs(squeeze(channel_f(4,1,:))))
% hold on
% plot(abs(squeeze(est_channel_f(4,1,:))))
% hold off
% title("User 1- Tx 4")
% legend("Real Channel","Estimated Channel")

% subplot(6,1,4)
% plot(abs(squeeze(channel_f(2,1,:))))
% hold on
% plot(abs(squeeze(est_channel_f(2,1,:))))
% hold off
% title("User 2- Tx 1")
% legend("Real Channel","Estimated Channel")
% 
% subplot(6,1,5)
% plot(abs(squeeze(channel_f(2,2,:))))
% hold on
% plot(abs(squeeze(est_channel_f(2,2,:))))
% hold off
% title("User 2- Tx 2")
% 
% subplot(6,1,6)
% plot(abs(squeeze(channel_f(2,3,:))))
% hold on
% plot(abs(squeeze(est_channel_f(2,3,:))))
% hold off
% title("User 2- Tx 3")
% legend("Real Channel","Estimated Channel")

% figure
% subplot(4,1,1)
% plot(phase(squeeze(channel_f(1,1,:))))
% hold on
% plot(phase(squeeze(est_channel_f(1,1,:))))
% hold off
% title("User 1- Tx 1")
% 
% subplot(4,1,2)
% plot(phase(squeeze(channel_f(1,2,:))))
% hold on
% plot(phase(squeeze(est_channel_f(1,2,:))))
% hold off
% title("User 1- Tx 2")
% 
% 
% subplot(4,1,3)
% plot(phase(squeeze(channel_f(1,3,:))))
% hold on
% plot(phase(squeeze(est_channel_f(1,3,:))))
% hold off
% title("User 1- Tx 3")

%% Doing the precoding and switching once more based on estimated channels

% Generating OFDM symbols for different precodings
[tx_data, tx_payload_vec, tx_syms]=generate_ofdm_data(prm,est_channel_f, phi);

% Preappending preamble
mimo_signals=[preamble tx_payload_vec];

% Generating antenna streams for different precodings
if(prm.BF=="SB" || prm.BF=="PB")
    TxAntStreams=BinarySwitching(prm,mimo_signals,phi);
% TxAntStreams=;
elseif(prm.BF=="DB")
    TxAntStreams=mimo_signals;
end

% Passing mimo signals through the channel for channel estimation
[RxAntStreams]=ApplyChannel(prm,channel_t,TxAntStreams);

% Adding noise to the received signals
if(prm.noise_flag==1)
    RxAntStreamsNoisy=AWGN(prm, RxAntStreams);
else
    RxAntStreamsNoisy=RxAntStreams;
end


% receiver decoding
decodedBits=OFDMReceiver(prm,RxAntStreamsNoisy);

% displaying the BER per user
[~,BER]=biterr(tx_data, decodedBits, 'row-wise');
end
