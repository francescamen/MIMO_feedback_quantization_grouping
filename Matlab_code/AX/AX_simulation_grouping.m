clear all
clear classes
clear java
close all

addpath('utilities/')
addpath('../beamforming_utilities/')

%% Francesca Meneghello
% This is a modification of the HECompressedBeamformingExample.m provided by MathWorks, Inc.

% Copyright 2023 Francesca Meneghello

%% 802.11ax Compressed Beamforming Packet Error Rate Emulation
%
% This example shows how to measure the packet error rate of a beamformed
% IEEE(R) 802.11ax(TM) high efficiency single user (HE SU) format link with
% different beamforming feedback quantization levels.

%% Emulation Parameters
% 802.11ax specifies only two sets of quantization resolution for single 
% user beamforming. Set the value of index_bf to slect the number of bits
% used to quantize the beamforming feedback angles (phi and psi). 
% To disable compression set index_bf to Inf| This table shows the 
% quantization levels for each value of index_bf:
%
%    index_bf              Compression Configuration
% --------------------------------------------------------
%      0                   NumBitsphi = 4; NumBitspsi = 2
%      1                   NumBitsphi = 6; NumBitspsi = 4
%      2                   No compression
% --------------------------------------------------------

snr           = 20;                  % SNR in dB
index_bf = 1; % Index for the beamforming feedback
N_grouping = 0;

switch index_bf
      case 0
        NumBitsPsi = 2; % Number of bits for psi
        NumBitsPhi = 4; % Number of bits for phi
        disp('End-to-End emulation with compressed beamforming quantization with');
        disp(['Number of Bits for phi = ' num2str(NumBitsPhi) ...
              ' and Number of Bits for psi = ' num2str(NumBitsPsi)]);
      case 1
        NumBitsPsi = 4; % Number of bits for psi
        NumBitsPhi = 6; % Number of bits for phi
        disp('End-to-End emulation with compressed beamforming quantization with');
        disp(['Number of Bits for phi = ' num2str(NumBitsPhi) ...
              ' and Number of Bits for psi = ' num2str(NumBitsPsi)]);
      otherwise
        disp('End-to-End emulation with non-compressed beamforming');
end

%% Waveform Configuration
% s = rng(20);                           % Set RNG seed for repeatability

numUsers    = 1;                   % Number of users

NumTxAnts = 2; % Number of transmit antennas
NumSTSVec = 2;    % Number of space-time streams
NumRxAntsVec = NumSTSVec; % Number of receive antennas
chanBW      = 'CBW160';             % Channel bandwidth
userPos     = 0;               % User positions (inside the Group ID Management frame together with the Membership Status, i.e., GroupID inside the wlanVHTConfig)
mcsVec      = 3;               % MCS per user
apepVec     = 1e3;         % Payload per user, in bytes
chCodingVec = 'LDPC'; % Channel coding per user

cfgHEBase = wlanHESUConfig;
cfgHEBase.ChannelBandwidth = chanBW;      % Channel bandwidth
cfgHEBase.NumSpaceTimeStreams = NumSTSVec;    % Number of space-time streams
cfgHEBase.NumTransmitAntennas = NumTxAnts; % Number of transmit antennas
cfgHEBase.APEPLength = apepVec;                % Payload length in bytes
cfgHEBase.ExtendedRange = false;           % Do not use extended range format
cfgHEBase.Upper106ToneRU = false;          % Do not use upper 106 tone RU
cfgHEBase.PreHESpatialMapping = false;     % Spatial mapping of pre-HE fields
cfgHEBase.GuardInterval = 3.2;             % Guard interval duration
cfgHEBase.HELTFType = 4;                   % HE-LTF compression mode
cfgHEBase.ChannelCoding = chCodingVec;          % Channel coding
cfgHEBase.MCS = mcsVec;                         % Modulation and coding scheme
cfgHEBase.SpatialMapping = 'Custom';       % Custom for beamforming

% Precoding and equalization parameters
precodingType = 'ZF';                % Precoding type; ZF or MMSE
eqMethod      = 'MMSE';                 % Equalization method

%% Parameters
% Files of channel parameters of the users
standard = 'AX';
configuration = '/4x2/';
name_environment = 'simulation';
stream_considered = '1';

%% Folder to save the results
name_save_folder = strcat('../simulation_results/', standard, '/', name_environment);
name_save_folder = strcat(name_save_folder, '/', ...
    'NumTxAnts_', num2str(NumTxAnts), ...
    '_NumRxAntsVec_', num2str(NumSTSVec), ...
    '_BW_', chanBW(4:5), ...
    '_idx_grouping_', num2str(N_grouping), ... 
    '_idx_bf_', num2str(index_bf), ...
    '_snr_', num2str(snr), ...
    '_mcs_', num2str(mcsVec), '/');

if ~isfolder(name_save_folder)
    mkdir(name_save_folder);
% else
%     return
end

%% Channel Configuration
% This example uses a TGax NLOS indoor channel model with delay profile
% Model-B. The Model-B profile is considered NLOS when the distance between
% transmitter and receiver is greater than or equal to five meters. For
% more information, see
% <docid:wlan_ref#mw_43b5900e-69e1-4636-b084-1e72dbd46293 wlanTGaxChannel>.

% Create and configure the TGax channel
chanBW = cfgHEBase.ChannelBandwidth;
tgaxChannel = wlanTGaxChannel;
tgaxChannel.DelayProfile = 'Model-B';
tgaxChannel.NumTransmitAntennas = NumTxAnts;
tgaxChannel.NumReceiveAntennas = NumSTSVec;
tgaxChannel.TransmitReceiveDistance = 5; % Distance in meters for NLOS
tgaxChannel.ChannelBandwidth = chanBW;
tgaxChannel.LargeScaleFadingEffect = 'None';
tgaxChannel.NormalizeChannelOutputs = false;
fs = wlanSampleRate(cfgHEBase);
tgaxChannel.SampleRate = fs;

end_time = 600;
end_iter = 4; 

% Check if already processed entirely
time_idx_expanded = sprintf('%05d', end_time);
name_save_file = strcat(name_save_folder, '/results_time', ...
                num2str(time_idx_expanded), '_iteration', ...
                num2str(end_iter), '.mat');

% if isfile(name_save_file)
%     return
% end

%% Null Data Packet (NDP) Configuration
% Configure the NDP transmission to have data length of zero. Since the NDP
% is used to obtain the channel state information, set the number of
% space-time streams equal to the number of transmit antennas and directly
% map each space-time stream to a transmit antenna.

cfgNDP = cfgHEBase;
cfgNDP.APEPLength = 0;                  % NDP has no data
cfgNDP.NumSpaceTimeStreams = NumTxAnts; % For feedback matrix calculation
cfgNDP.SpatialMapping = 'Direct';       % Each TxAnt carries a STS
cfgNDP.MCS = 0; 

% Generate the null data packet, with no dataAPEP
txNDPSig = wlanWaveformGenerator([],cfgNDP);
NPDSigLen = size(txNDPSig, 1);

fs = wlanSampleRate(cfgHEBase);

%% Define common variables for the loop

tx_ant_vect = linspace(0, NumTxAnts-1, NumTxAnts).';
numPadZeros = 50;

%% Run the emulation

% Get occupied subcarrier indices and OFDM parameters
ofdmInfo = wlanHEOFDMInfo('HE-Data',cfgHEBase);

% Indices to extract fields from the PPDU
ind = wlanFieldIndices(cfgHEBase);

% Account for noise energy in nulls so the SNR is defined per
% active subcarrier
packetSNR = snr-10*log10(ofdmInfo.FFTLength/ofdmInfo.NumTones);

% Set random substream index per iteration to ensure that each
% iteration uses a repeatable set of random numbers
stream = RandStream('combRecursive','Seed',100);
stream.Substream = 1;
RandStream.setGlobalStream(stream);

% bertot=0;
for time_idx=1:end_time
    disp(time_idx)
    % Create an instance of the HE configuration object per sample
    % emulated. This will enable to use parfor.
    cfgHE = cfgHEBase;

    % disp(['time_idx ' num2str(time_idx)]);
    time_idx_expanded = sprintf('%05d', time_idx);

    for num_iter=1:end_iter

        num_iter_expanded = sprintf('%02d', num_iter);
        name_save_file = strcat(name_save_folder, '/results_time', ...
            num2str(time_idx_expanded), '_iteration', ...
            num2str(num_iter_expanded), '.mat');

        % if isfile(name_save_file)
        %      continue
        % end

         %% Sounding
        % Add trailing zeros to allow for channel delay
        txPad = [txNDPSig; zeros(numPadZeros,cfgNDP.NumTransmitAntennas)];

        % Pass through a fading indoor TGax channel
        reset(tgaxChannel); % Reset channel for different realization

        % Sound the emulated channel
        chanOutNDP{1} = tgaxChannel(txPad); 

        % Add AWGN
        rxNDPSig = awgn(chanOutNDP{1}, packetSNR);
            
        %% Channel State Information Feedback
        % Calculate the feedback matrix at the beamformee
        V = heUserBeamformingFeedback(rxNDPSig,cfgNDP,true);  % Nst-by-Nsts-by-Ntx
    
        if isempty(V)
            continue; % Go to next loop iteration
        end
        numST = length(V);

        %% Beamforming feedback grouping, quantization and reconstruction
        %
        V_complete = V;
        if N_grouping~=0
            freq_vector = N_grouping:N_grouping:size(V_complete, 1)/2-N_grouping;
            freq_vector = [flip(-freq_vector), freq_vector];
            freq_vector = freq_vector + size(V_complete, 1)/2;
            V = V_complete(freq_vector, :, :);
        end

        [~,Nc,Nr] = size(V(1,1:NumSTSVec,:));
        if index_bf~=2
            % Find quantized angles of the beamforming feedback matrix
            angidx = bfCompressQuantize(V(:,1:NumSTSVec,:),NumBitsPhi,NumBitsPsi);
    
            % Calculate precoding matrix from the quantized angles at
            % beamformer:
            % Assuming zero delay in transmitting the quantized angles
            % from beamformee to beamformer, the precoding matrix is
            % calculated from the quantized angles and is used in the
            % data transmission of beamformer.
    
            V_rec = bfDecompress(angidx,Nr,Nc,NumBitsPhi,NumBitsPsi); % Nst-by-Nsts-by-Nt
        else
            V_rec = V;
        end

        V_rec_complete = zeros(size(V_complete));
        if N_grouping~=0
            for idx_g = 1:N_grouping
                indices_low = freq_vector(1:(size(freq_vector, 2)/2)) + idx_g - 1;
                V_rec_complete(indices_low, :, :) = V_rec(1:(size(freq_vector, 2))/2, :, :);

                index_center = freq_vector(size(freq_vector, 2)/2) + N_grouping + idx_g - 1;
                V_rec_complete(index_center, :, :) = V_rec((size(freq_vector, 2))/2, :, :);

                indices_up = freq_vector(size(freq_vector, 2)/2+1:size(freq_vector, 2)) + idx_g - 1;
                V_rec_complete(indices_up, :, :) = V_rec((size(freq_vector, 2))/2+1:size(freq_vector, 2), :, :);
            end
            V_rec_complete(1:freq_vector(1), :, :) = repmat(V_rec(1, :, :), freq_vector(1), 1);
            V_rec_complete(freq_vector(end)+1:size(V_complete, 1), :, :) = repmat(V_rec(end, :, :), freq_vector(1), 1);
        else
            V_rec_complete = V_rec;
        end

        V_rec_complete = V_rec_complete(1:size(V_complete, 1), :, :);

        precodingMat = V_rec_complete(:,1:NumSTSVec,:);
     
        %% Zero-forcing or MMSE precoding solution
        if strcmp(precodingType, 'ZF')
            delta = 0; % Zero-forcing
        else
            delta = (NumTxAnts/(10^(snr/10))) * eye(NumTxAnts); % MMSE
        end
        for i = 1:numST
            % Channel inversion precoding
            v = squeeze(precodingMat(i,:,:)).';
            precodingMat(i,:,:) = (v/(v'*v + delta)).';
        end

        % Set the spatial mapping based on the precoding matrix
        cfgHE.SpatialMappingMatrix = precodingMat;

        %% Data Transmission
        % Generate payload
        psduLength = getPSDULength(cfgHE);
        % txDataBits = randi([0 1],cfgHE.APEPLength*8,1,'int8');
        txDataBits = randi([0 1],psduLength*8,1);
            
        % Pad payload with zeros to form a PSDU
        % psduDataBits = [txDataBits; ...
        %     zeros((psduLength-cfgHE.APEPLength)*8,1,'int8')];
        psduDataBits = txDataBits;
      
        %%
        % Using the format configuration, |cfgHE|, with the steering matrix, to
        % generate the multiuser HE waveform.

        txSig = wlanWaveformGenerator(psduDataBits,cfgHE);
    
        % Add trailing zeros to allow for channel delay
        txPad = [txSig; zeros(numPadZeros,cfgHE.NumTransmitAntennas)];
    
        %%
        % As we restart the channel, we want
        % it to re-process the NDP before the waveform so as to accurately mimic
        % the channel continuity. Only the waveform portion of the channel's output
        % is extracted for the subsequent processing of each user.
        
        % Transmit through the emulated channel for all users, with 10 all-zero
        % samples appended to account for channel filter delay
        
        chanOut{1} = tgaxChannel([txSig; zeros(numPadZeros,NumTxAnts)]); 
        
        % Extract the waveform output for each user
        % chanOut = cellfun(@(x) x(NPDSigLen+numPadZeros+1:end,:),chanOut,'UniformOutput',false);
        
        % Add AWGN
        rxSig = awgn(chanOut{1}, packetSNR);
    
        %% Data Recovery Per User
        % Packet detect and determine coarse packet offset
        coarsePktOffset = wlanPacketDetect(rxSig,chanBW);
        if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
            continue; % Go to next loop iteration
        end
    
        % Extract L-STF and perform coarse frequency offset correction
        lstf = rxSig(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
        coarseFreqOff = wlanCoarseCFOEstimate(lstf,chanBW);
        rxSig = frequencyOffset(rxSig,fs,-coarseFreqOff);
    
        % Extract the non-HT fields and determine fine packet offset
        nonhtfields = rxSig(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
        finePktOffset = wlanSymbolTimingEstimate(nonhtfields,chanBW);
    
        % Determine final packet offset
        pktOffset = coarsePktOffset+finePktOffset;
    
        % If packet detected outwith the range of expected delays from
        % the channel modeling; packet error
        if pktOffset>numPadZeros
            continue; % Go to next loop iteration
        end
    
        % Extract L-LTF and perform fine frequency offset correction
        rxLLTF = rxSig(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
        fineFreqOff = wlanFineCFOEstimate(rxLLTF,chanBW);
        rxSig = frequencyOffset(rxSig,fs,-fineFreqOff);
    
        % HE-LTF demodulation and channel estimation
        rxHELTF = rxSig(pktOffset+(ind.HELTF(1):ind.HELTF(2)),:);
        heltfDemod = wlanHEDemodulate(rxHELTF,'HE-LTF',cfgHE);
        [chanEst,pilotEst] = wlanHELTFChannelEstimate(heltfDemod,cfgHE);
    
        % Data demodulate
        rxData = rxSig(pktOffset+(ind.HEData(1):ind.HEData(2)),:);
        demodSym = wlanHEDemodulate(rxData,'HE-Data',cfgHE);
    
        % Pilot phase tracking
        % Average single-stream pilot estimates over symbols (2nd dimension)
        pilotEstTrack = mean(pilotEst,2);
        demodSym = wlanHETrackPilotError(demodSym,pilotEstTrack,cfgHE,'HE-Data');
    
        % Estimate noise power in HE fields
        nVarEst = wlanHEDataNoiseEstimate(demodSym(ofdmInfo.PilotIndices,:,:),pilotEstTrack,cfgHE);
    
        % Extract data subcarriers from demodulated symbols and channel
        % estimate
        demodDataSym = demodSym(ofdmInfo.DataIndices,:,:);
        chanEstData = chanEst(ofdmInfo.DataIndices,:,:);
    
        % Equalization and STBC combining
        [eqDataSym,csi] = wlanHEEqualize(demodDataSym,chanEstData,nVarEst,cfgHE,'HE-Data');
    
        % Recover data
        rxDataBits = wlanHEDataBitRecover(eqDataSym,nVarEst,csi,cfgHE,'LDPCDecodingMethod','norm-min-sum');
        
        % Compare recovered bits against APEPLength information bits
        idx = (1:cfgHE.APEPLength*8).';
        [~,ber] = biterr(txDataBits(idx),rxDataBits(idx));
    
        %% Save data
        cellV{1} = V;
        cellV_rec{1} = V_rec;
        dict.packetSNR = packetSNR;
        dict.V_mat = cellV;
        dict.V_mat_rec = cell(cellV_rec);
        dict.ber = ber;
        % bertot = bertot+ber;
        disp(ber)
        
        save(name_save_file, '-struct', 'dict');
    end
end
