clear all
clear classes
clear java
close all

addpath('utilities/')
addpath('../beamforming_utilities/')

%% Francesca Meneghello
% This is a modification of the WINNERVHTMUMIMOExample.m provided by MathWorks, Inc.

% Copyright 2023 Francesca Meneghello

%% 802.11ac Compressed Beamforming Packet Error Rate Emulation
%
% This example shows how to measure the packet error rate of a beamformed
% IEEE(R) 802.11ac(TM) very high throughput single user (VHT SU) format link with
% different beamforming feedback quantization levels.

%% Emulation Parameters
% 802.11ac specifies only two sets of quantization resolution for single 
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
chanBW      = 'CBW80';             % Channel bandwidth
userPos     = 0;               % User positions (inside the Group ID Management frame together with the Membership Status, i.e., GroupID inside the wlanVHTConfig)
mcsVec      = 3;               % MCS per user
apepVec     = 1e3;         % Payload per user, in bytes
chCodingVec = 'LDPC'; % Channel coding per user

cfgVHTBase = wlanVHTConfig;
cfgVHTBase.ChannelBandwidth = chanBW;      % Channel bandwidth
cfgVHTBase.NumUsers = numUsers;
cfgVHTBase.NumSpaceTimeStreams = NumSTSVec;    % Number of space-time streams
cfgVHTBase.NumTransmitAntennas = NumTxAnts; % Number of transmit antennas
cfgVHTBase.APEPLength = apepVec;                % Payload length in bytes
cfgVHTBase.ChannelCoding = chCodingVec;          % Channel coding
cfgVHTBase.MCS = mcsVec;                         % Modulation and coding scheme
cfgVHTBase.SpatialMapping = 'Custom';       % Custom for beamforming
cfgVHTBase.GroupID = 0; % GroupID must be 0 or 63 when NumUsers is 1 and between 1 and 62 inclusive when NumUsers is 2, 3 or 4.
cfgVHTBase.UserPositions = userPos;

% Precoding and equalization parameters
precodingType = 'ZF';                % Precoding type; ZF or MMSE
eqMethod      = 'MMSE';                 % Equalization method

%% Parameters
standard = 'AC';
configuration = '/4x2/';
folder_environment = 'Anechoic_2';
name_environment = 'anechoic2';
stream_considered = '1';

%% Folder to save the results
name_save_folder = strcat('../emulation_results/', standard, '/', name_environment);
name_save_folder = strcat(name_save_folder, '/', ...
    'NumTxAnts_', num2str(NumTxAnts), ...
    '_NumRxAntsVec_', num2str(NumRxAntsVec), ...
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

%% Upload the channel obtained from ray tracing of experimental data
folder_name_base = strcat(standard, configuration, folder_environment);
folder_name_data = strcat('../../Data/CSI/', folder_name_base, ...
    '/CSI/Antenna_Separated/decomposed/');
file_name_base = {strcat(folder_name_data, name_environment, '_', ...
    chanBW(4:5), '_', lower(standard), '_ss', stream_considered)};

amplitudes_users = cell(numUsers, 1);
aoas_users = cell(numUsers, 1);
toas_users = cell(numUsers, 1);

for us=1:numUsers
    % Upload the experimental channel parameters
            
    amplitudes_users{us} = importdata(strcat(file_name_base{us}, '_paths_amplitude_list.m'));
    aoas_users{us} = importdata(strcat(file_name_base{us}, '_paths_aoa_list.m'));
    toas_user = importdata(strcat(file_name_base{us}, '_paths_toa_list.m'));
    for ii=1:size(toas_user, 2)
        toas_user{ii} = toas_user{ii} - toas_user{ii}(:, 1);
    end
    toas_users{us} = toas_user;

end

end_time = size(amplitudes_users{1}, 2); 
disp(strcat('end_time ', num2str(end_time)))
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

cfgNDP = cfgVHTBase;
cfgNDP.APEPLength = 0;                  % NDP has no data
cfgNDP.NumSpaceTimeStreams = NumTxAnts; % For feedback matrix calculation
cfgNDP.SpatialMapping = 'Direct';       % Each TxAnt carries a STS
cfgNDP.MCS = 0; 
cfgNDP.GroupID = 0;
cfgNDP.NumUsers = 1;

% Generate the null data packet, with no dataAPEP
txNDPSig = wlanWaveformGenerator([],cfgNDP);
NPDSigLen = size(txNDPSig, 1);

fs = wlanSampleRate(cfgVHTBase);

%% Define common variables for the loop

tx_ant_vect = linspace(0, NumTxAnts-1, NumTxAnts).';
numPadZeros = 50;

%% Run the emulation

% Get occupied subcarrier indices and OFDM parameters
ofdmInfo = wlanVHTOFDMInfo('VHT-Data',cfgVHTBase);

% Indices to extract fields from the PPDU
ind = wlanFieldIndices(cfgVHTBase);

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
    % Create an instance of the VHT configuration object per sample
    % emulated. This will enable to use parfor.
    cfgVHT = cfgVHTBase;

    % disp(['time_idx ' num2str(time_idx)]);
    time_idx_expanded = sprintf('%05d', time_idx);

    for num_iter=1:end_iter

        num_iter_expanded = sprintf('%02d', num_iter);
        name_save_file = strcat(name_save_folder, '/results_time', ...
            num2str(time_idx_expanded), '_iteration', ...
            num2str(num_iter_expanded), '.mat');

        if isfile(name_save_file)
             continue
        end

        %% Sounding
        % Add trailing zeros to allow for channel delay
        txPad = [txNDPSig; zeros(numPadZeros,cfgNDP.NumTransmitAntennas)];

        chanOutNDP = cell(numUsers, 1);
        chanFiltUsers = cell(numUsers, 1);
        amplitudesChannelUsers = cell(numUsers, 1);
    
        for us=1:numUsers
            % Upload the experimental channel parameters
                
            amplitudes = amplitudes_users{us}{time_idx};
            aoas = aoas_users{us}{time_idx};
            toas = toas_users{us}{time_idx};
            
            amplitudes_channel = 20*log10(abs(amplitudes));

            chanFilt = comm.MIMOChannel( ...
                            'SampleRate',fs, ...
                            'PathDelays',toas, ...
                            'AveragePathGains',amplitudes_channel, ...
                            'MaximumDopplerShift',0, ...
                            'SpatialCorrelationSpecification','none', ...
                            'NumTransmitAntennas',NumTxAnts, ...
                            'NumReceiveAntennas',NumRxAntsVec, ...
                            'ChannelFiltering',true);
                
            chanFiltUsers{us} = chanFilt;
            amplitudesChannelUsers{us} = amplitudes_channel;
            
            % Sound the emulated channel
            chanOutNDP{us} = chanFilt(txPad); 
        end
        % Add AWGN
        % rxNDPSig = cellfun(@awgn,chanOutNDP, ...
        %     num2cell(packetSNR*ones(numUsers,1)),'UniformOutput',false);
        rxNDPSig = cell(1);
        rxNDPSig{1} = awgn(chanOutNDP{1}, packetSNR);
                    
        %% Channel State Information Feedback
        % Calculate the feedback matrix at the beamformee
        V_mat = cell(numUsers,1);  % V matrix
        csi_feed = cell(numUsers,1); % chanEstMinusCSD: channel estimate without cyclic prefix effect
        corrupted = false;
        for uIdx = 1:numUsers
            % Compute the feedback matrix based on received signal per user
            [V, csi_feed{uIdx}] = vhtCSIFeedback(rxNDPSig{uIdx}, ...
                cfgNDP,uIdx,NumSTSVec(uIdx));  % Nst-by-Nr-by-Nsts
            if isempty(V)
                corrupted = true; % Go to next loop iteration
                continue; 
            end
            V = permute(V,[1 3 2]);
            V_mat{uIdx} = V;
        end
        if corrupted
            continue; % Go to next loop iteration
        end

        numST = length(V_mat{1});         % Number of subcarriers
        
        %% Beamforming feedback grouping, quantization and reconstruction
        %
        V_complete_mat = V_mat;
        if N_grouping~=0
            for uIdx = 1:numUsers
                V_complete = V_complete_mat{uIdx};
                freq_vector = - size(V_complete, 1)/2+2:N_grouping:size(V_complete, 1)/2-N_grouping;
                freq_vector((size(freq_vector, 2)+1)/2) = [];
                freq_vector = freq_vector + size(V_complete, 1)/2 - 1;
                V = V_complete(freq_vector, :, :);
                V_mat{uIdx} = V; % Nst-by-Nsts-by-Nt
            end
        end

        V_mat_rec = cell(numUsers,1);  % V tilde matrix reconstructed
        if index_bf~=2
            for uIdx = 1:numUsers
                V = V_mat{uIdx};
                % Find quantized angles of the beamforming feedback matrix
                % V must be Nst-by-Nc-by-Nr
                angidx = bfCompressQuantize(V(:,1:NumSTSVec,:),NumBitsPhi,NumBitsPsi);
        
                % Calculate precoding matrix from the quantized angles at
                % beamformer:
                % Assuming zero delay in transmitting the quantized angles
                % from beamformee to beamformer, the precoding matrix is
                % calculated from the quantized angles and is used in the
                % data transmission of beamformer.
        
                [~,Nc,Nr] = size(V(1,1:NumSTSVec,:));
                V_rec = bfDecompress(angidx,Nr,Nc,NumBitsPhi,NumBitsPsi);
    
                V_mat_rec{uIdx} = V_rec; % Nst-by-Nsts-by-Nt
            end
        else
            V_mat_rec = V_mat;
        end

        precodingMat = zeros(numST, sum(NumSTSVec), NumTxAnts);
        for uIdx = 1:numUsers
            V_rec = V_mat_rec{uIdx};
            V_rec_complete = zeros(size(V_complete_mat{uIdx}));
            if N_grouping~=0
                for idx_g = 1:N_grouping
                    indices_low = freq_vector(1:(size(freq_vector, 2)/2)) + idx_g - 1;
                    V_rec_complete(indices_low, :, :) = V_rec(1:(size(freq_vector, 2))/2, :, :);
    
                    index_center = freq_vector(size(freq_vector, 2)/2) + N_grouping + idx_g - 1;
                    V_rec_complete(index_center, :, :) = V_rec((size(freq_vector, 2))/2, :, :);
    
                    indices_up = freq_vector(size(freq_vector, 2)/2+1:size(freq_vector, 2)) + idx_g - 1;
                    V_rec_complete(indices_up, :, :) = V_rec((size(freq_vector, 2))/2+1:size(freq_vector, 2), :, :);
                end
                V_rec_complete = V_rec_complete(1:numST, :, :);
            else
                V_rec_complete = V_rec;
            end
    
            stsIdx = sum(NumSTSVec(1:uIdx-1))+(1:NumSTSVec(uIdx));
            precodingMat(:,stsIdx, :) = V_rec_complete;     % Nst-by-Nsts-by-Nt
        end
        
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
        cfgVHT.SpatialMappingMatrix = precodingMat; % Nst-by-Nsts-by-Nt
       
        %% Data Transmission
        % Create data sequences, one for each user
        txDataBits = cell(numUsers,1);
        psduDataBits = cell(numUsers,1);
        for uIdx = 1:numUsers
            % Generate payload for each user
            psduLength = cfgVHT.PSDULength(uIdx);
            % txDataBits{uIdx} = randi([0 1],cfgVHT.APEPLength(uIdx)*8,1,'int8');
            txDataBits{uIdx} = randi([0 1],psduLength*8,1);
            
            % Pad payload with zeros to form a PSDU
            % psduDataBits{uIdx} = [txDataBits{uIdx}; ...
            %     zeros((cfgVHT.PSDULength(uIdx)-cfgVHT.APEPLength(uIdx))*8,1,'int8')];
            psduDataBits{uIdx} = txDataBits{uIdx};
        end
        
        %%
        % Using the format configuration, |cfgVHT|, with the steering matrix, to
        % generate the multiuser VHT waveform.
        
        txSig = wlanWaveformGenerator(psduDataBits,cfgVHT);
        
        % Add trailing zeros to allow for channel delay
        txPad = [txSig; zeros(numPadZeros,cfgVHT.NumTransmitAntennas)];

        %%
        % As we restart the channel, we want
        % it to re-process the NDP before the waveform so as to accurately mimic
        % the channel continuity. Only the waveform portion of the channel's output
        % is extracted for the subsequent processing of each user.
        
        % Transmit through the emulated channel for all users, with 10 all-zero
        % samples appended to account for channel filter delay
        
        chanOut = cell(numUsers, 1);
        rxSig = cell(numUsers, 1);
        for us=1:numUsers
            % chanOut{us} = chanFiltUsers{us}([txNDPSig; zeros(numPadZeros,NumTxAnts); ...
            %     txSig; zeros(numPadZeros,NumTxAnts)]); 
            chanOut{us} = chanFiltUsers{us}([txSig; zeros(numPadZeros,NumTxAnts)]); 
            
            % Add AWGN
            rxSig{us} = awgn(chanOut{us}, packetSNR);
        end
        
        % Extract the waveform output for each user
        % chanOut = cellfun(@(x) x(NPDSigLen+numPadZeros+1:end,:),chanOut,'UniformOutput',false);
        
        % Add AWGN
        % rxSig = cellfun(@awgn,chanOut, ...
        %     num2cell(packetSNR*ones(numUsers,1)),'UniformOutput',false);
        
        %% Data Recovery Per User
        %
        % The receive signals for each user are processed individually. The example
        % assumes that there are no front-end impairments and that the transmit
        % configuration is known by the receiver for simplicity.
        %
        % A user number specifies the user of interest being decoded for the
        % transmission. This is also used to index into the vector properties of
        % the configuration object that are user-specific.
                
        % Single-user receivers recover payload bits
        rxDataBits = cell(numUsers,1);
        GMatrices = cell(numUsers,1);
        corrupted = false;
        for uIdx = 1:numUsers
            % rxNSig = rxSig{uIdx}(chanDelay(uIdx)+1:end, :);
            rxNSig = rxSig{uIdx};
            coarsePktOffset = wlanPacketDetect(rxNSig,chanBW);
            if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
                corrupted = true;
                continue; % Go to next loop iteration
            end

            % Extract L-STF and perform coarse frequency offset correction
            lstf = rxNSig(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
            coarseFreqOff = wlanCoarseCFOEstimate(lstf,chanBW);
            rxNSig = frequencyOffset(rxNSig,fs,-coarseFreqOff);

            % Extract the non-HT fields and determine fine packet offset
            nonhtfields = rxNSig(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
            finePktOffset = wlanSymbolTimingEstimate(nonhtfields,chanBW);

            % Determine final packet offset
            pktOffset = coarsePktOffset+finePktOffset;

            % If packet detected outwith the range of expected delays from
            % the channel modeling; packet error
            if pktOffset>numPadZeros
                corrupted = true;
                continue; % Go to next loop iteration
            end

            % Extract L-LTF and perform fine frequency offset correction
            rxLLTF = rxNSig(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
            fineFreqOff = wlanFineCFOEstimate(rxLLTF,chanBW);
            rxNSig = frequencyOffset(rxNSig,fs,-fineFreqOff);

            % Estimate noise power in VHT fields
            lltf = rxNSig(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
            demodLLTF = wlanLLTFDemodulate(lltf,chanBW);
            nVar = helperNoiseEstimate(demodLLTF,chanBW,sum(NumSTSVec(uIdx)));
            % fprintf('noise %.5f\n', nVar);
            
            % VHT-LTF demodulation and channel estimation
            rxVHTLTF  = rxNSig(pktOffset+(ind.VHTLTF(1):ind.VHTLTF(2)),:);
            % demodVHTLTF = wlanVHTLTFDemodulate(rxVHTLTF,chanBW,NumSTSVec(uIdx));
            demodVHTLTF = wlanVHTLTFDemodulate(rxVHTLTF,cfgVHT);
            % chanEst = wlanVHTLTFChannelEstimate(demodVHTLTF,chanBW,NumSTSVec(uIdx)); % H tilde Nst-by-Nsts-by-Nr
            chanEst = wlanVHTLTFChannelEstimate(demodVHTLTF,cfgVHT); % H tilde Nst-by-Nsts-by-Nr
        
            % Recover information bits in VHT Data field
            rxVHTData = rxNSig(pktOffset+(ind.VHTData(1):ind.VHTData(2)),:);
            [rxDataBits{uIdx},~,eqsym,GMatrices{uIdx}] = wlanVHTDataRecover(rxVHTData, ...
                chanEst,nVar,cfgVHT,uIdx, ...
                'EqualizationMethod',eqMethod,'PilotPhaseTracking','None', ...
                'LDPCDecodingMethod','layered-bp','MaximumLDPCIterationCount',6);
        end
        if corrupted
            continue
        end
        
        %%
        % Per-stream equalized symbol constellation plots validate the simulation
        % parameters and convey the effectiveness of the technique. Note the
        % discernible 16QAM, 64QAM and QPSK constellations per user as specified on
        % the transmit end. Also observe the EVM degradation over the different
        % streams for an individual user. This is a representative characteristic
        % of the channel inversion technique.
        %
        % The recovered data bits are compared with the transmitted payload bits to
        % determine the bit error rate.
        
        % Compare recovered bits against per-user APEPLength information bits
        ber = inf(1, numUsers);
        for uIdx = 1:numUsers
            idx = (1:cfgVHT.APEPLength(uIdx)*8).';
            [~,ber(uIdx)] = biterr(txDataBits{uIdx}(idx),rxDataBits{uIdx}(idx));
            % disp(['Bit Error Rate for User ' num2str(uIdx) ': ' num2str(ber(uIdx))]);
        end
        
        % rng(s); % Restore RNG state
        
        %% Save data
        dict.packetSNR = packetSNR;
        dict.V_mat = V_mat;
        dict.V_mat_rec = V_mat_rec;
        dict.ber = ber;
        % bertot = bertot+ber;
        
        save(name_save_file, '-struct', 'dict');
    end
end

% for i=1:NumTxAnts
%     figure(); 
%     for ii=1:numUsers
%         plot(abs(mat_rec_no_attack{ii}(:, i, 1)), 'DisplayName', strcat('user', num2str(ii)))
%         hold on;
%     end
%     legend();
%     title(['V matrix attack - TX ', num2str(i)])
% end
% 
% for i=1:NumTxAnts
%     figure(); 
%     for ii=1:numUsers
%         plot(abs(mat_rec{ii}(:, i, 1)), 'DisplayName', strcat('user', num2str(ii)))
%         hold on;
%     end
%     legend();
%     title(['V matrix attack - TX ', num2str(i)])
% end
