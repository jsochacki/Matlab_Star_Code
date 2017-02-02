%%%%%%%%%%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

%SELECT ADAPTIVE ALGORITHM
ALG='rls';

%SET UP PA
pa_coefficients=[(1.1500+j*0.7715) (0.1482-j*1.6597) (-0.5407+j*1.2630) (0.3631-j*0.4015) (-0.0903+j*0.0342)];
GINIT=0.1;

%SET UP FILTER
NSYMBOLS_LONG_FILTER=24;
%24 is fine for 25 percent rolloff
%4*24 is fine for 5 percent rolloff
USAMPR=8; ROLLOFF=0.25; ORDER=USAMPR*NSYMBOLS_LONG_FILTER;
SYMBOL_RATE=1; Fc=SYMBOL_RATE/2;
h=(firrcos(ORDER,Fc,ROLLOFF,USAMPR,'rolloff','sqrt'));

%SET UP FRAME
%NSYMBOLS=2994*10; %10 D&I++ Frames
%as the length of the transmitted sequence will be
%seqlen=length([upsample(ones(1,NSYMBOLS+2*length(PS)),USAMPR) ones(1,length(h)-1)])
%seqlen=length([upsample(ones(1,NSYMBOLS+(2*length(PS))),USAMPR) ones(1,NSYMBOLS_LONG_FILTER*USAMPR)])
%seqlen=((NSYMBOLS+(2*length(PS)))*USAMPR)+(NSYMBOLS_LONG_FILTER*USAMPR)
%if mod(seqlen/NTAPS_OR_PREVEQ,1)~=0
%the short loop will throw out the last batch of samples as
%1~=mod(seqlen/NTAPS_OR_PREVEQ,1).  This is not recommended as your
%processing for the last frame will be done with an old channel model.
%To avoid this, NSYMBOLS is carefully selected below

N_PILOT_SYMBOLS=power(2,nextpow2(ceil(NSYMBOLS_LONG_FILTER/USAMPR)+1));
%must be greater than or equal to the number of symbols the filter takes to
%kick in divided by the upsample rate to prevent erroneous ber results due
%to filter roll in and roll out distortion (i.e. due to non continuous
%signal with the first
% (filterlength-1)/2 symbols not actually being fully filtered and the last
% (filterlength-1)/2 symbols not being fully filtered
%The +1 is just for safety
%The 2^nextpow2 is added to make the NSYMBOLS calculation below always work

%NSYMBOLS Selected to make mod(seqlen/NTAPS_OR_PREVEQ,1)=0
NSYMBOLS=2201-(N_PILOT_SYMBOLS*2); %(ALL USAMPR OK)
                                                     
%SET UP CONSTELLATION
MOD_COD=1;
[complex_mapping Binary_Alphabet Decimal_Alphabet BITS_PER_WORD PARR_Remapping_Vector]=dvbs2_CBAM(MOD_COD);
Binary_Alphabet=custom_data_stream_to_words(Binary_Alphabet,BITS_PER_WORD);

%SET UP TRANSMITTER
AGC_SET_POWER_TX_dBm=10;

%SET UP RECEIVER
AGC_SET_POWER_RX_dBm=10; RECEIVER_AGC_OVERHEAD_DB=2;

%SET UP FOR BER
ber_iteration=1; ERROR_LIMIT=100; EbNo=15; ERRORS=0;

%CHANNEL MODEL SETUP
NTAPS_OR_PREVEQ=25;
switch ALG
    case 'lms'
        %For use with LMS
        mu=0.01; cm_itterations=10;
    case 'nlms'
        %For use with LMS
        mu=0.01; cm_itterations=10;
    case 'rls'
        %For use with RLS
        mu=0.9; Previous_P=[]; epsilon=1e-6; cm_itterations=3;
    otherwise        
end
DATAPATH_TX_TO_RX_POWER_RATIO=40; DATAPATH_TX_TO_RX_VOLTAGE_RATIO=power(10,DATAPATH_TX_TO_RX_POWER_RATIO/20);
SELF_TX_TO_RX_POWER_RATIO=20; SELF_TX_TO_RX_VOLTAGE_RATIO=power(10,SELF_TX_TO_RX_POWER_RATIO/20);
%%%%%simulated_channel_model=[linspace(0,1,13)+j*linspace(1,0,13) linspace(1,0,12)+j*linspace(0,1,12)]./25;
simulated_channel_model_initial_vector=[0.0145981258594804 - 0.000768322413656866i,0.0178719299204647 - 0.000940627890550775i,0.0144289845280510 - 0.000759420238318473i,0.0201719061149929 - 0.00106167926921015i,-0.0219921749080518 + 0.00115748288989746i,0.675441434652535 - 0.0355495491922387i,-0.0219921749080518 + 0.00115748288989746i,0.0201719061149929 - 0.00106167926921015i,0.0144289845280510 - 0.000759420238318473i,0.0178719299204647 - 0.000940627890550775i,0.0145981258594804 - 0.000768322413656866i];
%%%%%simulated_channel_model=simulated_channel_model_initial_vector.*AGC_1_Ohm_System(lconv(h,upsample(simulated_channel_model_initial_vector,USAMPR),'full'),(-SELF_TX_TO_RX_POWER_RATIO));
simulated_channel_model=simulated_channel_model_initial_vector.*AGC_1_Ohm_System(lconv(h,simulated_channel_model_initial_vector,'full'),(-SELF_TX_TO_RX_POWER_RATIO));

%TRAINING LOOP SETUP
tl_itterations=1;

%GRAPHING SET UP
AVERAGING_FILTER_POINTS=2*ceil((NSYMBOLS/USAMPR)/100);
averaging_filter=(1/AVERAGING_FILTER_POINTS).*ones(1,AVERAGING_FILTER_POINTS);

%PREALLOCATE COMPOSITE SIGNAL VECTORS
sout_NP_GDO_EC_COM=[]; s_at_receiver_Post_AGC_COM=[]; s_at_receiver_PCAEC_COM=[];
RX_NOISE_POWER_FROM_EC_COM=[]; RX_NOISE_FROM_EC_dB_COM=[]; SNR_AFTER_EC_COM=[];
ECHO_CANCELATION_COM=[]; MSE_COM=[];
%%%%%%%%%%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% SET UP TRAINING LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tl_n=0;
for tl_n=1:1:tl_itterations
    clearvars -except ALG ...
        pa_coefficients GINIT ...
        NSYMBOLS_LONG_FILTER USAMPR ROLLOFF ORDER SYMBOL_RATE Fc h ...    
        NSYMBOLS N_PILOT_SYMBOLS ...     
        MOD_COD complex_mapping Binary_Alphabet Decimal_Alphabet BITS_PER_WORD PARR_Remapping_Vector Binary_Alphabet ...
        AGC_SET_POWER_TX_dBm ...
        AGC_SET_POWER_RX_dBm RECEIVER_AGC_OVERHEAD_DB ...
        ber_iteration ERROR_LIMIT EbNo ERRORS ...
        NTAPS_OR_PREVEQ cm_itterations ...
        mu ...
        Previous_P epsilon ...
        DATAPATH_TX_TO_RX_POWER_RATIO DATAPATH_TX_TO_RX_VOLTAGE_RATIO SELF_TX_TO_RX_POWER_RATIO SELF_TX_TO_RX_VOLTAGE_RATIO simulated_channel_model_initial_vector simulated_channel_model ...
        tl_itterations ...
        AVERAGING_FILTER_POINTS averaging_filter ...
        sout_NP_GDO_EC_COM s_at_receiver_Post_AGC_COM s_at_receiver_PCAEC_COM RX_NOISE_POWER_FROM_EC_COM RX_NOISE_FROM_EC_dB_COM SNR_AFTER_EC_COM ECHO_CANCELATION_COM MSE_COM
%%%%%%%%%%%%% SET UP TRAINING LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%% GENERATE RANDOM SYMBOLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    srx=randsrc(1,NSYMBOLS,complex_mapping);
    mapped_symbol_stream_receiving=srx;
    stx=randsrc(1,NSYMBOLS,complex_mapping);
    mapped_symbol_stream_transmitting=stx;
    %%%%%%%%%%%%% GENERATE RANDOM SYMBOLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%% ADD PILOT SYMBOLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    PS=randsrc(1,N_PILOT_SYMBOLS,exp(1j*([0,3,1,2].'*pi/2+pi/4)).');
    DATA_GOOD_AT_NP=(length(PS)*USAMPR)+(length(h)-1)/2;
    mapped_symbol_stream_receiving=[PS mapped_symbol_stream_receiving PS];
    mapped_symbol_stream_transmitting=[PS mapped_symbol_stream_transmitting PS];
    %%%%%%%%%%%%% ADD PILOT SYMBOLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%% UPSAMPLE AND FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sup_receiving=upsample(mapped_symbol_stream_receiving,USAMPR);
    sin=lconv(sup_receiving,h,'full');
    sin_NP=sin./max(abs(sin));
    sup_transmitting=upsample(mapped_symbol_stream_transmitting,USAMPR);
    sout=lconv(sup_transmitting,h,'full');
    sout_NP=sout./max(abs(sout));
    %%%%%%%%%%%%% UPSAMPLE AND FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%% TX AGC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [AGC_VOLTAGE_GAIN_OUT sout_NP_GDO]=AGC_1_Ohm_System(sout_NP,AGC_SET_POWER_TX_dBm);
    [TEMP sin_NP_GDO]=AGC_1_Ohm_System(sin_NP,AGC_SET_POWER_TX_dBm); TEMP=[];
    sin_NP_GDO=sin_NP_GDO./DATAPATH_TX_TO_RX_VOLTAGE_RATIO;
    %CHECK IF YOU LIKE
    %One_Ohm_System_Power_dBm(sout_NP_GDO)
    %One_Ohm_System_Power_dBm(sin_NP_GDO)
    %%%%%%%%%%%%% TX AGC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%% CHANNEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[sin_NP_GDO_PAWGN]=AWNG_Generator(sin_NP_GDO,EbNo,USAMPR,h,BITS_PER_WORD);
    [sin_NP_GDO_PAWGN]=sin_NP_GDO;
    sout_NP_GDO_PC=lconv(sout_NP_GDO,simulated_channel_model,'full');
    sout_NP_GDO_PC=sout_NP_GDO_PC((((length(simulated_channel_model)-1)/2)+1):(end-((length(simulated_channel_model)-1)/2)));
    [TEMP sout_NP_GDO_PC]=AGC_1_Ohm_System(sout_NP_GDO_PC,AGC_SET_POWER_TX_dBm); TEMP=[];
    sout_NP_GDO_PC=sout_NP_GDO_PC./SELF_TX_TO_RX_VOLTAGE_RATIO;
    s_at_receiver=sin_NP_GDO_PAWGN+sout_NP_GDO_PC;
    %%%%%%%%%%%%% CHANNEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%% SET UP ITTERATION PER NTAPS_OR_PREVEQ SAMPLES %%%%%%%%%%%%%
    sl_n=0; sl_n_end=[]; STEP_SIZE=[]; START_OFFSET=[]; FILTER_OFFSET=[];
    seqlen=((NSYMBOLS+(2*length(PS)))*USAMPR)+(NSYMBOLS_LONG_FILTER*USAMPR);
    if length(NTAPS_OR_PREVEQ)==1
        sl_n_end=(seqlen-mod(seqlen,NTAPS_OR_PREVEQ));
        STEP_SIZE=NTAPS_OR_PREVEQ;
        FILTER_OFFSET=((NTAPS_OR_PREVEQ-1)/2);
    else
        sl_n_end=(seqlen-mod(seqlen,length(NTAPS_OR_PREVEQ)));
        STEP_SIZE=length(NTAPS_OR_PREVEQ);
        FILTER_OFFSET=((length(NTAPS_OR_PREVEQ)-1)/2);
    end
    START_OFFSET=FILTER_OFFSET+1;
    %%%%%%%%%%%%% SET UP ITTERATION PER NTAPS_OR_PREVEQ SAMPLES %%%%%%%%%%%%%
    
    %%%%%%%%%%%%% ITTERATE PER NTAPS_OR_PREVEQ SAMPLES %%%%%%%%%%%%%%%%%%%%%%
    for sl_n=START_OFFSET:STEP_SIZE:sl_n_end
        vector_values=((sl_n-FILTER_OFFSET):1:(sl_n+FILTER_OFFSET));
        %%%%%%%%%%%%% RX AGC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [s_at_receiver_PAPR]=custom_unitless_PAPR(s_at_receiver(vector_values),[]);
        [COMPOSITE_RX_AGC_V_GAIN s_at_receiver_Post_AGC]=AGC_1_Ohm_System(s_at_receiver(vector_values),AGC_SET_POWER_RX_dBm-s_at_receiver_PAPR-RECEIVER_AGC_OVERHEAD_DB);
        %CHECK IF YOU LIKE
        %One_Ohm_System_Power_dBm(s_at_receiver_Post_AGC)
        %%%%%%%%%%%%% RX AGC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%% FORM CHANNEL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch ALG
            case 'lms'
            %LMS
                [channel_model MSE]=LMS_Channel_Model(s_at_receiver_Post_AGC,sout_NP_GDO(vector_values),NTAPS_OR_PREVEQ,mu,1);
                NTAPS_OR_PREVEQ=channel_model;
            case 'nlms'
            %Normalized LMS
                [channel_model MSE]=NLMS_Channel_Model(s_at_receiver_Post_AGC,sout_NP_GDO(vector_values),NTAPS_OR_PREVEQ,mu,1);
                NTAPS_OR_PREVEQ=channel_model;
            case 'ls'
            %LS
                [channel_model MSE RES]=LS_Channel_Model(s_at_receiver_Post_AGC,sout_NP_GDO(vector_values),NTAPS_OR_PREVEQ);
            case 'rls'
            %RLS
                [channel_model MSE Previous_P]=RLS_Channel_Model(s_at_receiver_Post_AGC,sout_NP_GDO(vector_values),NTAPS_OR_PREVEQ,mu,epsilon,Previous_P,1);
                NTAPS_OR_PREVEQ=channel_model; epsilon=[];
            otherwise
        end
        MSE_COM=[MSE_COM MSE];
        %%%%%%%%%%%%% FORM CHANNEL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%% CANCEL THE ECHO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sout_NP_GDO_EC=lconv(sout_NP_GDO(vector_values),channel_model,'full');
        sout_NP_GDO_EC=sout_NP_GDO_EC((((length(channel_model)-1)/2)+1):(end-((length(channel_model)-1)/2)));
        
        %     %%%%%%%%%%%%% EC AGC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     [TEMP sout_NP_GDO_EC]=AGC_1_Ohm_System(sout_NP_GDO_EC,AGC_SET_POWER_RX_dBm); TEMP=[];
        %     %CHECK IF YOU LIKE
        %     %One_Ohm_System_Power_dBm(sout_NP_GDO_EC)
        %     %%%%%%%%%%%%% EC AGC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        s_at_receiver_PCAEC=s_at_receiver_Post_AGC-sout_NP_GDO_EC;
        
        RX_NOISE_POWER_FROM_EC=s_at_receiver_PCAEC-(sin_NP_GDO_PAWGN(vector_values).*COMPOSITE_RX_AGC_V_GAIN);
        RX_NOISE_FROM_EC_dB=10*log10(((RX_NOISE_POWER_FROM_EC)*(RX_NOISE_POWER_FROM_EC)')/length(RX_NOISE_POWER_FROM_EC));
        SNR_AFTER_EC=One_Ohm_System_Power_dBm(s_at_receiver_PCAEC)-One_Ohm_System_Power_dBm(RX_NOISE_POWER_FROM_EC);
        ECHO_CANCELATION=SNR_AFTER_EC+(DATAPATH_TX_TO_RX_POWER_RATIO-SELF_TX_TO_RX_POWER_RATIO);
        
        %Concatenate results for plotting later
        sout_NP_GDO_EC_COM=[sout_NP_GDO_EC_COM sout_NP_GDO_EC];
        s_at_receiver_Post_AGC_COM=[s_at_receiver_Post_AGC_COM s_at_receiver_Post_AGC];
        s_at_receiver_PCAEC_COM=[s_at_receiver_PCAEC_COM s_at_receiver_PCAEC];
        RX_NOISE_POWER_FROM_EC_COM=[RX_NOISE_POWER_FROM_EC_COM RX_NOISE_POWER_FROM_EC];
        RX_NOISE_FROM_EC_dB_COM=[RX_NOISE_FROM_EC_dB_COM RX_NOISE_FROM_EC_dB];
        SNR_AFTER_EC_COM=[SNR_AFTER_EC_COM SNR_AFTER_EC];
        ECHO_CANCELATION_COM=[ECHO_CANCELATION_COM ECHO_CANCELATION];
    end
    %%%%%%%%%%%%% ITTERATE PER NTAPS_OR_PREVEQ SAMPLES %%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%% CLOSE THE TRAINING LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%% CLOSE THE TRAINING LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TX_SIGNAL_SPECTRUM_TRANSMITTED=lconv(averaging_filter,20*log10(fftshift(abs(fft(sout_NP_GDO,length(sout_NP_GDO))))),'full')-10*log10(0.001);
TX_SIGNAL_POWER_TRANSMITTED=sum(TX_SIGNAL_SPECTRUM_TRANSMITTED)/length(TX_SIGNAL_SPECTRUM_TRANSMITTED);
TX_SIGNAL_SPECTRUM_SELF_RECEIVED=lconv(averaging_filter,20*log10(fftshift(abs(fft(sout_NP_GDO_PC,length(sout_NP_GDO_PC))))),'full')-10*log10(0.001);
TX_SIGNAL_POWER_SELF_RECEIVED=sum(TX_SIGNAL_SPECTRUM_SELF_RECEIVED)/length(TX_SIGNAL_SPECTRUM_SELF_RECEIVED);
TX_SIGNAL_SPECTRUM_DATA_RECEIVED=lconv(averaging_filter,20*log10(fftshift(abs(fft(sin_NP_GDO_PAWGN,length(sin_NP_GDO_PAWGN))))),'full')-10*log10(0.001);
TX_SIGNAL_POWER_DATA_RECEIVED=sum(TX_SIGNAL_SPECTRUM_DATA_RECEIVED)/length(TX_SIGNAL_SPECTRUM_DATA_RECEIVED);
RX_SIGNAL_SPECTRUM_COMPOSITE_POST_AGC_PRE_EC=lconv(averaging_filter,20*log10(fftshift(abs(fft(s_at_receiver_Post_AGC_COM,length(s_at_receiver_Post_AGC_COM))))),'full')-10*log10(0.001);
RX_SIGNAL_POWER_COMPOSITE_POST_AGC_PRE_EC=sum(RX_SIGNAL_SPECTRUM_COMPOSITE_POST_AGC_PRE_EC)/length(RX_SIGNAL_SPECTRUM_COMPOSITE_POST_AGC_PRE_EC);
RX_SIGNAL_SPECTRUM_COMPOSITE_POST_AGC_AND_EC=lconv(averaging_filter,20*log10(fftshift(abs(fft(s_at_receiver_PCAEC_COM,length(s_at_receiver_PCAEC_COM))))),'full')-10*log10(0.001);
RX_SIGNAL_POWER_COMPOSITE_POST_AGC_AND_EC=sum(RX_SIGNAL_SPECTRUM_COMPOSITE_POST_AGC_AND_EC)/length(RX_SIGNAL_SPECTRUM_COMPOSITE_POST_AGC_AND_EC);

RX_SIGNAL_NOISE_SPECTRUM_COMPOSITE_POST_AGC_AND_EC=lconv(averaging_filter,20*log10(fftshift(abs(fft(RX_NOISE_POWER_FROM_EC_COM,length(RX_NOISE_POWER_FROM_EC_COM))))),'full')-10*log10(0.001);
RX_SIGNAL_NOISE_POWER_COMPOSITE_POST_AGC_AND_EC=sum(RX_SIGNAL_NOISE_SPECTRUM_COMPOSITE_POST_AGC_AND_EC)/length(RX_SIGNAL_NOISE_SPECTRUM_COMPOSITE_POST_AGC_AND_EC);

[SNR_AFTER_EC_COM ECHO_CANCELATION_COM]

figure(1), plot(TX_SIGNAL_SPECTRUM_TRANSMITTED), hold on, plot(TX_SIGNAL_SPECTRUM_SELF_RECEIVED,'r'), plot(TX_SIGNAL_SPECTRUM_DATA_RECEIVED,'y'), plot(RX_SIGNAL_SPECTRUM_COMPOSITE_POST_AGC_PRE_EC,'g'), plot(RX_SIGNAL_SPECTRUM_COMPOSITE_POST_AGC_AND_EC,'k'), plot(RX_SIGNAL_NOISE_SPECTRUM_COMPOSITE_POST_AGC_AND_EC,'m'), grid on
axis_1=gca;
axis([0 length(RX_SIGNAL_NOISE_SPECTRUM_COMPOSITE_POST_AGC_AND_EC) -80 70])
legend_1=legend('S1: Self Transmited Signal Spectrum Transmitted','S2: Self Transmited Signal Spectrum Received R-S1','S3: Inbound Transmited Signal Spectrum Received R-S1','S4: Composite Signal Spectrum Received post-AGC pre-EC','S5: Composite Signal Spectrum Received post-AGC post-EC R-S4','S6: Noise Signal Spectrum Due To Echo Cancellaton Relative to R-S5');
title_1=title('Power Spectral Density Of Signal At Various Points');
xlabel_1=xlabel('Digital Frequency Relative To Sampling Frequency');
ylabel_1=ylabel('Power (dBm)');
set(axis_1,'FontSize',12); set(axis_1,'FontName','TimesNewRoman'); 
set(axis_1,'XTick',[0:length(RX_SIGNAL_NOISE_SPECTRUM_COMPOSITE_POST_AGC_AND_EC)/USAMPR:length(RX_SIGNAL_NOISE_SPECTRUM_COMPOSITE_POST_AGC_AND_EC)]);
set(axis_1,'XTickLabel',{'-Fs/2','-Fs*3/8','-Fs/4','-Fs/8','0','Fs/8','Fs/4','Fs*3/8','Fs/2'});
set(legend_1,'Location','south'); set(legend_1,'FontSize',16); set(legend_1,'FontName','TimesNewRoman'); 
set(title_1,'FontSize',20); set(title_1,'FontName','TimesNewRoman');
set(xlabel_1,'FontSize',20); set(xlabel_1,'FontName','TimesNewRoman');
set(ylabel_1,'FontSize',20); set(ylabel_1,'FontName','TimesNewRoman');

figure(2), plot(abs(sout_NP_GDO)), hold on, plot(abs(sout_NP_GDO_PC),'r'), plot(abs(sin_NP_GDO_PAWGN),'y'), plot(abs(s_at_receiver_Post_AGC_COM),'g'), plot(abs(s_at_receiver_PCAEC_COM),'k'), plot(abs(RX_NOISE_POWER_FROM_EC_COM),'m'), grid on
axis_2=gca;
axis([0 length(RX_SIGNAL_NOISE_SPECTRUM_COMPOSITE_POST_AGC_AND_EC) 0 0.2])
legend_2=legend('S1: Self Transmited Signal Transmitted','S2: Self Transmited Signal Received R-S1','S3: Inbound Transmited Signal Received R-S1','S4: Composite Signal Received post-AGC pre-EC','S5: Composite Signal Received post-AGC post-EC R-S4','S6: Noise Signal Due To Echo Cancellaton Relative to R-S5');
title_2=title('Signals Vs Sample Time');
xlabel_2=xlabel('Sample');
ylabel_2=ylabel('Signal (V)');
set(axis_2,'FontSize',12); set(axis_2,'FontName','TimesNewRoman'); 
set(axis_2,'XTick',[0:length(RX_SIGNAL_NOISE_SPECTRUM_COMPOSITE_POST_AGC_AND_EC)/USAMPR:length(RX_SIGNAL_NOISE_SPECTRUM_COMPOSITE_POST_AGC_AND_EC)]);
%set(axis_2,'XTickLabel',{'-Fs/2','-Fs*3/8','-Fs/4','-Fs/8','0','Fs/8','Fs/4','Fs*3/8','Fs/2'});
set(legend_2,'Location','north'); set(legend_2,'FontSize',20); set(legend_2,'FontName','TimesNewRoman');
set(title_2,'FontSize',20); set(title_2,'FontName','TimesNewRoman');
set(xlabel_2,'FontSize',20); set(xlabel_2,'FontName','TimesNewRoman');
set(ylabel_2,'FontSize',20); set(ylabel_2,'FontName','TimesNewRoman');
%%%%%%%%%%%%% CANCEL THE ECHO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% RX FILTER, ALIGN, AND DOWNSAMPLE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%RECEIVE FILTER
expected_bbdata_rx=lconv(s_at_receiver_PCAEC_COM,fliplr(h),'full');

%Align the processed vector
expected_bbdata_rx=circ_shift_2_to_1(sup_receiving,expected_bbdata_rx);

%Chop off ring up and down to pretend it is continuous signal
%Chop needs to be dynamic so it changes with the length of h
expected_bbdata_rx=expected_bbdata_rx(1:1:(end-(length(expected_bbdata_rx)-length(sup_receiving))));
%Chop off ring up and down to pretend it is continuous signal

expected_bbdata_rx=downsample(expected_bbdata_rx./max(abs(expected_bbdata_rx)),USAMPR);

%GetBBdata befor you chop ring off for power measurements
%%%%%%%%%%%% RX FILTER, ALIGN, AND DOWNSAMPLE %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% REMOVE PILOT SYMBOLS AND NORMALIZE %%%%%%%%%%%%%%%%%%%%%%%%%
%Remove pilot symbols
expected_bbdata_rx=expected_bbdata_rx((1+length(PS)):1:(end-length(PS)));

%Normalize
expected_bbdata_rx=expected_bbdata_rx./max(abs(expected_bbdata_rx));
%%%%%%%%%%%% REMOVE PILOT SYMBOLS AND NORMALIZE %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% MAXIMUM LIKLEYHOOD HARD DECISION DECODER %%%%%%%%%%%%%%%%%%%
%Decode
expected_decode_mapping=(mean(abs(expected_bbdata_rx))).*complex_mapping./max(abs(complex_mapping));
[expected_unmapped_symbol_stream]=AWGN_maximum_likelyhood_decoder(expected_bbdata_rx,expected_decode_mapping,complex_mapping);
[Ideal_unmapped_symbol_stream]=AWGN_maximum_likelyhood_decoder(srx,complex_mapping,complex_mapping);
expected_decoded_binary_word_stream=custom_mapper(expected_unmapped_symbol_stream.',complex_mapping,Binary_Alphabet);
Ideal_decoded_binary_word_stream=custom_mapper(Ideal_unmapped_symbol_stream.',complex_mapping,Binary_Alphabet);
%%%%%%%%%%%% MAXIMUM LIKLEYHOOD HARD DECISION DECODER %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% ERROR METRICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate Bit Errors and BER
Bit_Errors=sum(sum(abs(expected_decoded_binary_word_stream-Ideal_decoded_binary_word_stream)));
BER=Bit_Errors/(size(Ideal_decoded_binary_word_stream,1)*size(Ideal_decoded_binary_word_stream,2));
[BER Bit_Errors]

ERRORS=ERRORS+Bit_Errors;
BERCUR=ERRORS/(ber_iteration*NSYMBOLS*BITS_PER_WORD);
[BERCUR ERRORS]

%[SDM SDP SDC Point_Averaged_Frame_Constellation_Soft_Decision]=Constellation_Distortion_Detector(MOD_COD,expected_bbdata_rx,s);
%%%%%%%%%%%% ERROR METRICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%