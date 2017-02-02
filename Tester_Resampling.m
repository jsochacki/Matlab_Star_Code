clear all

%Test setup control
NTAPS_BY_2_LESS_1=12; OTHER_SCALING_dB=110; SELF_SCALING_dB=20; SDOSAMPR=16;
SD_osamper_h_interp=fir1(SDOSAMPR*8,1/SDOSAMPR);

%Algortihm specifics
itterations=1; mu_NLMS=0.1; mu_RLS=0.9; epsilon=1e-6; Previous_P=[];

%Other parameters
USAMPR=4; NSYMBOLS=power(2,13)/SDOSAMPR; alphabet=[exp(j*(pi/4.*[1 3 5 7]))];

%Graphing Set Up
AVERAGING_FILTER_POINTS=2*ceil((NSYMBOLS/USAMPR)/100);
averaging_filter=(1/AVERAGING_FILTER_POINTS).*ones(1,AVERAGING_FILTER_POINTS);

%Generate a random channel model for self path
half_channel=randsrc(1,NTAPS_BY_2_LESS_1,linspace(-10,10,1024))+j.*randsrc(1,NTAPS_BY_2_LESS_1,linspace(-10,10,1024));
channel_model=[half_channel 1 fliplr(half_channel)]; NTAPS_OR_PREVEQ=length(channel_model);

%Generate the tx/rx filter
NSYMBOLS_LONG_FILTER=24; 
ROLLOFF=0.25; ORDER=USAMPR*NSYMBOLS_LONG_FILTER; SYMBOL_RATE=1; Fc=SYMBOL_RATE/2;
h=(firrcos(ORDER,Fc,ROLLOFF,USAMPR,'rolloff','sqrt'));

%Generate the user 2 signal
OTHER_SCALING=power(10,-OTHER_SCALING_dB/20);
other=randsrc(1,NSYMBOLS,alphabet).*OTHER_SCALING;

%Generate the user 1 signal
SELF_SCALING=power(10,-SELF_SCALING_dB/20);
self=randsrc(1,NSYMBOLS,alphabet).*SELF_SCALING;

%Filter for Transmission
other=upsample(other,USAMPR);
self=upsample(self,USAMPR);
other=lconv(other,h,'full');
self=lconv(self,h,'full');

%Create final signals
received_self_full=lconv(self,channel_model,'full');
received_noisefree=other+received_self_full((1+((length(channel_model)-1)/2)):1:(end-((length(channel_model)-1)/2)));
clear received_self_full
received=received_noisefree;

%Filter for Receiption
received=lconv(received,fliplr(h),'full');
self=lconv(self,fliplr(h),'full');
other=lconv(other,fliplr(h),'full');

%Oversample
received_sd=lconv(upsample(received,SDOSAMPR),SD_osamper_h_interp,'full');
self_sd=lconv(upsample(self,SDOSAMPR),SD_osamper_h_interp,'full');
other_sd=lconv(upsample(other,SDOSAMPR),SD_osamper_h_interp,'full');

received_sd=received_sd((1+((length(SD_osamper_h_interp)-1)/2)):1:(end-((length(SD_osamper_h_interp)-1)/2)));
self_sd=self_sd((1+((length(SD_osamper_h_interp)-1)/2)):1:(end-((length(SD_osamper_h_interp)-1)/2)));
other_sd=other_sd((1+((length(SD_osamper_h_interp)-1)/2)):1:(end-((length(SD_osamper_h_interp)-1)/2)));

Units=0; received_il=[]; self_il=[]; other_il=[]; 
CURRENT_NLMS=[]; CURRENT_RLS=[];
ERROR_NLMS=[]; ERROR_RLS=[];
LIMIT=3;
for Units=0:1:LIMIT
    %Add circular shift to immitate non integer sampling clock delays from tx
    %to rx for modified z transform
    desired=[]; current=[];
    received_il(Units+1,:)=downsample(circshift(received_sd,[0 -Units]),SDOSAMPR);
    self_il(Units+1,:)=downsample(self_sd,SDOSAMPR);
    other_il(Units+1,:)=downsample(circshift(other_sd,[0 -Units]),SDOSAMPR);

    received_il(Units+1,:)=((received_il(Units+1,:))./max(abs(received_il(Units+1,:)))).*max(abs(received));
    self_il(Units+1,:)=((self_il(Units+1,:))./max(abs(self_il(Units+1,:)))).*max(abs(self));
    other_il(Units+1,:)=((other_il(Units+1,:))./max(abs(other_il(Units+1,:)))).*max(abs(other));
    
    %Pre algorithm signal selection
    desired=received_il(Units+1,:); current=self_il(Units+1,:);

    NTAPS_OR_PREVEQ=length(channel_model);
    %NLMS
        MSE_NLMS=[]; OCF_Level=1; TAPS_NLMS=[];

        if size(desired,2) < size(desired,1), desired=desired.';, end;
        if size(current,2) < size(current,1), current=current.';, end;

        if sum(size(NTAPS_OR_PREVEQ))==length(size(NTAPS_OR_PREVEQ))
            channel_model_NLMS=ones(1,NTAPS_OR_PREVEQ);
        else
            channel_model_NLMS=NTAPS_OR_PREVEQ;
        end

        if size(channel_model_NLMS,2) < size(channel_model_NLMS,1), channel_model_NLMS=channel_model_NLMS.';, end;
        NTAPS_OR_PREVEQ=size(channel_model_NLMS,2);

        for n=1:1:itterations
            nn=0; nnn=0;
            for nn=(((NTAPS_OR_PREVEQ-1)/2)+1):1:(length(current)-((NTAPS_OR_PREVEQ-1)/2))
                nnn=nnn+1;
                CURRENT=channel_model_NLMS*current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2))).';
                CURRENT_NLMS(Units+1,nnn)=CURRENT;
                DESIRED=desired(nn);
                e=DESIRED-CURRENT;
                if (current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)))*current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)))')>OCF_Level
                    channel_model_NLMS=channel_model_NLMS+2*mu_NLMS.*e.*current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)))'.'/abs(current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)))'.'*current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2))).');
                else
                    channel_model_NLMS=channel_model_NLMS+2*mu_NLMS.*e.*current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)))'.';
                end
                ERROR_NLMS(Units+1,nnn)=e;
                TAPS_NLMS=[TAPS_NLMS channel_model_NLMS.'];
            end
            MSE_NLMS=[MSE_NLMS 10*log10((ERROR_NLMS(Units+1,:)*ERROR_NLMS(Units+1,:)')/((length(current)-((NTAPS_OR_PREVEQ-1)/2))-(((NTAPS_OR_PREVEQ-1)/2))))];
            %NORMALIZED MSE
            %MSE_NLMS=[MSE_NLMS 10*log10(((ERROR_NLMS*ERROR_NLMS')/max(ERROR_NLMS*ERROR_NLMS'))/((length(current)-((NTAPS_OR_PREVEQ-1)/2))-(((NTAPS_OR_PREVEQ-1)/2))))];
        end
        %figure(1)
        NLMS_PLOT_VECTOR=(((NTAPS_OR_PREVEQ-1)/2)+1):1:(length(current)-((NTAPS_OR_PREVEQ-1)/2));
        NLMS_onestep=lconv(current,channel_model_NLMS,'full');
        NLMS_onestep_out(Units+1,:)=NLMS_onestep((((length(channel_model_NLMS)-1)/2)+1):(end-((length(channel_model_NLMS)-1)/2)));
        MSE_NLMS_out(Units+1)=MSE_NLMS;
%         plot(abs(CURRENT_NLMS),'b') %Per itteration signal
%         hold on , plot(abs(ERROR_NLMS),'r') %Per itteration error
%         plot(abs(desired(NLMS_PLOT_VECTOR)),'k') %Received Signal
%         plot(abs(NLMS_onestep(NLMS_PLOT_VECTOR)),'m') %One step signal
%         plot(abs(desired(NLMS_PLOT_VECTOR)-NLMS_onestep(NLMS_PLOT_VECTOR)),'g') %One step error
%         plot(abs(other(NLMS_PLOT_VECTOR)),'y') %Desired signal from user 2
%         plot(abs(other(NLMS_PLOT_VECTOR)-ERROR_NLMS),'c') %per itteration noise
%         plot(abs(other(NLMS_PLOT_VECTOR)-(desired(NLMS_PLOT_VECTOR)-NLMS_onestep(NLMS_PLOT_VECTOR))),'.c') %one step noise
%         hold off, grid on
    %NLMS

    NTAPS_OR_PREVEQ=length(channel_model);
    %RLS
        MSE_RLS=[]; TAPS_RLS=[];
        e=[]; CURRENT=[]; DESIRED=[]; current_vec=[]; gain=[]; PI_subo=[]; P=[];

        %mu is the forgetting factor here and should be 0 << mu <=1
        mu_inv=1/mu_RLS;

        %PI_subo is the weighting matrix which is initialized to
        %P=PI_subo=I*epsilon where epsilon is a very small positive number
        %(this makes PI_subo^-1 large which is good for an initial estimate
        %(yeilds lms for first itteration roughly))

        if size(desired,2) < size(desired,1), desired=desired.';, end;
        if size(current,2) < size(current,1), current=current.';, end;

        if sum(size(NTAPS_OR_PREVEQ))==length(size(NTAPS_OR_PREVEQ))
            channel_model_RLS=ones(1,NTAPS_OR_PREVEQ);
            P=eye(NTAPS_OR_PREVEQ)*epsilon; Previous_P=[];
        else
            channel_model_RLS=NTAPS_OR_PREVEQ;
            P=Previous_P; epsilon=[];
        end

        if size(channel_model_RLS,2) < size(channel_model_RLS,1), channel_model_RLS=channel_model_RLS.';, end;
        NTAPS_OR_PREVEQ=size(channel_model_RLS,2);

        for n=1:1:itterations
            nn=0; nnn=0;
            for nn=(((NTAPS_OR_PREVEQ-1)/2)+1):1:(length(current)-((NTAPS_OR_PREVEQ-1)/2))
                nnn=nnn+1;
    %STANDARD
    %             current_vec=current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)));
    %             CURRENT=linear_equalizer*current_vec.';
    %             DESIRED=desired(nn);
    %             e=DESIRED-CURRENT;
    %             gain=(mu_inv*P*current_vec.')/(1+mu_inv*(current_vec*(P*current_vec.')));
    %             linear_equalizer=linear_equalizer+gain.'*e;
    %             P=mu_inv*(P-gain*(current_vec*P));
    %             ERROR=[ERROR e];
    %COMPLEX
                current_vec=current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)));
                CURRENT=channel_model_RLS*current_vec.';
                CURRENT_RLS(Units+1,nnn)=CURRENT;
                DESIRED=desired(nn);
                e=DESIRED-CURRENT;
                gain=(mu_inv*P*current_vec.')/(1+mu_inv*(current_vec'.'*(P*current_vec.')));
                channel_model_RLS=channel_model_RLS+gain'*e;
                P=mu_inv*(P-gain*(current_vec'.'*P));
                ERROR_RLS(Units+1,nnn)=e;
                TAPS_RLS=[TAPS_RLS channel_model_RLS.'];
            end
            MSE_RLS=[MSE_RLS 10*log10((ERROR_RLS(Units+1,:)*ERROR_RLS(Units+1,:)')/((length(current)-((NTAPS_OR_PREVEQ-1)/2))-(((NTAPS_OR_PREVEQ-1)/2))))];
            %NORMALIZED MSE
            %MSE_RLS=[MSE_RLS 10*log10(((ERROR_RLS*ERROR_RLS')/max(ERROR_RLS*ERROR_RLS'))/((length(current)-((NTAPS_OR_PREVEQ-1)/2))-(((NTAPS_OR_PREVEQ-1)/2))))];
        end
        %figure(2)
        RLS_PLOT_VECTOR=(((NTAPS_OR_PREVEQ-1)/2)+1):1:(length(current)-((NTAPS_OR_PREVEQ-1)/2));
        RLS_onestep=lconv(current,channel_model_RLS,'full');
        RLS_onestep_out(Units+1,:)=RLS_onestep((((length(channel_model_RLS)-1)/2)+1):(end-((length(channel_model_RLS)-1)/2)));
        MSE_RLS_out(Units+1)=MSE_RLS;
%         plot(abs(CURRENT_RLS),'b') %Per itteration signal
%         hold on , plot(abs(ERROR_RLS),'r') %Per itteration error
%         plot(abs(desired(RLS_PLOT_VECTOR)),'k') %Received Signal
%         plot(abs(RLS_onestep(RLS_PLOT_VECTOR)),'m') %One step signal
%         plot(abs(desired(RLS_PLOT_VECTOR)-NLMS_onestep(RLS_PLOT_VECTOR)),'g') %One step error
%         plot(abs(other(RLS_PLOT_VECTOR)),'y') %Desired signal from user 2
%         plot(abs(other(RLS_PLOT_VECTOR)-ERROR_RLS),'c') %per itteration noise
%         plot(abs(other(RLS_PLOT_VECTOR)-(desired(RLS_PLOT_VECTOR)-RLS_onestep(RLS_PLOT_VECTOR))),'.c') %one step noise
%         hold off, grid on
    %RLS

    NTAPS_OR_PREVEQ=length(channel_model);
    %LS
    [channel_model_LS MSE_LS RES]=LS_Channel_Model(desired,current,NTAPS_OR_PREVEQ);
    %figure(3)
    LS_onestep=lconv(current,channel_model_LS,'full');
    LS_onestep_out(Units+1,:)=LS_onestep((((length(channel_model_LS)-1)/2)+1):(end-((length(channel_model_LS)-1)/2)));
%     plot(abs(desired(RLS_PLOT_VECTOR)),'k') %Received Signal
%     hold on
%     plot(abs(LS_onestep(RLS_PLOT_VECTOR)),'m') %One step signal
%     plot(abs(desired(RLS_PLOT_VECTOR)-LS_onestep(RLS_PLOT_VECTOR)),'g') %One step error
%     plot(abs(other(RLS_PLOT_VECTOR)),'y') %Desired signal from user 2
%     plot(abs(other(RLS_PLOT_VECTOR)-(desired(RLS_PLOT_VECTOR)-LS_onestep(RLS_PLOT_VECTOR))),'.c') %one step noise
%     hold off, grid on
    %LS
end

desired_vector=repmat(desired((((NTAPS_OR_PREVEQ-1)/2)+1):1:(end-((NTAPS_OR_PREVEQ-1)/2))),1,itterations);

% for n=1:1:LIMIT+1
%     figure(n)
%     semilogy(abs(ERROR_NLMS(n,:)),'r')
%     hold on
%     semilogy(abs(ERROR_RLS(n,:)),'g')
%     semilogy(abs(received_il(n,:)-LS_onestep_out(n,:)),'k')
%     semilogy(abs(other_il(n,:)),'y')
%     grid on
%     hold off
%     axis([0 length(other) 1e-10 1])
% end
 
for n=1:1:LIMIT+1
    figure(n)
    semilogy(abs(other(NLMS_PLOT_VECTOR)-ERROR_NLMS(n,:)),'r')
    hold on
    semilogy(abs(other(RLS_PLOT_VECTOR)-ERROR_RLS(n,:)),'g')
    semilogy(abs(other(RLS_PLOT_VECTOR)-(received_il(n,RLS_PLOT_VECTOR)-LS_onestep_out(n,RLS_PLOT_VECTOR))),'k')
    semilogy(abs(other),'y')
    grid on
    hold off
    axis([0 length(other) 1e-10 1])
end

% figure(5)
% hold on
% for n=1:1:size(TAPS_NLMS,1)
%     plot(abs(TAPS_NLMS(n,:)),'r')
% end
% hold on
% for n=1:1:size(TAPS_RLS,1)
%     plot(abs(TAPS_RLS(n,:)),'g')
% end
% grid on
% hold off
% 
% desired_vector=repmat(desired((((NTAPS_OR_PREVEQ-1)/2)+1):1:(end-((NTAPS_OR_PREVEQ-1)/2))),1,itterations);
% 
% figure(9)
% plot(abs((desired_vector-CURRENT_NLMS(1,:))),'b')
% hold on
% plot(abs((desired_vector-CURRENT_RLS(1,:))),'k')
% grid on
% hold off
% 
% figure(6)
% plot(abs(desired_vector-CURRENT_NLMS(1,:)),'r')
% hold on
% plot(abs(desired_vector-CURRENT_RLS(1,:)),'k')
% grid on
% hold off
% 
% MSE_NLMS=10*log10(((CURRENT_NLMS(1,:)-desired_noisefree_vector)*(CURRENT_NLMS(1,:)-desired_noisefree_vector)')/length(desired_noisefree_vector));
% MSE_RLS=10*log10(((CURRENT_RLS(1,:)-desired_noisefree_vector)*(CURRENT_RLS(1,:)-desired_noisefree_vector)')/length(desired_noisefree_vector));
% 
% figure(7)
% plot(abs(CURRENT_NLMS(1,:)),'r')
% hold on
% plot(abs(CURRENT_RLS(1,:)),'g')
% hold on
% plot(abs(desired_vector),'y')
% grid on
% hold off
% 
% CURRENT_NLMS(1,:)=CURRENT_NLMS(1,:)((end-length(CURRENT_NLMS(1,:))/5):end);
% CURRENT_RLS(1,:)=CURRENT_RLS(1,:)((end-length(CURRENT_RLS(1,:))/5):end);
% desired_noisefree_vector=desired_noisefree_vector((end-length(desired_noisefree_vector)/5):end);
% CURRENT_NLMS_SPECTRUM=lconv(averaging_filter,20*log10(fftshift(abs(fft(CURRENT_NLMS(1,:),length(CURRENT_NLMS(1,:)))))),'full')-10*log10(0.001);
% CURRENT_NLMS_SPECTRUM_ERROR=lconv(averaging_filter,20*log10(fftshift(abs(fft(CURRENT_NLMS(1,:)-desired_noisefree_vector,length(CURRENT_NLMS(1,:)-desired_noisefree_vector))))),'full')-10*log10(0.001);
% CURRENT_RLS_SPECTRUM=lconv(averaging_filter,20*log10(fftshift(abs(fft(CURRENT_RLS(1,:),length(CURRENT_RLS(1,:)))))),'full')-10*log10(0.001);
% CURRENT_RLS_SPECTRUM_ERROR=lconv(averaging_filter,20*log10(fftshift(abs(fft(CURRENT_RLS(1,:)-desired_noisefree_vector,length(CURRENT_RLS(1,:)-desired_noisefree_vector))))),'full')-10*log10(0.001);
% 
% figure(8)
% plot((CURRENT_NLMS_SPECTRUM),'k')
% hold on
% plot((CURRENT_NLMS_SPECTRUM_ERROR),'r')
% hold on
% plot((CURRENT_RLS_SPECTRUM),'b')
% hold on
% plot((CURRENT_RLS_SPECTRUM_ERROR),'g')
% grid on
% hold off
% 
% figure(1)
% plot(channel_model_LMS,'b')
% hold on
% plot(channel_model_NLMS,'r')
% hold on
% plot(channel_model_LS,'k')
% hold on
% plot(channel_model_RLS,'g')
% hold on
% plot(channel_model,'y')
% grid on
% 
% channel_model_SPECTRUM=10*log10(abs(fftshift(fft(channel_model,length(channel_model)))));
% channel_model_LMS_SPECTRUM=10*log10(abs(fftshift(fft(channel_model_LMS,length(channel_model_LMS)))));
% channel_model_NLMS_SPECTRUM=10*log10(abs(fftshift(fft(channel_model_NLMS,length(channel_model_NLMS)))));
% channel_model_LS_SPECTRUM=10*log10(abs(fftshift(fft(channel_model_LS,length(channel_model_LS)))));
% channel_model_RLS_SPECTRUM=10*log10(abs(fftshift(fft(channel_model_RLS,length(channel_model_RLS)))));
% 
% figure(2)
% plot(channel_model_LMS_SPECTRUM,'b')
% hold on
% plot(channel_model_NLMS_SPECTRUM,'r')
% hold on
% plot(channel_model_LS_SPECTRUM,'k')
% hold on
% plot(channel_model_RLS_SPECTRUM,'g')
% hold on
% plot(channel_model_SPECTRUM,'y')
% grid on
% 
% fvtool(channel_model_NLMS)
% fvtool(channel_model_LS)
% fvtool(channel_model_RLS)
% fvtool(channel_model)
% 
% figure(3)
% plot(abs(lconv(current,channel_model_LMS,'full')-desired_noisefree_full),'b')
% hold on
% plot(abs(lconv(current,channel_model_NLMS,'full')-desired_noisefree_full),'r')
% hold on
% plot(abs(lconv(current,channel_model_LS,'full')-desired_noisefree_full),'k')
% hold on
% plot(abs(lconv(current,channel_model_RLS,'full')-desired_noisefree_full),'g')
% grid on
