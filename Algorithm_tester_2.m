clear all

%SUMMARY:
%ALL ARE PERFECT WHEN RUNNING AT THE SAME RATE AND NUMBER OF TAPS
%(MSE=-250)
%DIFFERING NUMBER OF TAPS AND SAME RATE YEILDS all MSE=-260 so long as taps
%are greater then the number of taps required to model the channel
%i.e. (length(channel_model)
%DIFFERING RATE AND SAME NUMBER OF TAPS YEILDS all MSE=-260 so long as taps
%are greater then the number of taps required to model the channel
%i.e. (length(channel_model)
%DIFFERING RATE AND NUMBER OF TAPS YEILDS all MSE=-260 so long as taps
%are greater then the number of taps required to model the channel
%i.e. (length(channel_model)
%SAME RATE AND NUMBER OF TAPS WITH NOISE YEILDS

NTAPS_BY_2_LESS_1=12; NOISE_ON=1; FILTER_ON=1;
half_channel=randsrc(1,NTAPS_BY_2_LESS_1,linspace(-10,10,1024))+j.*randsrc(1,NTAPS_BY_2_LESS_1,linspace(-10,10,1024));
channel_model=[half_channel 1 fliplr(half_channel)]; NTAPS_OR_PREVEQ=length(channel_model);
Previous_P=[]; USAMPR=4;
NSYMBOLS=512; alphabet=[exp(j*(pi/4.*[1 3 5 7]))];

itterations=5; mu_NLMS=0.1; mu_RLS=0.9; epsilon=1e-6; Previous_P=[];

NSYMBOLS_LONG_FILTER=24; 
ROLLOFF=0.25; ORDER=USAMPR*NSYMBOLS_LONG_FILTER; SYMBOL_RATE=1; Fc=SYMBOL_RATE/2;
h=(firrcos(ORDER,Fc,ROLLOFF,USAMPR,'rolloff','sqrt'));

noise=NOISE_ON.*randsrc(1,NSYMBOLS+(NSYMBOLS_LONG_FILTER*FILTER_ON),alphabet); NOISE_SCALING_dB=60;
NOISE_SCALING=power(10,-NOISE_SCALING_dB/20);
noise=upsample(noise,USAMPR).*NOISE_SCALING;
current=randsrc(1,NSYMBOLS,alphabet);
current=upsample(current,USAMPR);
if FILTER_ON
    current=lconv(current,h,'full');
end
desired_noisefree_full=lconv(current,channel_model,'full');
desired_noisefree=desired_noisefree_full((1+((length(channel_model)-1)/2)):1:(end-((length(channel_model)-1)/2)));
desired=lconv(current+noise,channel_model,'full');
desired=desired((1+((length(channel_model)-1)/2)):1:(end-((length(channel_model)-1)/2)));

%GRAPHING SET UP
AVERAGING_FILTER_POINTS=2*ceil((NSYMBOLS/USAMPR)/100);
averaging_filter=(1/AVERAGING_FILTER_POINTS).*ones(1,AVERAGING_FILTER_POINTS);

%NLMS
    MSE_NLMS=[]; ERROR_NLMS=[]; OCF_Level=1; TAPS_NLMS=[]; CURRENT_NLMS=[];
    
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
        nn=0;
        for nn=(((NTAPS_OR_PREVEQ-1)/2)+1):1:(length(current)-((NTAPS_OR_PREVEQ-1)/2))
            CURRENT=channel_model_NLMS*current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2))).';
            CURRENT_NLMS=[CURRENT_NLMS CURRENT];
            DESIRED=desired(nn);
            e=DESIRED-CURRENT;
            if (current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)))*current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)))')>OCF_Level
                channel_model_NLMS=channel_model_NLMS+2*mu_NLMS.*e.*current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)))'.'/abs(current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)))'.'*current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2))).');
            else
                channel_model_NLMS=channel_model_NLMS+2*mu_NLMS.*e.*current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)))'.';
            end
            ERROR_NLMS=[ERROR_NLMS e];
            TAPS_NLMS=[TAPS_NLMS channel_model_NLMS.'];
        end
        MSE_NLMS=[MSE_NLMS 10*log10((ERROR_NLMS*ERROR_NLMS')/((length(current)-((NTAPS_OR_PREVEQ-1)/2))-(((NTAPS_OR_PREVEQ-1)/2))))];
        %NORMALIZED MSE
        %MSE_NLMS=[MSE_NLMS 10*log10(((ERROR_NLMS*ERROR_NLMS')/max(ERROR_NLMS*ERROR_NLMS'))/((length(current)-((NTAPS_OR_PREVEQ-1)/2))-(((NTAPS_OR_PREVEQ-1)/2))))];
    end
%NLMS

%RLS
    MSE_RLS=[]; ERROR_RLS=[]; TAPS_RLS=[]; CURRENT_RLS=[];
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
        nn=0;
        for nn=(((NTAPS_OR_PREVEQ-1)/2)+1):1:(length(current)-((NTAPS_OR_PREVEQ-1)/2))
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
            CURRENT_RLS=[CURRENT_RLS CURRENT];
            DESIRED=desired(nn);
            e=DESIRED-CURRENT;
            gain=(mu_inv*P*current_vec.')/(1+mu_inv*(current_vec'.'*(P*current_vec.')));
            channel_model_RLS=channel_model_RLS+gain'*e;
            P=mu_inv*(P-gain*(current_vec'.'*P));
            ERROR_RLS=[ERROR_RLS e];
            TAPS_RLS=[TAPS_RLS channel_model_RLS.'];
        end
        MSE_RLS=[MSE_RLS 10*log10((ERROR_RLS*ERROR_RLS')/((length(current)-((NTAPS_OR_PREVEQ-1)/2))-(((NTAPS_OR_PREVEQ-1)/2))))];
        %NORMALIZED MSE
        %MSE_RLS=[MSE_RLS 10*log10(((ERROR_RLS*ERROR_RLS')/max(ERROR_RLS*ERROR_RLS'))/((length(current)-((NTAPS_OR_PREVEQ-1)/2))-(((NTAPS_OR_PREVEQ-1)/2))))];
    end
%RLS

%LS
[channel_model_LS MSE_LS RES]=LS_Channel_Model(desired,current,NTAPS_OR_PREVEQ);
%LS

figure(4)
semilogy(abs(ERROR_NLMS),'r')
hold on
semilogy(abs(ERROR_RLS),'g')
grid on
hold off

figure(5)
hold on
for n=1:1:size(TAPS_NLMS,1)
    plot(abs(TAPS_NLMS(n,:)),'r')
end
hold on
for n=1:1:size(TAPS_RLS,1)
    plot(abs(TAPS_RLS(n,:)),'g')
end
grid on
hold off

desired_noisefree_vector=repmat(desired_noisefree((((NTAPS_OR_PREVEQ-1)/2)+1):1:(end-((NTAPS_OR_PREVEQ-1)/2))),1,itterations);
desired_vector=repmat(desired((((NTAPS_OR_PREVEQ-1)/2)+1):1:(end-((NTAPS_OR_PREVEQ-1)/2))),1,itterations);

figure(9)
plot(abs((desired_vector-CURRENT_NLMS)-(desired_vector-desired_noisefree_vector)),'b')
hold on
plot(abs((desired_vector-CURRENT_RLS)-(desired_vector-desired_noisefree_vector)),'k')
grid on
hold off

figure(6)
plot(abs(desired_noisefree_vector-CURRENT_NLMS),'r')
hold on
plot(abs(desired_vector-CURRENT_NLMS),'b')
hold on
plot(abs(desired_noisefree_vector-CURRENT_RLS),'g')
hold on
plot(abs(desired_vector-CURRENT_RLS),'k')
hold on
plot(abs(desired_vector-desired_noisefree_vector),'m') %this is just the noise
grid on
hold off

MSE_NLMS=10*log10(((CURRENT_NLMS-desired_noisefree_vector)*(CURRENT_NLMS-desired_noisefree_vector)')/length(desired_noisefree_vector));
MSE_RLS=10*log10(((CURRENT_RLS-desired_noisefree_vector)*(CURRENT_RLS-desired_noisefree_vector)')/length(desired_noisefree_vector));

figure(7)
plot(abs(CURRENT_NLMS),'r')
hold on
plot(abs(CURRENT_RLS),'g')
hold on
plot(abs(desired_noisefree_vector),'y')
grid on
hold off

CURRENT_NLMS=CURRENT_NLMS((end-length(CURRENT_NLMS)/5):end);
CURRENT_RLS=CURRENT_RLS((end-length(CURRENT_RLS)/5):end);
desired_noisefree_vector=desired_noisefree_vector((end-length(desired_noisefree_vector)/5):end);
CURRENT_NLMS_SPECTRUM=lconv(averaging_filter,20*log10(fftshift(abs(fft(CURRENT_NLMS,length(CURRENT_NLMS))))),'full')-10*log10(0.001);
CURRENT_NLMS_SPECTRUM_ERROR=lconv(averaging_filter,20*log10(fftshift(abs(fft(CURRENT_NLMS-desired_noisefree_vector,length(CURRENT_NLMS-desired_noisefree_vector))))),'full')-10*log10(0.001);
CURRENT_RLS_SPECTRUM=lconv(averaging_filter,20*log10(fftshift(abs(fft(CURRENT_RLS,length(CURRENT_RLS))))),'full')-10*log10(0.001);
CURRENT_RLS_SPECTRUM_ERROR=lconv(averaging_filter,20*log10(fftshift(abs(fft(CURRENT_RLS-desired_noisefree_vector,length(CURRENT_RLS-desired_noisefree_vector))))),'full')-10*log10(0.001);

figure(8)
plot((CURRENT_NLMS_SPECTRUM),'k')
hold on
plot((CURRENT_NLMS_SPECTRUM_ERROR),'r')
hold on
plot((CURRENT_RLS_SPECTRUM),'b')
hold on
plot((CURRENT_RLS_SPECTRUM_ERROR),'g')
grid on
hold off

figure(1)
plot(channel_model_LMS,'b')
hold on
plot(channel_model_NLMS,'r')
hold on
plot(channel_model_LS,'k')
hold on
plot(channel_model_RLS,'g')
hold on
plot(channel_model,'y')
grid on

channel_model_SPECTRUM=10*log10(abs(fftshift(fft(channel_model,length(channel_model)))));
channel_model_LMS_SPECTRUM=10*log10(abs(fftshift(fft(channel_model_LMS,length(channel_model_LMS)))));
channel_model_NLMS_SPECTRUM=10*log10(abs(fftshift(fft(channel_model_NLMS,length(channel_model_NLMS)))));
channel_model_LS_SPECTRUM=10*log10(abs(fftshift(fft(channel_model_LS,length(channel_model_LS)))));
channel_model_RLS_SPECTRUM=10*log10(abs(fftshift(fft(channel_model_RLS,length(channel_model_RLS)))));

figure(2)
plot(channel_model_LMS_SPECTRUM,'b')
hold on
plot(channel_model_NLMS_SPECTRUM,'r')
hold on
plot(channel_model_LS_SPECTRUM,'k')
hold on
plot(channel_model_RLS_SPECTRUM,'g')
hold on
plot(channel_model_SPECTRUM,'y')
grid on

fvtool(channel_model_NLMS)
fvtool(channel_model_LS)
fvtool(channel_model_RLS)
fvtool(channel_model)

figure(3)
plot(abs(lconv(current,channel_model_LMS,'full')-desired_noisefree_full),'b')
hold on
plot(abs(lconv(current,channel_model_NLMS,'full')-desired_noisefree_full),'r')
hold on
plot(abs(lconv(current,channel_model_LS,'full')-desired_noisefree_full),'k')
hold on
plot(abs(lconv(current,channel_model_RLS,'full')-desired_noisefree_full),'g')
grid on
