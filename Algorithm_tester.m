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
half_channel=randsrc(1,12,linspace(-10,10,1024))+j.*randsrc(1,12,linspace(-10,10,1024));
channel_model=[half_channel 1 fliplr(half_channel)]; NTAPS_OR_PREVEQ=length(channel_model);
Previous_P=[]; USAMPR=1;
NSYMBOLS=1024; alphabet=[exp(j*(pi/4.*[1 3 5 7]))];
noise=randsrc(1,NSYMBOLS,alphabet); NOISE_SCALING=0.1;
noise=upsample(noise,USAMPR).*NOISE_SCALING;
current=randsrc(1,NSYMBOLS,alphabet);
current=upsample(current,USAMPR);
desired_noisefree=lconv(current,channel_model,'full');
desired=lconv(current+noise,channel_model,'full');
desired=desired((1+((length(channel_model)-1)/2)):1:(end-((length(channel_model)-1)/2)));

[channel_model_LMS MSE_LMS]=LMS_Channel_Model(desired,current,NTAPS_OR_PREVEQ,0.1,5);
[channel_model_NLMS MSE_NLMS]=NLMS_Channel_Model(desired,current,NTAPS_OR_PREVEQ,0.1,5);
[channel_model_LS MSE_LS RES]=LS_Channel_Model(desired,current,NTAPS_OR_PREVEQ);
[channel_model_RLS MSE_RLS P_out]=RLS_Channel_Model(desired,current,NTAPS_OR_PREVEQ,0.9,1e-6,Previous_P,5);

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
plot(abs(lconv(current,channel_model_LMS,'full')-desired_noisefree),'b')
hold on
plot(abs(lconv(current,channel_model_NLMS,'full')-desired_noisefree),'r')
hold on
plot(abs(lconv(current,channel_model_LS,'full')-desired_noisefree),'k')
hold on
plot(abs(lconv(current,channel_model_RLS,'full')-desired_noisefree),'g')
grid on
