clear all
%%%%%%%%%%%%% DATA SIGNAL SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NSYMBOLS=128;

%SET UP FOR DATA FILTER
USAMPR=128; NSYMBOLS_LONG_FILTER=24; 
ROLLOFF=0.25; ORDER=USAMPR*NSYMBOLS_LONG_FILTER; SYMBOL_RATE=1; Fc=SYMBOL_RATE/2;
h_c=firrcos(ORDER,Fc,ROLLOFF,USAMPR,'rolloff','sqrt');

%GENERATE DATA SIGNAL
complex_mapping=[exp(-j*pi/4.*[1 3 5 7])];
data=randsrc(1,NSYMBOLS,complex_mapping);
data_2=randsrc(1,NSYMBOLS,complex_mapping);
data_up=upsample(data,USAMPR);
data_2_up=upsample(data_2,USAMPR);
x=lconv(data_up,h_c,'full');
x_2=lconv(data_2_up,h_c,'full');
%%%%%%%%%%%%% DATA SIGNAL SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% FILTER SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAMPR_F=8; NSYMBOLS_LONG_FILTER_F=4; 
% ROLLOFF_F=0.25; ORDER_F=USAMPR_F*NSYMBOLS_LONG_FILTER_F; SYMBOL_RATE_F=1; Fc_F=SYMBOL_RATE_F/2;
% h_f=firrcos(ORDER_F,Fc_F,ROLLOFF_F,USAMPR_F,'rolloff','sqrt');
% ORDER_F=ORDER/(USAMPR/4);
% h_f=fir1(ORDER_F,4/USAMPR,hamming(1+ORDER_F));
%%%%%%%%%%%%% FILTER SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%x_ec=x.*exp(-j*0*(2*pi/360)); % 0 degrees results in a difference spectrum
%that is perfectly overlaid (no change) while 180 degrees results in a
%difference spectrum that is 3dB higher than either separate spectrum when
%x is identical to x_ec
%%%
%x_ec=circshift(x,[0 -16]); %circular shifts( or simply difference in time
%of arrival) result in very large changes in the receive power spectrum
%that seem to bee cyclical according to the amount of delay but also have a
%asymptotic approach (like a sin wave times decaying exp) but the most
%notable and important feature is that the noise transfer function becomes
%cyclical and only overlays the correct noise transfer function when there
%is no shift between the two signals and has becomes very aliased as the
%amount of delay increases
%%%
%x_ec=0.3*x; %Changes in the amplitude have zero effect on the difference
%signal if the output of the MOD1s have the same voltage limits (1 and zero
%or vH and Vl) and the feedback dacs also have the same voltage levels
%(vfbH and vfbL) to the input of the MOD1.  They do have a difference if
%each MOD is connected to it's own set of limits.  THis is seen in the form
%of the STF going up and down in magnitude cyclically as the ratio of the
%voltage is altered in a direction while the NTF remains almost exactly the
%same until the point where the signal level starts to be too far from the
%FS range of the MOD and the SQNR is too large and then the NTF starts to
%look a little weird but still not that much.
x_ec=x;

x_2=x_2*(power(10,-30/20));
x=x+x_2; %Add in user 2 signal
EPBOR=mean(abs(x_2)/mean(abs(x-x_2)));
%NOTE: there is a big difference in oversampling (ADC) a nyquist converter
%waveform that has been oversampled on tx ADC (say by 4x) and a SD converter
%with the same OSRTx as the nyquist will sit a vout for the OSR samples
% (where OSR = OSRRx/OSRTx) the SD converter output will not and will be
% interpolated between the points

% %FIRST ORDER ADC
% x_r_q=0; x_r_in=real(x); x_r_Feedback_Output=0; x_r_Integrator_Output=0;
% x_i_q=0; x_i_in=imag(x); x_i_Feedback_Output=0; x_i_Integrator_Output=0;
% xq=[];
% x_ec_r_q=0; x_ec_r_in=real(x_ec); x_ec_r_Feedback_Output=0; x_ec_r_Integrator_Output=0;
% x_ec_i_q=0; x_ec_i_in=imag(x_ec); x_ec_i_Feedback_Output=0; x_ec_i_Integrator_Output=0;
% x_ec_q=[];
% diff_r=0; diff_i=0; diff_c=[];
% for n=1:1:length(x)
%     [x_r_q x_r_Feedback_Output x_r_Integrator_Output]=single_sample_MOD1(x_r_in(n),x_r_Feedback_Output,x_r_Integrator_Output,mean(x_r_in),1,0,max(x_r_in),-max(x_r_in));
%     [x_i_q x_i_Feedback_Output x_i_Integrator_Output]=single_sample_MOD1(x_i_in(n),x_i_Feedback_Output,x_i_Integrator_Output,mean(x_i_in),1,0,max(x_i_in),-max(x_i_in));
%     [x_ec_r_q x_ec_r_Feedback_Output x_ec_r_Integrator_Output]=single_sample_MOD1(x_ec_r_in(n),x_ec_r_Feedback_Output,x_ec_r_Integrator_Output,mean(x_ec_r_in),1,0,max(x_ec_r_in),-max(x_ec_r_in));
%     [x_ec_i_q x_ec_i_Feedback_Output x_ec_i_Integrator_Output]=single_sample_MOD1(x_ec_i_in(n),x_ec_i_Feedback_Output,x_ec_i_Integrator_Output,mean(x_ec_i_in),1,0,max(x_ec_i_in),-max(x_ec_i_in));
%     %[x_ec_r_q x_ec_r_Feedback_Output x_ec_r_Integrator_Output]=single_sample_MOD1(x_ec_r_in(n),x_ec_r_Feedback_Output,x_ec_r_Integrator_Output,mean(x_r_in),1,0,max(x_r_in),-max(x_r_in));
%     %[x_ec_i_q x_ec_i_Feedback_Output x_ec_i_Integrator_Output]=single_sample_MOD1(x_ec_i_in(n),x_ec_i_Feedback_Output,x_ec_i_Integrator_Output,mean(x_i_in),1,0,max(x_i_in),-max(x_i_in));
%     xq(n)=x_r_q+j.*x_i_q;
%     x_ec_q(n)=x_ec_r_q+j.*x_ec_i_q;
%     diff_r=x_r_q-x_ec_r_q; diff_i=x_i_q-x_ec_i_q;
%     diff_c(n)=diff_r+j.*diff_i;
% end

% %FIRST ORDER ADC
% x_r_q=0; x_r_in=real(x); x_r_Feedback_Output=0; x_r_Integrator_Output=0;
% x_i_q=0; x_i_in=imag(x); x_i_Feedback_Output=0; x_i_Integrator_Output=0;
% xq=[];
% x_ec_r_q=0; x_ec_r_in=real(x_ec); x_ec_r_Feedback_Output=0; x_ec_r_Integrator_Output=0;
% x_ec_i_q=0; x_ec_i_in=imag(x_ec); x_ec_i_Feedback_Output=0; x_ec_i_Integrator_Output=0;
% x_ec_q=[];
% diff_r=0; diff_i=0; diff_c=[];
% for n=1:1:length(x)
%     [x_ec_r_q x_ec_r_Feedback_Output x_ec_r_Integrator_Output]=single_sample_MOD1(x_ec_r_in(n),x_ec_r_Feedback_Output,x_ec_r_Integrator_Output,mean(x_ec_r_in),1,0,max(x_ec_r_in),-max(x_ec_r_in));
%     [x_ec_i_q x_ec_i_Feedback_Output x_ec_i_Integrator_Output]=single_sample_MOD1(x_ec_i_in(n),x_ec_i_Feedback_Output,x_ec_i_Integrator_Output,mean(x_ec_i_in),1,0,max(x_ec_i_in),-max(x_ec_i_in));
%     %[x_ec_r_q x_ec_r_Feedback_Output x_ec_r_Integrator_Output]=single_sample_MOD1(x_ec_r_in(n),x_ec_r_Feedback_Output,x_ec_r_Integrator_Output,mean(x_r_in),1,0,max(x_r_in),-max(x_r_in));
%     %[x_ec_i_q x_ec_i_Feedback_Output x_ec_i_Integrator_Output]=single_sample_MOD1(x_ec_i_in(n),x_ec_i_Feedback_Output,x_ec_i_Integrator_Output,mean(x_i_in),1,0,max(x_i_in),-max(x_i_in));
%     [x_r_q x_r_Feedback_Output x_r_Integrator_Output]=single_sample_MOD1(x_r_in(n)-x_ec_r_Feedback_Output,x_r_Feedback_Output,x_r_Integrator_Output,mean(x_r_in),1,0,max(x_r_in),-max(x_r_in));
%     [x_i_q x_i_Feedback_Output x_i_Integrator_Output]=single_sample_MOD1(x_i_in(n)-x_ec_i_Feedback_Output,x_i_Feedback_Output,x_i_Integrator_Output,mean(x_i_in),1,0,max(x_i_in),-max(x_i_in));
%     xq(n)=x_r_q+j.*x_i_q;
%     x_ec_q(n)=x_ec_r_q+j.*x_ec_i_q;
%     diff_r=x_r_q-x_ec_r_q; diff_i=x_i_q-x_ec_i_q;
%     diff_c(n)=diff_r+j.*diff_i;
% end

%FIRST ORDER ADC
x_r_q=0; x_r_in=real(x); x_r_Feedback_Output=0; x_r_Integrator_Output=0;
x_i_q=0; x_i_in=imag(x); x_i_Feedback_Output=0; x_i_Integrator_Output=0;
xq=[];
x_ec_r_q=0; x_ec_r_in=real(x_ec); x_ec_r_Feedback_Output=0; x_ec_r_Integrator_Output=0;
x_ec_i_q=0; x_ec_i_in=imag(x_ec); x_ec_i_Feedback_Output=0; x_ec_i_Integrator_Output=0;
x_ec_q=[];
diff_r=0; diff_i=0; diff_c=[];
for n=1:1:length(x)
    [x_ec_r_q x_ec_r_Feedback_Output x_ec_r_Integrator_Output]=single_sample_MOD1(x_ec_r_in(n),x_ec_r_Feedback_Output,x_ec_r_Integrator_Output,mean(x_ec_r_in),1,0,max(x_ec_r_in),-max(x_ec_r_in));
    [x_ec_i_q x_ec_i_Feedback_Output x_ec_i_Integrator_Output]=single_sample_MOD1(x_ec_i_in(n),x_ec_i_Feedback_Output,x_ec_i_Integrator_Output,mean(x_ec_i_in),1,0,max(x_ec_i_in),-max(x_ec_i_in));
    %[x_ec_r_q x_ec_r_Feedback_Output x_ec_r_Integrator_Output]=single_sample_MOD1(x_ec_r_in(n),x_ec_r_Feedback_Output,x_ec_r_Integrator_Output,mean(x_r_in),1,0,max(x_r_in),-max(x_r_in));
    %[x_ec_i_q x_ec_i_Feedback_Output x_ec_i_Integrator_Output]=single_sample_MOD1(x_ec_i_in(n),x_ec_i_Feedback_Output,x_ec_i_Integrator_Output,mean(x_i_in),1,0,max(x_i_in),-max(x_i_in));
    [x_r_q x_r_Feedback_Output x_r_Integrator_Output]=single_sample_MOD1(x_r_in(n)-x_ec_r_Feedback_Output,x_r_Feedback_Output,x_r_Integrator_Output,mean(x_r_in),1,0,max(x_r_in),-max(x_r_in));
    [x_i_q x_i_Feedback_Output x_i_Integrator_Output]=single_sample_MOD1(x_i_in(n)-x_ec_i_Feedback_Output,x_i_Feedback_Output,x_i_Integrator_Output,mean(x_i_in),1,0,max(x_i_in),-max(x_i_in));
    xq(n)=x_r_q+j.*x_i_q;
    x_ec_q(n)=x_ec_r_q+j.*x_ec_i_q;
    diff_r=x_r_q-x_ec_r_q; diff_i=x_i_q-x_ec_i_q;
    diff_c(n)=diff_r+j.*diff_i;
end
x_p=x(1+((length(h_c)-1)/2):1:(end-((length(h_c)-1)/2)));
xq_p=xq(1+((length(h_c)-1)/2):1:(end-((length(h_c)-1)/2)));
x_ec_q_p=x_ec_q(1+((length(h_c)-1)/2):1:(end-((length(h_c)-1)/2)));
diff_c_p=diff_c(1+((length(h_c)-1)/2):1:(end-((length(h_c)-1)/2)));
x_2_p=x_2(1+((length(h_c)-1)/2):1:(end-((length(h_c)-1)/2)));

figure(1)
plot(10*log10(fftshift(power(abs(fft(x_p,length(x_p))),2))),'b')
hold on
plot(10*log10(fftshift(power(abs(fft(xq_p,length(xq_p))),2))),'r')
plot(10*log10(fftshift(power(abs(fft(x_ec_q_p,length(x_ec_q_p))),2))),'g')
plot(10*log10(fftshift(power(abs(fft(diff_c_p,length(diff_c_p))),2))),'k')
plot(10*log10(fftshift(power(abs(fft(x_2_p,length(x_2_p))),2))),'c')
plot(10*log10(fftshift(power(abs(fft(filtfilt(fir1(USAMPR*8,2/USAMPR),1,xq_p),length(xq_p))),2))),'r')
hold off

figure(2)
plot(abs(x_p-x_2_p),'b')
hold on
%plot(abs(xq(((length(h_f)-1)/2)+1:1:end)),'r')
plot(abs(xq_p),'r')
plot(abs(x_ec_q_p),'g')
plot(abs(diff_c_p),'k')
plot(abs(x_2_p),'c')
hold off

DFXQ=filtfilt(fir1(USAMPR*8,2/USAMPR),1,xq);
DFC=filtfilt(fir1(USAMPR*8,2/USAMPR),1,diff_c);
DFXQ=DFXQ(1+((length(h_c)-1)/2):1:(end-((length(h_c)-1)/2)));
DFC=DFC(1+((length(h_c)-1)/2):1:(end-((length(h_c)-1)/2)));
DFXQ=(DFXQ-mean(DFXQ));
DFXQ=(DFXQ./max(abs(DFXQ))).*max(abs(x_2_p));
figure(3)
plot(abs(x_p-x_2_p),'b')
hold on
%plot(abs(xq(((length(h_f)-1)/2)+1:1:end)),'r')
plot(abs(DFXQ./max(abs(DFXQ(floor(length(DFXQ))/4:floor(length(DFXQ))*3/4)))),'r')
plot(abs(filtfilt(fir1(USAMPR*8,2/USAMPR),1,x_ec_q_p)),'g')
plot(abs(DFC./max(abs(DFC(1:floor(length(DFC))/2)))),'k')
plot(abs(x_2_p./max(abs(x_2_p))),'c')
hold off

expected_decode_mapping_DFXQ=(mean(abs(DFXQ))).*complex_mapping./max(abs(complex_mapping));
[expected_unmapped_symbol_stream_DFXQ]=AWGN_maximum_likelyhood_decoder(downsample(DFXQ,USAMPR),expected_decode_mapping_DFXQ,complex_mapping);
expected_decode_mapping_x_2=(mean(abs(x_2_p))).*complex_mapping./max(abs(complex_mapping));
[expected_unmapped_symbol_stream_x_2]=AWGN_maximum_likelyhood_decoder(downsample(x_2_p,USAMPR),expected_decode_mapping_x_2,complex_mapping);

figure(4)
plot(expected_unmapped_symbol_stream_DFXQ,'bo')
hold on
plot(expected_unmapped_symbol_stream_x_2,'ro')
plot(expected_unmapped_symbol_stream_DFXQ-expected_unmapped_symbol_stream_x_2,'ko')

% figure(3)
% plot(10*log10(fftshift(power(abs(fft(x,length(x))),2))),'b')
% hold on
% plot(10*log10(fftshift(power(abs(fft(xq,length(xq))),2))),'r')
% hold off

% figure(4)
% plot(abs(filtfilt(fir1(128,2/USAMPR),1,x)),'b')
% hold on
% %plot(abs(xq(((length(h_f)-1)/2)+1:1:end)),'r')
% plot(abs(filtfilt(fir1(128,2/USAMPR),1,xq)),'r')
% hold off
% 
% figure(5)
% plot(10*log10(fftshift(power(abs(fft(filtfilt(fir1(128,2/USAMPR),1,x),length(filtfilt(fir1(128,2/USAMPR),1,x)))),2))),'b')
% hold on
% plot(10*log10(fftshift(power(abs(fft(filtfilt(fir1(128,2/USAMPR),1,xq),length(filtfilt(fir1(128,2/USAMPR),1,xq)))),2))),'r')
% hold off