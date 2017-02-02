clear all
%%%%%%%%%%%%% DATA SIGNAL SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NSYMBOLS=128;

%SET UP FOR DATA FILTER
USAMPR=128; NSYMBOLS_LONG_FILTER=24; 
ROLLOFF=0.25; ORDER=USAMPR*NSYMBOLS_LONG_FILTER; SYMBOL_RATE=1; Fc=SYMBOL_RATE/2;
h_c=firrcos(ORDER,Fc,ROLLOFF,USAMPR,'rolloff','sqrt');

%GENERATE DATA SIGNAL
data=randsrc(1,NSYMBOLS,[exp(-j*pi/4.*[1 3 5 7])]);
data_up=upsample(data,USAMPR);
x=lconv(data_up,h_c,'full');
%%%%%%%%%%%%% DATA SIGNAL SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% FILTER SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAMPR_F=8; NSYMBOLS_LONG_FILTER_F=4; 
% ROLLOFF_F=0.25; ORDER_F=USAMPR_F*NSYMBOLS_LONG_FILTER_F; SYMBOL_RATE_F=1; Fc_F=SYMBOL_RATE_F/2;
% h_f=firrcos(ORDER_F,Fc_F,ROLLOFF_F,USAMPR_F,'rolloff','sqrt');
% ORDER_F=ORDER/(USAMPR/4);
% h_f=fir1(ORDER_F,4/USAMPR,hamming(1+ORDER_F));
%%%%%%%%%%%%% FILTER SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(3)
% plot(10*log10(abs(fftshift(fft(x,length(x))))))
% hold on
% plot(10*log10(abs(fftshift(fft(filtfilt(h_f,1,x),length(filtfilt(h_f,1,x)))))),'r')
% hold off

%FIRST ORDER ADC
x_r_q=0; x_r_in=real(x); r_Feedback_Output=0; r_Integrator_Output=0;
x_i_q=0; x_i_in=imag(x); i_Feedback_Output=0; i_Integrator_Output=0;
xq=[];
for n=1:1:length(x)
    [x_r_q r_Feedback_Output r_Integrator_Output]=single_sample_MOD1(x_r_in(n),r_Feedback_Output,r_Integrator_Output,mean(x_r_in),1,0,max(x_r_in),-max(x_r_in));
    [x_i_q i_Feedback_Output i_Integrator_Output]=single_sample_MOD1(x_i_in(n),i_Feedback_Output,i_Integrator_Output,mean(x_i_in),1,0,max(x_i_in),-max(x_i_in));
    xq(n)=x_r_q+j.*x_i_q;
end

%FIRST ORDER ADC
w_r_n1=0; x_r_s1_n1=0; x_r_q=[]; x_r_s1=[]; x_r_i1=[]; x_r_in=real(x);
w_i_n1=0; x_i_s1_n1=0; x_i_q=[]; x_i_s1=[]; x_i_i1=[]; x_i_in=imag(x);
xq=[];
for n=1:1:length(x)
    x_r_s1(n)=x_r_in(n)-w_r_n1;
    x_i_s1(n)=x_i_in(n)-w_i_n1;
    x_r_i1(n)=Custom_Integrator(x_r_s1(n),x_r_s1_n1); x_r_s1_n1=x_r_i1(n);
    x_i_i1(n)=Custom_Integrator(x_i_s1(n),x_i_s1_n1); x_i_s1_n1=x_i_i1(n);
    x_r_q(n)=Custom_1_Bit_Comparator(x_r_i1(n),mean(x_r_in),1,0);
    x_i_q(n)=Custom_1_Bit_Comparator(x_i_i1(n),mean(x_i_in),1,0);
    w_r_n1=Custom_1_Bit_DAC(x_r_q(n),1,0,max(x_r_in),-max(x_r_in));
    w_i_n1=Custom_1_Bit_DAC(x_i_q(n),1,0,max(x_i_in),-max(x_i_in));
    xq(n)=x_r_q(n)+j.*x_i_q(n);
end

%FIRST ORDER DAC
w_r_n1=0; x_r_q=[]; x_r_s1=[]; x_r_in=real(x);
w_i_n1=0; x_i_q=[]; x_i_s1=[]; x_i_in=imag(x);
xq=[];
for n=1:1:length(x)
    x_r_s1(n)=Custom_Integrator(x_r_in(n),w_r_n1);
    x_i_s1(n)=Custom_Integrator(x_i_in(n),w_i_n1);
    x_r_q(n)=Custom_1_Bit_Comparator(x_r_s1(n),mean(x_r_in),max(abs(x_r_in)),-max(abs(x_r_in)));
    x_i_q(n)=Custom_1_Bit_Comparator(x_i_s1(n),mean(x_i_in),max(abs(x_i_in)),-max(abs(x_i_in)));
    w_r_n1=x_r_s1(n)-x_r_q(n);
    w_i_n1=x_i_s1(n)-x_i_q(n);
    xq(n)=x_r_q(n)+j.*x_i_q(n);
end

figure(2)
plot(10*log10(fftshift(power(abs(fft(x,length(x))),2))),'b')
hold on
plot(10*log10(fftshift(power(abs(fft(xq,length(xq))),2))),'r')

%SECOND ORDER
%This is programmed properly but it is unstable for inputs that are close
%to full scale so will appear not to work (as it has gone unstable and is
%oscillating) so you must not run the feedback full scale i.e. modify the
%w_(r/i)_n1 value by making it larger then Max(abs(x)) (the velue that comes
%out of the 1 bid dac
w_r_n1=0; x_r_s1_n1=0; x_r_s2_n1=0; x_r_q=[]; x_r_s1=[]; x_r_s2=[]; x_r_i1=[]; x_r_i2=[]; x_r_in=real(x);
w_i_n1=0; x_i_s1_n1=0; x_i_s2_n1=0; x_i_q=[]; x_i_s1=[]; x_i_s2=[]; x_i_i1=[]; x_i_i2=[]; x_i_in=imag(x);
xq=[];
for n=1:1:length(x)
    x_r_s1(n)=x_r_in(n)-w_r_n1;
    x_i_s1(n)=x_i_in(n)-w_i_n1;
    x_r_i1(n)=Custom_Integrator(x_r_s1(n),x_r_s1_n1);
    x_i_i1(n)=Custom_Integrator(x_i_s1(n),x_i_s1_n1);
    
    x_r_s1_n1=x_r_i1(n);
    x_i_s1_n1=x_i_i1(n);
    x_r_s2(n)=x_r_s1_n1-w_r_n1;
    x_i_s2(n)=x_i_s1_n1-w_i_n1;

    x_r_i2(n)=Custom_Integrator(x_r_s2(n),x_r_s2_n1); x_r_s2_n1=x_r_i2(n);
    x_i_i2(n)=Custom_Integrator(x_i_s2(n),x_i_s2_n1); x_i_s2_n1=x_i_i2(n);
    x_r_q(n)=Custom_1_Bit_Comparator(x_r_i2(n),0,1,0);
    x_i_q(n)=Custom_1_Bit_Comparator(x_i_i2(n),0,1,0);
    w_r_n1=Custom_1_Bit_DAC(x_r_q(n),1,0,max(x_r_in),-max(x_r_in));
    w_i_n1=Custom_1_Bit_DAC(x_i_q(n),1,0,max(x_i_in),-max(x_i_in));
    xq(n)=x_r_q(n)+j.*x_i_q(n);
end

%SECOND ORDER DAC
w_r_n1=0; x_r_q=[]; x_r_q1=[]; x_r_s1=[]; x_r_s2=[]; x_r_in=real(x);
w_i_n1=0; x_i_q=[]; x_i_q2=[]; x_i_s1=[]; x_i_s2=[]; x_i_in=imag(x);
xq=[];
%w_r_n2=0; w_i_n2=0;
for n=1:1:length(x)
    x_r_s1(n)=Custom_Integrator(x_r_in(n),w_r_n1);
    x_i_s1(n)=Custom_Integrator(x_i_in(n),w_i_n1);
    x_r_q(n)=Custom_1_Bit_Comparator(x_r_s1(n),mean(x_r_in),max(abs(x_r_in)),-max(abs(x_r_in)));
    x_i_q(n)=Custom_1_Bit_Comparator(x_i_s1(n),mean(x_i_in),max(abs(x_i_in)),-max(abs(x_i_in)));
    %w_r_n1=x_r_s1(n)-x_r_q(n);
    %w_i_n1=x_i_s1(n)-x_i_q(n);
    x_r_s2(n)=Custom_Integrator(x_r_q(n),w_r_n1);
    x_i_s2(n)=Custom_Integrator(x_i_q(n),w_i_n1);
    %x_r_s2(n)=Custom_Integrator(x_r_q(n),w_r_n2);
    %x_i_s2(n)=Custom_Integrator(x_i_q(n),w_i_n2);
    x_r_q1(n)=Custom_1_Bit_Comparator(x_r_s2(n),mean(x_r_in),max(abs(x_r_in)),-max(abs(x_r_in)));
    x_i_q1(n)=Custom_1_Bit_Comparator(x_i_s2(n),mean(x_i_in),max(abs(x_i_in)),-max(abs(x_i_in)));
    w_r_n1=x_r_s2(n)-x_r_q1(n);
    w_i_n1=x_i_s2(n)-x_i_q1(n);
    %w_r_n2=x_r_s2(n)-x_r_q1(n);
    %w_i_n2=x_i_s2(n)-x_i_q1(n);
    xq(n)=x_r_q1(n)+j.*x_i_q1(n);
end

plot(10*log10(fftshift(power(abs(fft(xq,length(xq))),2))),'k')
hold off
%hold on, plot(10*log10(fftshift(power(abs(fft(xq,length(xq))),2))),'k')

% figure(1)
% plot(abs(x),'b')
% hold on
% %plot(abs(xq(((length(h_f)-1)/2)+1:1:end)),'r')
% plot(abs(xq),'r')
% hold off

figure(6)
plot(10*log10(fftshift(power(abs(fft(x,length(x))),2))),'b')
hold on
plot(10*log10(fftshift(power(abs(fft(xq,length(xq))),2))),'r')
hold off

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