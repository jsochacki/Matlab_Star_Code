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
ORDER_F=ORDER/(USAMPR/4);
h_f=fir1(ORDER_F,4/USAMPR,hamming(1+ORDER_F));
%%%%%%%%%%%%% FILTER SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% h_f=zeros(1,54);
% h_f(1)=1;
% h_f(2)=-1.279827863;
% h_f(8)=0.3489001603;
% h_f(22)=-0.090299256;
% h_f(39)=0.0260998396;
% h_f(54)=-0.004872880;

% h_f=h_f./2.57;


% h_f=zeros(1,8);
% h_f(1)=1;
% h_f(2)=-1.700333;
% h_f(4)=0.653712;
% h_f(5)=0.196287;
% h_f(8)=-0.149667;
% 
% h_f=h_f./3.7;

h_f=fir1(3,0.5);

%h_f=zeros(1,54);

% figure(3)
% plot(10*log10(abs(fftshift(fft(x,length(x))))))
% hold on
% plot(10*log10(abs(fftshift(fft(filtfilt(h_f,1,x),length(filtfilt(h_f,1,x)))))),'r')
% hold off

%IN PATH LOW PASS FILTER
% e_c=0; x_c=0; e=[]; xvec=[]; x_q=[]; x_in=x; 
% for n=1:1:length(x)
%     x_in(n)=x(n)-e_c;
%     if n < length(h_f)
%         xvec=[x_in(1:1:n) zeros(1,length(h_f)-n)];
%     else
%         xvec=x_in(n-((length(h_f)-1)):1:n);
%     end
%     x_c=xvec*h_f.';
%     xq(n)=Sample_Quantizer(x_c,max(abs(x)),4,'mid-tread');
%     e(n)=xq(n)-x_c; e_c=e(n);
% end

en_m1=0; xq=[]; e=[]; evec=[]; x_in=x;
for n=1:1:length(x)
    x_in(n)=x(n)-en_m1;
    xq(n)=Sample_Quantizer(x_in(n),max(abs(x)),4,'mid-tread');
    e(n)=xq(n)-x(n);
    if n < length(h_f)
        evec=[e(1:1:n) zeros(1,length(h_f)-length(e))];
    else
        evec=e(n-((length(h_f)-1)):1:n);
    end
    en_m1=evec*h_f.';
end

%hold on, plot(10*log10(fftshift(power(abs(fft(xq,length(xq))),2))),'k')

% figure(1)
% plot(abs(x),'b')
% hold on
% %plot(abs(xq(((length(h_f)-1)/2)+1:1:end)),'r')
% plot(abs(xq),'r')
% hold off

figure(2)
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