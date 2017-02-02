simulated_channel_model=[linspace(0,1,13)+j*linspace(1,0,13) linspace(1,0,12)+j*linspace(0,1,12)];
srand=randsrc(1,1000,[exp(j*[pi/4 3*pi/4 5*pi/4 7*pi/4])]);
desired=conv(simulated_channel_model,srand);
current=srand;
NTAPS_OR_PREVEQ=25; mu=0.01; itterations=100;
[channel_model MSE]=LMS_Channel_Model(desired,current,NTAPS_OR_PREVEQ,mu,itterations);
plot(desired);
hold on;
plot(conv(channel_model,current),'r')