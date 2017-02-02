function [ONE_OHM_SYSTEM_POWER_dBm]=One_Ohm_System_Power_dBm(signal)
    
    %FOR A 1 OHM SYSTEM
    %SIGNAL_POWER_dBm=10*log10((sin_NP*sin_NP')/(length(sin_NP)*50*0.001));
    %FOR A 50 OHM SYSTEM
    %SIGNAL_POWER_dBm=10*log10((sin_NP*sin_NP')/(length(sin_NP)*50*0.001));
    %VOLTAGE_GAIN_OUT=power(10,(DESIRED_POWER_dB-SIGNAL_POWER_dBm)/20);
    
    if size(signal,2) < size(signal,1), signal=signal.';, end;
    ONE_OHM_SYSTEM_POWER_dBm=10*log10((signal*signal')/(size(signal,2)*0.001));

end