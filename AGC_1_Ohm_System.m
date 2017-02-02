function [VOLTAGE_GAIN_OUT, gain_controlled_signal]=AGC_1_Ohm_System(signal,DESIRED_POWER_dB)
    
    %FOR A 1 OHM SYSTEM
    %SIGNAL_POWER_dBm=10*log10((sin_NP*sin_NP')/(length(sin_NP)*50*0.001));
    %FOR A 50 OHM SYSTEM
    %SIGNAL_POWER_dBm=10*log10((sin_NP*sin_NP')/(length(sin_NP)*50*0.001));
    %VOLTAGE_GAIN_OUT=power(10,(DESIRED_POWER_dB-SIGNAL_POWER_dBm)/20);
    
    if size(signal,2) < size(signal,1), signal=signal.';, end;
    VOLTAGE_GAIN_OUT=power(10,(DESIRED_POWER_dB-(10*log10((signal*signal')/(size(signal,2)*0.001))))/20);
    gain_controlled_signal=VOLTAGE_GAIN_OUT.*signal;
    
end