function [output]=Custom_1_Bit_DAC(input,VTH_P_IN,VTH_N_IN,VTH_P_OUT,VTH_N_OUT)
    output=((input==VTH_P_IN).*VTH_P_OUT)+((input==VTH_N_IN).*VTH_N_OUT);
end