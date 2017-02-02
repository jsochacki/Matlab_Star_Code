function [output]=Custom_1_Bit_Comparitor(input,VTH,VTH_P,VTH_N)
    output=((input>=VTH).*VTH_P)+((input<VTH).*VTH_N);
end