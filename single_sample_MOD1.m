function [output Feedback_Output Integrator_Output]=single_sample_MOD1(Signal_Input,Feedback_Input,Old_Integrator_Input,COMPARITOR_THRESHOLD,VH_COMPARITOR,VL_COMPARITOR,VH_DAC,VL_DAC)
    Integrator_Input=Signal_Input-Feedback_Input;
    Integrator_Output=Custom_Integrator(Old_Integrator_Input,Integrator_Input);
    output=Custom_1_Bit_Comparator(Integrator_Output,COMPARITOR_THRESHOLD,VH_COMPARITOR,VL_COMPARITOR);
    Feedback_Output=Custom_1_Bit_DAC(output,VH_COMPARITOR,VL_COMPARITOR,VH_DAC,VL_DAC);
end