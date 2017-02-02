function [quantized_output]=Sample_Quantizer(input,MAX_AMPLITUDE,WORD_LENGTH,ROUNDING_METHOD)
%input_range should be [minimum_value maximum_value]
    delta=MAX_AMPLITUDE/(power(2,WORD_LENGTH)-1);
    switch ROUNDING_METHOD
        case 'mid-tread'
            quantized_output=delta*floor((input/delta)+(1/2));
    end
end