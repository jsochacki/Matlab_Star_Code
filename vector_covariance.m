function [output]=vector_covariance(input_1,input_2)
    if size(input_1,1) < size(input_1,2), input_1=input_1.';, end;
    if size(input_2,1) < size(input_2,2), input_2=input_2.';, end;
    if size(input_1,1)==size(input_2,1)
        output=(1/(size(input_1,1)-1))*...
            ((input_1-((ones(1,size(input_1,1))*input_1)/size(input_1,1)))'*...
            (input_2-((ones(1,size(input_2,1))*input_2)/size(input_2,1))));
    else
        output='INPUT VECTORS MUST BE THE SAME LENGTH';
    end
    
end