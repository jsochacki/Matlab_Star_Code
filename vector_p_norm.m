function [output]=vector_p_norm(vector,p)
    if p=='cheb'
        %p='cheb' yeilds the max of the absolute values 
        %Also known as the chebyshev norm
        output=max(abs(vector));
    elseif p==1
        %p=1 yeilds the sum of the absolute values
        output=sum(abs(vector));
    else
        %Minskowski norm
        %when p=2 this is commonly known as the euclidian norm or
        %euclidian distance/length
        output=power(sum(power(abs(vector),p)),1/p);
    end
end