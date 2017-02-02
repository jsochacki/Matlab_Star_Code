function [output]=covariance_matrix(input_1,input_2,override)
    if isempty(override), OVERRIDE=0;,else, OVERRIDE=1;, end
    if ((size(input_1,1)==1)|(size(input_1,2)==1))&((size(input_2,1)==1)|(size(input_2,2)==1))
        if size(input_1,1) < size(input_1,2), input_1=input_1.';, end;
        if size(input_2,1) < size(input_2,2), input_2=input_2.';, end;
        if size(input_1,1)==size(input_2,1)
            output=[vector_covariance(input_1,input_1) vector_covariance(input_1,input_2)...
                ;vector_covariance(input_2,input_1) vector_covariance(input_2,input_2)];
        else
            output='INPUT VECTORS MUST BE THE SAME LENGTH';
        end
    elseif ((size(input_1,1)~=size(input_1,2))&~(size(input_2,1)~=size(input_2,2)))|...
            (~(size(input_1,1)~=size(input_1,2))&(size(input_2,1)~=size(input_2,2)))
        output='BOTH INPUTS MUST BE MATRICIES';
    else
        if (size(input_1,1)==size(input_2,1))&(size(input_1,2)==size(input_2,2))
            %FYI THIS HANDLES DIFFERENTLY THAN MATLABS COV FUNCTION AND
            %DOES A(:) B(:) IF A~=B|A==B UNLESS YOU OVERRIDE, THEN IT
            %DOES COVARIANCE(A,B)(NOT IN MATLAB) IF A~=B AND DOES COV(A) IF A==B
            %if (sum(sum((input_1==input_2)))==(size(input_1,1)*size(input_1,2)))&OVERRIDE
            if OVERRIDE
                n=0;
                for n=1:1:size(input_1,2)
                    nn=0;
                    for nn=1:1:size(input_1,2)
                        output(n,nn)=vector_covariance(input_1(:,n),input_2(:,nn));
                    end
                end
            else
                input_1=input_1(:); input_2=input_2(:);
                output=[vector_covariance(input_1,input_1) vector_covariance(input_1,input_2)...
                ;vector_covariance(input_2,input_1) vector_covariance(input_2,input_2)];
            end
        else
            output='BOTH INPUTS MARICIES MUST BE THE SAME SIZE';
        end
    end
end