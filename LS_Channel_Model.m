function [channel_model MSE RES]=LS_Channel_Model(desired,current,NTAPS)
    
    channel_model=[]; MSE=[]; RES=[];
    channel_model_1=[]; channel_model_2=[]; channel_1=[]; channel_2=[];
    residual_1=[]; residual_2=[]; e_1=[]; e_2=[]; MSE_1=[]; MSE_2=[];
    
    if size(desired,2) < size(desired,1), desired=desired.';, end;
    if size(current,2) < size(current,1), current=current.';, end;
    
    %Create the Input Matrix A
    %NOTE: Only works with odd number of taps as I don't ever plan to run
    %even
    A=[]; n=0;
    for n=(1+((NTAPS-1)/2)):1:(length(current)-((NTAPS-1)/2))
        A=[A;current((n-((NTAPS-1)/2)):(n+((NTAPS-1)/2)))];
    end
    
    %Create the Sample Matrix B
    B=desired((1+((NTAPS-1)/2)):1:(length(current)-((NTAPS-1)/2)));
    
    %Perform the Left Division
    %This uses matlabs built in QR decomposition solver
    %In the end you need to implement the QR decomposition on your own
    channel_model_1=(A\B.').';
    
    %Perform the Left Moore-Penrose Matrix Pseudoinverse
    A_plus=[];
    [A_plus]=left_mp_matrix_pseudoinverse(A);
    channel_model_2=(A_plus*B.').';
    
    %Find the error, norm, and residuals of the solutions
    channel_1=lconv(current,channel_model_1,'full');
    channel_1=channel_1((((length(channel_model_1)-1)/2)+1):(end-((length(channel_model_1)-1)/2)));
    e_1=channel_1-desired;
    residual_1=vector_p_norm(e_1,2);
    %MSE_1=10*log10((e_1*e_1')/((length(current)-((NTAPS-1)/2))-(1+((NTAPS-1)/2))));
    %NORMALIZED MSE
    MSE_1=10*log10(((e_1*e_1')/max(e_1*e_1'))/((length(current)-((NTAPS-1)/2))-(1+((NTAPS-1)/2))));
    
    channel_2=lconv(current,channel_model_2,'full');
    channel_2=channel_2((((length(channel_model_2)-1)/2)+1):(end-((length(channel_model_2)-1)/2)));
    e_2=channel_2-desired;
    residual_2=vector_p_norm(channel_2-desired,2);
    %MSE_2=10*log10((e_2*e_2')/((length(current)-((NTAPS-1)/2))-(1+((NTAPS-1)/2))));
    %NORMALIZED MSE
    MSE_2=10*log10(((e_2*e_2')/max(e_2*e_2'))/((length(current)-((NTAPS-1)/2))-(1+((NTAPS-1)/2))));
    
    if residual_1 < residual_2
        channel_model=channel_model_1; MSE=MSE_1; RES=residual_1;
    else
        channel_model=channel_model_2; MSE=MSE_2; RES=residual_2;
    end

end