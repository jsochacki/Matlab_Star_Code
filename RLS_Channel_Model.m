function [channel_model MSE P_out]=RLS_Channel_Model(desired,current,NTAPS_OR_PREVEQ,mu,epsilon,Previous_P,itterations)
    
    MSE=[]; ERROR=[];
    e=[]; CURRENT=[]; DESIRED=[]; current_vec=[]; gain=[]; PI_subo=[]; P=[];

    %mu is the forgetting factor here and should be 0 << mu <=1
    mu_inv=1/mu;
    
    %PI_subo is the weighting matrix which is initialized to
    %P=PI_subo=I*epsilon where epsilon is a very small positive number
    %(this makes PI_subo^-1 large which is good for an initial estimate
    %(yeilds lms for first itteration roughly))
    
    if size(desired,2) < size(desired,1), desired=desired.';, end;
    if size(current,2) < size(current,1), current=current.';, end;
    
    if sum(size(NTAPS_OR_PREVEQ))==length(size(NTAPS_OR_PREVEQ))
        linear_equalizer=ones(1,NTAPS_OR_PREVEQ);
        P=eye(NTAPS_OR_PREVEQ)*epsilon; Previous_P=[];
    else
        linear_equalizer=NTAPS_OR_PREVEQ;
        P=Previous_P; epsilon=[];
    end

    if size(linear_equalizer,2) < size(linear_equalizer,1), linear_equalizer=linear_equalizer.';, end;
    NTAPS_OR_PREVEQ=size(linear_equalizer,2);

    for n=1:1:itterations
        nn=0;
        for nn=(((NTAPS_OR_PREVEQ-1)/2)+1):1:(length(current)-((NTAPS_OR_PREVEQ-1)/2))
%STANDARD
%             current_vec=current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)));
%             CURRENT=linear_equalizer*current_vec.';
%             DESIRED=desired(nn);
%             e=DESIRED-CURRENT;
%             gain=(mu_inv*P*current_vec.')/(1+mu_inv*(current_vec*(P*current_vec.')));
%             linear_equalizer=linear_equalizer+gain.'*e;
%             P=mu_inv*(P-gain*(current_vec*P));
%             ERROR=[ERROR e];
%COMPLEX
            current_vec=current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)));
            CURRENT=linear_equalizer*current_vec.';
            DESIRED=desired(nn);
            e=DESIRED-CURRENT;
            gain=(mu_inv*P*current_vec.')/(1+mu_inv*(current_vec'.'*(P*current_vec.')));
            linear_equalizer=linear_equalizer+gain'*e;
            P=mu_inv*(P-gain*(current_vec'.'*P));
            ERROR=[ERROR e];
        end
        MSE=[MSE 10*log10((ERROR*ERROR')/((length(current)-((NTAPS_OR_PREVEQ-1)/2))-(((NTAPS_OR_PREVEQ-1)/2))))]; ERROR=[];
        %NORMALIZED MSE
        %MSE=[MSE 10*log10(((ERROR*ERROR')/max(ERROR*ERROR'))/((length(current)-((NTAPS_OR_PREVEQ-1)/2))-(((NTAPS_OR_PREVEQ-1)/2))))]; ERROR=[];
    end
    channel_model=linear_equalizer; P_out=P;
end