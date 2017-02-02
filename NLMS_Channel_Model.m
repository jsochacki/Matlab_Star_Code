function [channel_model MSE CURRENT]=LMS_Channel_Model(desired,current,NTAPS_OR_PREVEQ,mu,itterations)
    %Note: 0<mu<=1 as i have the 2 in the update equation already as
    %opposed the the usual 0<mu<=2
    
    %This uses the orthogonal correction factor method modified by me to
    %pulloff to the standard LMS when the signal power is less than or
    %equal to 1 and the NLMS algorithm for signal powers greater than 1
    
    MSE=[]; ERROR=[]; OCF_Level=1;
    
    if size(desired,2) < size(desired,1), desired=desired.';, end;
    if size(current,2) < size(current,1), current=current.';, end;
    
    if sum(size(NTAPS_OR_PREVEQ))==length(size(NTAPS_OR_PREVEQ))
        linear_equalizer=ones(1,NTAPS_OR_PREVEQ);
    else
        linear_equalizer=NTAPS_OR_PREVEQ;
    end
    
    if size(linear_equalizer,2) < size(linear_equalizer,1), linear_equalizer=linear_equalizer.';, end;
    NTAPS_OR_PREVEQ=size(linear_equalizer,2);
    
    for n=1:1:itterations
        nn=0;
        for nn=(((NTAPS_OR_PREVEQ-1)/2)+1):1:(length(current)-((NTAPS_OR_PREVEQ-1)/2))
            CURRENT=linear_equalizer*current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2))).';
            DESIRED=desired(nn);
            e=DESIRED-CURRENT;
            if (current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)))*current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)))')>OCF_Level
                linear_equalizer=linear_equalizer+2*mu.*e.*current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)))'.'/abs(current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)))'.'*current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2))).');
            else
                linear_equalizer=linear_equalizer+2*mu.*e.*current((nn-((NTAPS_OR_PREVEQ-1)/2)):(nn+((NTAPS_OR_PREVEQ-1)/2)))'.';
            end
            ERROR=[ERROR e];
        end
        MSE=[MSE 10*log10((ERROR*ERROR')/((length(current)-((NTAPS_OR_PREVEQ-1)/2))-(((NTAPS_OR_PREVEQ-1)/2))))]; ERROR=[];
        %NORMALIZED MSE
        %MSE=[MSE 10*log10(((ERROR*ERROR')/max(ERROR*ERROR'))/((length(current)-((NTAPS_OR_PREVEQ-1)/2))-(((NTAPS_OR_PREVEQ-1)/2))))]; ERROR=[];
    end
    channel_model=linear_equalizer;
end