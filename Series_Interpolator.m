function [output]=Series_Interpolator(x_data,y_data,desired_points,INTERP_CASE)

    OTURN=0;
    if size(x_data,2) < size(x_data,1), x_data=x_data.';, end;
    if size(y_data,2) < size(y_data,1), y_data=y_data.';, end;
    if size(desired_points,2) < size(desired_points,1), OTURN=1; desired_points=desired_points.';, end;

nn=1;
for n=1:1:length(desired_points)
    if n==1
        output(n)=y_data(n);
    else
        output(n)=Point_Interpolator(x_data(nn:1:nn+1),y_data(nn:1:nn+1),desired_points(n),INTERP_CASE);
        if desired_points(n)>=x_data(nn+1)
            nn=nn+1;
        end
    end
end

if OTURN, output=output.';, end;

end