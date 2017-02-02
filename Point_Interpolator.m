function [output]=Point_Interpolator(x_data,y_data,point,INTERP_CASE)

switch INTERP_CASE
    case 'Logarithmic'
        output=y_data(1)+((10*log10(point/x_data(1))*(y_data(end)-y_data(1)))...
            /(10*log10(x_data(end)/x_data(1))));
    case 'Linear'
        
    otherwise
end

end