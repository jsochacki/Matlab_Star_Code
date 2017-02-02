function [A_plus]=left_mp_matrix_pseudoinverse(A)
 
    %Perform the Left Moore-Penrose Matrix Pseudoinverse
    %NOTE: inv(X) uses singular value decomposition of X using
    %the QR solver as opposed to pinv(X) which uses svc with a
    %default tolerance of max(size(A))*norm(A)*eps and
    %pinv(X,tolerance) does the same but treates all values less than
    %tolerance as zero (i.e. it sets them equal to zero exactly)
    %This is why I do it this way
    A_plus=inv(A'*A)*A';
    
    %NOTE: The Q-less qr least squares soltion to A*x=b is
    %if issparse(A), R = qr(A); 
    % else R = triu(qr(A)); end
    % x = R\(R'\(A'*b));
    % r = b - A*x;
    % err = R\(R'\(A'*r));
    % x = x + err;

end