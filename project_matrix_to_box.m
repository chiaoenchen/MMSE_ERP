function X_out = project_matrix_to_box( X_in, X_lb, X_ub )
[dim1,dim2]=size(X_in);


X_out=X_in;
for I=1:dim1
    for J=1:dim2
        if X_in(I,J)>X_ub(I,J)
            X_out(I,J)=X_ub(I,J);
        elseif X_in(I,J)<X_lb(I,J)
            X_out(I,J)=X_lb(I,J);
        else
        end
        
    end
end

