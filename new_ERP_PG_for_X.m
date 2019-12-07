function X_out = new_ERP_PG_for_X(Gk_matrix, Sk_matrix, X, Mk_matrix, V, Nt, Nrk, K, L, P, var_w, rho_mpec, PG_precision)
alpha=0.25;
lambda=0.5;
X_n=X;

debug_old=eval_J_new_ERP_rho(X_n,Mk_matrix,V,Gk_matrix,Sk_matrix, Nrk, K, L, P, var_w, rho_mpec);
%debug_array=debug_old;
for n=0:10000    
    Z=X_n;
    
    Delta_n=zeros(2*Nt,L);
    for k=1:K
        MkGk=Mk_matrix(:,:,k)*Gk_matrix(:,:,k);
        Delta_n=Delta_n+2*MkGk.'*(MkGk*Z- kron( ones(Nrk,1) , Sk_matrix(:,:,k))  );
    end
    Delta_n=Delta_n-rho_mpec*V;
    t=1;
        
    X_np1=project_matrix_to_box( Z-t*Delta_n, -sqrt(P/2/Nt)*ones(2*Nt,L), sqrt(P/2/Nt)*ones(2*Nt,L) );
    Qt_n=(Z-X_np1)/t;
        
    m=0;
    while eval_J_new_ERP_rho(X_np1,Mk_matrix,V,Gk_matrix,Sk_matrix, Nrk, K, L, P, var_w, rho_mpec)> (eval_J_new_ERP_rho(Z,Mk_matrix,V,Gk_matrix,Sk_matrix, Nrk, K, L, P, var_w, rho_mpec)-alpha*t*trace(Delta_n.'*Qt_n))
        t=lambda*t;
        X_np1=project_matrix_to_box( Z-t*Delta_n, -sqrt(P/2/Nt)*ones(2*Nt,L), sqrt(P/2/Nt)*ones(2*Nt,L) );
        Qt_n=(Z-X_np1)/t;
        m=m+1;
    end
        
    X_n=X_np1;
        
    debug_new=eval_J_new_ERP_rho(X_np1,Mk_matrix,V,Gk_matrix,Sk_matrix, Nrk, K, L, P, var_w, rho_mpec);
    %debug_array=[debug_array;debug_new];

    
    if abs(debug_new-debug_old)/debug_new<PG_precision        
        break;
    else
        debug_old=debug_new;
    end    
    
    if n==10000
        keyboard
    else
    end
end
X_out=X_np1;


%figure
%plot(debug_array)
%debug_array
%keyboard










