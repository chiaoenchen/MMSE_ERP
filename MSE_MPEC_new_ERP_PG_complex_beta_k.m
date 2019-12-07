function [X,Mk_matrix]=MSE_MPEC_new_ERP_PG_complex_beta_k(Gk_matrix, Sk_matrix, Nt, Nrk, K, L, P, var_w, rho_mpec, kappa_mpec, Psi_mpec, PG_precision)
j=sqrt(-1);

%initialization
X=sqrt(P/2/Nt)*ones(2*Nt,L);
Mk_matrix=zeros(2*Nrk,2*Nrk,K);
beta_c_vector=ones(Nrk,K);

for k=1:K
    for ik=1:Nrk
        Mk_matrix((ik-1)*2+1:2*ik, (ik-1)*2+1:2*ik, k)=complex_scalar_to_real_matrix(beta_c_vector(ik,k));
    end
end



V=zeros(2*Nt,L);
%debug_xv=[];
for psi=0:1000
    FLAG=1;
    debug_old=eval_J_new_ERP_rho(X,Mk_matrix,V,Gk_matrix,Sk_matrix, Nrk, K, L, P, var_w, rho_mpec);
    %debug_array=debug_old;

    while FLAG
        X=new_ERP_PG_for_X(Gk_matrix, Sk_matrix, X, Mk_matrix, V, Nt, Nrk, K, L, P, var_w, rho_mpec, PG_precision);  

        
        beta_c_vector=ones(Nrk,K);
        v_r_c=zeros(2,Nrk,L,K);
        for k=1:K
            for ik=1:Nrk
                num_k_ik=0;
                den_k_ik=0;                
                for ell=1:L
                    v_r_c(:,ik,ell,k)=Gk_matrix(2*(ik-1)+1:2*ik,:,k)*X(:,ell);
                    r_c_ik_ell_k=v_r_c(1,ik,ell,k)+j*v_r_c(2,ik,ell,k);
                    num_k_ik=num_k_ik+r_c_ik_ell_k'*(Sk_matrix(1,ell,k)+j*Sk_matrix(2,ell,k));
                    den_k_ik=den_k_ik+r_c_ik_ell_k'*r_c_ik_ell_k;                   
                end
                beta_c_vector(ik,k)=num_k_ik/(L*var_w+den_k_ik);
            end
        end
        
         

                   
        for k=1:K
            for ik=1:Nrk
                Mk_matrix((ik-1)*2+1:2*ik, (ik-1)*2+1:2*ik, k)=complex_scalar_to_real_matrix(beta_c_vector(ik,k));
            end
        end
        
        
        debug_new=eval_J_new_ERP_rho(X,Mk_matrix,V,Gk_matrix,Sk_matrix, Nrk, K, L, P, var_w, rho_mpec);
        %debug_array=[debug_array; debug_new]              
        %keyboard
        
        if abs(debug_new-debug_old)/debug_new<PG_precision
            FLAG=0;
            
        else
            debug_old=debug_new;
            %debug_old
            %keyboard
        end 

    end

    
    %figure
    %plot(debug_array)
    
    V=sqrt(P*L)*X/norm(X,'fro');
    
    % test
    %debug_xv=[debug_xv; trace(X.'*V)];

    %keyboard
    
    if max(abs(vec(X)-vec(V)))<1e-5
        break;
    else        
    end
    
    if mod(psi,Psi_mpec)==0
        rho_mpec=kappa_mpec*rho_mpec;
    else
    end
        
    if psi==1000
        keyboard
    else
    end
    
    %psi
    
end


%rho_mpec
%figure
%plot(debug_xv)
%keyboard


