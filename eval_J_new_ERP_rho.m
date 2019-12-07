function value=eval_J_new_ERP_rho(X,Mk_matrix,V,Gk_matrix,Sk_matrix, Nrk, K, L, P, var_w, rho_mpec)
mse=0;
for k=1:K
    mse=mse+ norm(Mk_matrix(:,:,k)*Gk_matrix(:,:,k)*X- kron(ones(Nrk,1),Sk_matrix(:,:,k)),'fro')^2+L*var_w/2*norm(Mk_matrix(:,:,k),'fro')^2;
end

value=mse+rho_mpec*(P*L-trace(X.'*V));



