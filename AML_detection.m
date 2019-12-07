function [symbol_out,bits]=AML_detection(Nrk, OMEGA, BITS, yc_vector, beta_c_vector)

L_OMEGA=length(OMEGA);
distance_table=zeros(L_OMEGA,1);
for I=1:L_OMEGA
    s=OMEGA(I);
    for ik=1:Nrk
        distance_table(I)=distance_table(I)+abs(yc_vector(ik)-s/beta_c_vector(ik))^2;
    end    
end

[~,min_index]=min(distance_table(:));

bits=BITS(min_index,:);
symbol_out=OMEGA(min_index,:);
end

