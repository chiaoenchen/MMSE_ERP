% Run MSE based one-bit precoding
close all
clear all
clc
format long
j=sqrt(-1);

Nt=128;
K=16;
Nrk=1;
L=10;


rho_mpec=1; % penalty parameter
kappa_mpec=2; % multiply rho by kappa every L iterations 
Psi_mpec=10; 

% const_size=4; % QPSK
% a=sqrt(6/(const_size-1));
% OMEGA=0.5*a*[-1+j;1+j;-1-j;1-j]; % set of 4 QAM symbols
% BITS=[0,1;
%     1,1;
%     0,0;
%     1,0]; % the corresponding bits of OMEGA


% 16-QAM modulation
const_size=16;
a=sqrt(6/(const_size-1));
OMEGA=0.5*a*[-3+3j;-1+3j;1+3j;3+3j;-3+j;-1+j;1+j;3+j;-3-j;-1-j;1-j;3-j;-3-3j;-1-3j;1-3j;3-3j]; % set of 16 QAM symbols
BITS=[0,0,1,0;
    0,1,1,0;
    1,1,1,0;
    1,0,1,0;
    0,0,1,1;
    0,1,1,1;
    1,1,1,1;
    1,0,1,1;
    0,0,0,1;
    0,1,0,1;
    1,1,0,1;
    1,0,0,1;
    0,0,0,0;
    0,1,0,0;
    1,1,0,0;
    1,0,0,0]; % the corresponding bits of OMEGA



P=1;
NOE=1000;
SNR_DB_RANGE=0:2:16;

total_bits=zeros(size(SNR_DB_RANGE)); % total number of bits decoded
total_errors=zeros(size(SNR_DB_RANGE)); %total number of error bits detected
ber_array=zeros(size(SNR_DB_RANGE));
snr_index=1;
for snr_dB=SNR_DB_RANGE
    randn('state',0);
    rand('state',0);
    snr_dB
    var_w=P*10^(-0.1*snr_dB); % SNR defined as total power divided by noise variance
    while total_errors(1,snr_index)<=NOE        
        % all channel matrix
        Hc=sqrt(0.5)*randn(Nt,Nrk,K)+j*sqrt(0.5)*randn(Nt,Nrk,K); % channel matrix                       
        Gk_matrix=zeros(2*Nrk,2*Nt,K);                
        Wk_matrix=zeros(2*Nrk,L,K);
        Yk_matrix=zeros(2*Nrk,L,K);
        for k=1:K
            Hck_herm=Hc(:,:,k)';
            for n=1:Nt
                for ik=1:Nrk
                    Gk_matrix( (ik-1)*2+1:ik*2, (n-1)*2+1:n*2,k)=complex_scalar_to_real_matrix(Hck_herm(ik,n));
                end            
            end          
            Wk_matrix(:,:,k)=sqrt(var_w)*sqrt(0.5)*randn(2*Nrk,L);
        end

        % Generate symbol bits
        sc_index_matrix=randi(length(OMEGA),L,K); % Generate the symbol index for the transmitted vector s
        sc_matrix=zeros(K,L);
        sc_bit_matrix=zeros(log2(const_size),L,K);
        Sk_matrix=zeros(2,L,K);
        for k=1:K
            for ell=1:L
                sc_matrix(k,ell)=OMEGA(sc_index_matrix(ell,k).');
                sc_bit_matrix(:,ell,k)=BITS(sc_index_matrix(ell,k),:);
                Sk_matrix(:,ell,k)=[real(sc_matrix(k,ell)); imag(sc_matrix(k,ell))];
            end
        end
        



     

        % MSE_MPEC_ERP_APG_complex_beta_k
        PG_precision=1e-5;
        [X,Mk_matrix]=MSE_MPEC_new_ERP_PG_complex_beta_k(Gk_matrix, Sk_matrix, Nt, Nrk, K, L, P, var_w, rho_mpec, kappa_mpec, Psi_mpec, PG_precision);
        
        decoded_bits_matrix=zeros(log2(const_size),L,K);       
        Y_matrix=zeros(2*Nrk,L,K);
        Yc_matrix=zeros(Nrk,L,K);
        beta_c_matrix=zeros(Nrk,K);
        for k=1:K
            Y_matrix(:,:,k)=Gk_matrix(:,:,k)*X+Wk_matrix(:,:,k);
            Yc_matrix(:,:,k)=Y_matrix(1:2:2*Nrk-1,:,k)+j*Y_matrix(2:2:2*Nrk,:,k);
            for ik=1:Nrk
                beta_c_matrix(ik,k)=Mk_matrix((ik-1)*2+1,(ik-1)*2+1,k)+j*Mk_matrix(2*ik,(ik-1)*2+1,k);
            end
            
            % AML detection
            for ell=1:L
                [~,decoded_bits_matrix(:,ell,k)]=AML_detection(Nrk, OMEGA, BITS, Yc_matrix(:,ell,k), beta_c_matrix(:,k));
            end
        end
                      
        total_bits(snr_index)=total_bits(snr_index)+K*L*log2(const_size); % total number of bits decoded
        total_errors(snr_index)=total_errors(snr_index)+sum(sum(sum(abs(decoded_bits_matrix-sc_bit_matrix)))); %total number of error bits detected
 
        
        if mod(total_bits(snr_index),100000)==0
            total_errors
        else
        end
                
        total_bits
        total_errors
        
    end % end NOE
    (total_errors./total_bits).'
    snr_index=snr_index+1;
end
ber_array=total_errors./total_bits;