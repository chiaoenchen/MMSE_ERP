function matrix_out=complex_scalar_to_real_matrix(z)
matrix_out=[real(z), -imag(z); imag(z), real(z)];
end

