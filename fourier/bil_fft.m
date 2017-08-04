function [A PHI X_q] = bil_fft(x_n, NFFT)

  # numero campioni del segnale di ingresso
  N = length(x_n);
  
  # trasformata BILATERA di fourier: se N non è una potenza di 2, viene fatto zero-padding
  X_q = fft(x_n, NFFT);
  Re_q = real(X_q);
  Im_q = imag(X_q);

  # modulo fase
  MOD = abs(X_q);
  PHASE = arg(X_q);
  
  # ampiezze bilatere: l'ampiezza di ogni bin è da scalare sul numero di campioni
  # del segnale in ingresso, a meno che non sia già stato fatto zero-padding; in
  # tal caso il fattore di scala è errato e dopo aver ottenuto A da questa funzione
  # bisogna ri-moltiplicare tutto per N/length(x_n_senza_zero_padding)
  A = 1/N .* MOD;  
  PHI = rad2deg(PHASE);
  
  idx_q = [-NFFT/2+1 : 1 : NFFT/2 ];
  omega_q = 2*pi/NFFT .* idx_q;
  
  # grafico lo spettro bilatero
  figure;
  stem(omega_q, A);
  grid on;
  figure;  
  stem(omega_q, PHI);
  grid on;
  
endfunction