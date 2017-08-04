function [P PHI] = mon_psd(x_n, NFFT)

  # numero campioni del segnale di ingresso
  N = length(x_n);
  
  # trasformata MONOLATERA di fourier: se N non è una potenza di 2, viene fatto zero-padding
  X_q = fft(x_n, NFFT);
  Re_q = real(X_q);
  Im_q = imag(X_q);

  # modulo fase
  MOD = abs(X_q);
  PHASE = arg(X_q);
  
  MOD_SQ = X_q.*conj(X_q);
  
  # ampiezza monolatera: l'ampiezza di ogni bin è da scalare sul numero di campioni
  # del segnale in ingresso, a meno che non sia già stato fatto zero-padding; in
  # tal caso il fattore di scala è errato e dopo aver ottenuto A da questa funzione
  # bisogna ri-moltiplicare tutto per N/length(x_n_senza_zero_padding)
  # devo sommare le due metà delo spettro.
  P_left = MOD_SQ(1:NFFT/2+1);
  P_right = fliplr(MOD_SQ)(1:NFFT/2);
  P = 1/(N.*NFFT) .* ( P_left .+ [0 P_right]); 
  PHI = rad2deg(PHASE)(NFFT/2:end);
 
  idx_q = linspace(0,NFFT-1,NFFT);
  omega_q = 2*pi/NFFT .* idx_q;
  
  # grafico lo spettro bilatero  
  figure;
  stem(omega_q, MOD, ";MOD;", "color", "r", "markersize", 3);
  grid on;
  
  # grafico lo spettro monolatero
  
  idx_q = linspace(0,NFFT/2,NFFT/2+1);
  omega_q = 2*pi/NFFT .* idx_q;
  
  figure;
  stem(omega_q, P, ";P;", "markersize", 3);
  grid on;
  
  RMS_PSD = sqrt(P);
  
  figure;
  stem(omega_q, RMS_PSD, ";scaled_P;", "markersize", 3);
  grid on;
  
%  figure;  
%  stem(omega_q, PHI);
%  grid on;
  
endfunction