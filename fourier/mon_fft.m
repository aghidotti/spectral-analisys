function [A PHI] = mon_fft(x_n, NFFT)

  # numero campioni del segnale di ingresso
  N = length(x_n);
  
  # trasformata MONOLATERA di fourier: se N non � una potenza di 2, viene fatto zero-padding
  X_q = fft(x_n, NFFT);
  Re_q = real(X_q);
  Im_q = imag(X_q);

  # modulo fase
  MOD = abs(X_q);
  PHASE = arg(X_q);
  
  # ampiezza monolatera: l'ampiezza di ogni bin � da scalare sul numero di campioni
  # del segnale in ingresso, a meno che non sia gi� stato fatto zero-padding; in
  # tal caso il fattore di scala � errato e dopo aver ottenuto A da questa funzione
  # bisogna ri-moltiplicare tutto per N/length(x_n_senza_zero_padding)
  # devo sommare le due met� delo spettro.
  A_left = MOD(1:NFFT/2+1);
  A_right = fliplr(MOD)(1:NFFT/2);
  A = 1/N .* ( A_left .+ [0 A_right]); 
  PHI = rad2deg(PHASE)(NFFT/2:end);
 
  idx_q = linspace(0,NFFT/2,NFFT/2+1);
  omega_q = 2*pi/NFFT .* idx_q;
  
  # grafico lo spettro bilatero
  figure;
  stem(omega_q, A, "markersize", 3);
  grid on;
  title("spettro bilatero di ampiezza");
%  figure;  
%  stem(omega_q, PHI);
%  grid on;
  
endfunction