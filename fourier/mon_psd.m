function [P PHI] = mon_psd(x_n, NFFT)

  # numero campioni del segnale di ingresso
  N = length(x_n);
  
  
  # trasformata BILATERA di fourier: se N non è una potenza di 2, viene fatto zero-padding
  X_q = fft(x_n, NFFT);
  Re_q = real(X_q);
  Im_q = imag(X_q);

  # modulo fase
  MOD = abs(X_q);
  PHASE = arg(X_q);
  
  MOD_SQ = X_q.*conj(X_q);  
 
  idx_q = linspace(0,NFFT-1,NFFT);
  omega_q = 2*pi/NFFT .* idx_q;
  
  # genero finestra rettangolare --> voglio vedere la sinc che si sovrappone allo spettro
  # in caso di zero padding
  WIN_n = [window(@rectwin, N); zeros(NFFT-N,1)]';
  WIN_q = fft(WIN_n, NFFT);
  MOD_WIN_q = 1/NFFT .* abs(WIN_q);
  
  # grafico trasformata della finestra
  figure; 
  plot(omega_q, MOD_WIN_q, ";MOD_WIN_q;", "color", "g");
  #stem(omega_q, MOD_WIN_q, ";MOD_WIN_q;", "color", "g", "markersize", 3);
  grid on;
  title("modulo trasformata della window");
  
%  return;
  
  # grafico lo spettro bilatero  
  figure;  
  stem(omega_q, MOD, ";MOD;", "color", "r", "markersize", 3);
  grid on;
  title("modulo spettro fourier bilatero");
    
  # grafico lo spettro bilatero  
  figure;  
  stem(omega_q, MOD_SQ,";MOD_SQ;", "color", "r", "markersize", 3);
  grid on;
  title("spettro di pontenza non normalizzato");
  
  
  # ampiezza monolatera: l'ampiezza di ogni bin è da scalare sul numero di campioni
  # del segnale in ingresso, a meno che non sia già stato fatto zero-padding; in
  # tal caso il fattore di scala è errato e dopo aver ottenuto A da questa funzione
  # bisogna ri-moltiplicare tutto per N/length(x_n_senza_zero_padding)
  # devo sommare le due metà delo spettro.
  P_left = MOD_SQ(1:NFFT/2+1);
  P_right = fliplr(MOD_SQ)(1:NFFT/2);
  P = (P_left .+ [0 P_right]);
  PHI = rad2deg(PHASE)(NFFT/2:end);
  P_scaled = 1/(N*NFFT) .* P; 
  
  # trasformata della window monolatera
  WIN_left = MOD_WIN_q(1:NFFT/2+1);
  WIN_right = fliplr(MOD_WIN_q)(1:NFFT/2);
  WIN_mono = (WIN_left .+ [0 WIN_right]);
  
  # scalo/normalizzo lo spettro
  idx_q = linspace(0,NFFT/2,NFFT/2+1);
  omega_q = 2*pi/NFFT .* idx_q;
  
  figure;
  stem(omega_q, P, ";P;", "markersize", 3);
  grid on;
  title("potenza");
  
  figure;
  stem(omega_q, P_scaled, ";P_scaled;", "markersize", 3);
  grid on;
  title("potenza scalata");
  
  RMS_PSD = sqrt(P_scaled);
  
  figure;
  hold on;
  stem(omega_q, RMS_PSD, ";PSD_rms;", "markersize", 3);
  stem(omega_q, WIN_mono, ";WIN_mono;", "color", "g", "markersize", 3);
  grid on;
  title("psd (RMS)");
  
%  figure;  
%  stem(omega_q, PHI);
%  grid on;
  
endfunction