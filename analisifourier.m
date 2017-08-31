if(strcmp(graphics_toolkit,"gnuplot") == 0)
  graphics_toolkit gnuplot;
endif

clear;
close all;
clc;

addpath ("fourier");

# script per rifare mente locale sull'analisi dello spettro di fourier
# usiamo "t" per indicare tempo
# usiamo "s" per indicare spazio
# usiamo "f" per indicare frequenza
# usiamo "w" per indicare pulsazione (omega)
# usiamo "minuscola" per segnale nel dominio del tempo/spazio (e.g. x1s, x2t,...)
# usiamo "maiuscola" per segnale nel dominio di fourier (e.g. X1w, X2f,...)

#################################################################################


#durata temporale del segnale osservato in secondi
t = 1

# in questo caso il numero di punti del segnale e frequenza di campionamento
# coincidono perchè 

N = fs = 512; 
s = linspace(0,N-1,N);
f0 = 1; # frequenza normalizzata
w0 = 2*pi/N;

mu = 5;
x0_s = repmat(mu,1,N);  rms0 = sqrt(1/N*sumsq(x0_s));
x1_s = sin(w0*s);       rms1 = sqrt(1/N*sumsq(x1_s));
x2_s = 2*sin(2*w0*s);   rms2 = sqrt(1/N*sumsq(x2_s));
x3_s = 3*sin(3*w0*s);   rms3 = sqrt(1/N*sumsq(x3_s));
x4_s = 4*sin(4*w0*s);   rms4 = sqrt(1/N*sumsq(x4_s));

%mu = 15;
%x0_s = repmat(mu,1,N);   rms0 = sqrt(1/N*sumsq(x0_s));
%x1_s = 10*sin(w0*s);     rms1 = sqrt(1/N*sumsq(x1_s));
%x2_s = 12*sin(2*w0*s);   rms2 = sqrt(1/N*sumsq(x2_s));
%x3_s = 13*sin(3*w0*s);   rms3 = sqrt(1/N*sumsq(x3_s));
%x4_s = 14*sin(4*w0*s);   rms4 = sqrt(1/N*sumsq(x4_s));

# singole componenti del segnale
figure;
title("singole componenti del segnale");
plot(
  s,x0_s,strcat(";",num2str(mu),";"),
  s,x1_s,";sin(w0*s);",
  s,x2_s,";sin(2*w0*s);",
  s,x3_s,";sin(3*w0*s);",
  s,x4_s,";sin(4*w0*s);"
);


# segnale composto da più sinusoidi
x_s = x0_s .+ x1_s .+ x2_s .+ x3_s .+ x4_s;
rms = sqrt(1/N*sumsq(x_s));
figure;
plot(s, x_s, ";xs;");
title("segnale composto");

%bil_fft(x_s,N);
[P,PHI] = mon_psd(x_s,N*8);

return;

figure;
stem(Re_f);



mod_f = 1/N.*mod_X_f;

figure;
stem(bilmod_f);

## spettro monolatero (modulo/fase)
mod_f = bilmod_f(1:N/2) .+ [0, bilmod_f(2:N/2)];
figure;
stem(mod_f);
title("spettro monolatero");

## spettro di potenza
# calcolo periodogramma (in db)
figure;
[pxx,w] = periodogram(x_s,[],N,"onesided");
stem(w,pxx);

# calcolo periodogramma manuale
sq_abs_X_f = abs_X_f.^2;
psd = sq_abs_X_f(1:N/2) .+ [0, sq_abs_X_f(2:N/2)];

figure;
stem(psd);


