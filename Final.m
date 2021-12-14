
clear; close all;

fc = 1600; % carrier frequency
symbol_rate = 2400; % symbol rate


%% Read in 256 samples at 8229 Hz.

M = readmatrix('c0.txt')';
M = M(1: end-2)';

%% 1a.
Ck = zeros(1, 512);
Xq = 1:(8229/9600):256;
X = 1:1:256;
Vq = spline(X, M, Xq);

for k = 1:298
    Ck(k) = Vq(k);
end

%% 1b FFT to get the frequency Domain Representation
Cf = fft(Ck);

f_sampling = 9600; % Now the signal is sampled at 9600 Hz

Extended_X = linspace(-pi/2, pi/2, 512);
Extended_X = Extended_X .* 9600 ./ (2 .* pi);

% Plot the absolute values as Figure 1
figure();
plot(Extended_X, abs(Cf));
title('Figure 1: 9600 Hz and 512 Samples', FontSize=15);

%% 1c 
for i = 256:512
    Cf(i) = 0; % Zero out the negative frequencies
end

%% 1d

% Calculate the shift variable, or by how much do we have to shift
shifts = round((-fc .* 512)./f_sampling); 

% Truncate the negative frequencies
% truncated_Cf = Cf(257:end); 

% Cyclically shifted 
Shifted_Cf = circshift(abs(Cf), shifts);


% shifted_X = linspace(-pi/2, pi/2, 512);
% shifted_X = shifted_X .* 9600 ./ (2 .* pi);

X_2 = linspace(0, 2.*pi, 512); % X-axis for Figure2



figure();
plot(X_2, abs(Shifted_Cf)); % Shifted figure with the zeroed out negative frequencies
title('Figure 2: Cyclically shifted', FontSize=15);

%% 1e
% Compute a causal square root raised cosine filter

% s is a 1 * 129 long array
s = sqrtcos_v2(0.15, 33, 4)';


% Take the FFT
FFT_ZF_s = fft(s, 512);

%% 1f
% Multiply in frequency domain

% multiply in frequency domain
H = FFT_ZF_s .* Shifted_Cf; 
% For later use

% X_3 = 1:1:512;


% Create another figure
figure();
plot(Extended_X, abs(H));
title('Figure 3: Multiplication in frequency', FontSize=15);

%% 1g
% Take the IFFT

% Inverse fourier transform
h = ifft(H);

% Downsample it by 4 to 2400 Hz
h_down = downsample(h, 4);

% Take the absolute value, and plot it as Figure 4
abs_h_down = abs(h_down);
figure();
X_4 = 1:1:(length(abs_h_down));
plot(X_4, abs_h_down);
title('Figure 4: After IFFT and downsampling', FontSize=15);

% Compute the sum of squares
h_power = sum(abs_h_down.^2);




% Take the fft of the 128-point complex channel
FFT_down_h = fft(h_down);

% Invert it in the frequency domain
Inverted_freq = 1./FFT_down_h;

% Plot the inverted as Figure.5
X_5 = 1:1:128;
figure();
plot(X_5, abs(Inverted_freq));
grid on;
title('Figure5: Inverted complex channel', FontSize=15);


% Take the IFFT
ifft_inv_Freq = ifft(Inverted_freq);

% Convolve with the channel
Impulse_res = conv(h_down, ifft_inv_Freq);

% Plot the impulse response as Figure.6
X_6 = 1:1:255;
figure();
plot(X_6, abs(Impulse_res));
title('Figure6: Impulse response', FontSize=15);

% Find the maximum squared absolute value
max_sqr_abs = max(abs(Impulse_res).^2);
sum_of_sqr_abs = sum(abs(Impulse_res.^2));
SIR = max_sqr_abs./(sum_of_sqr_abs-max_sqr_abs);






% Noise variance
nvar = 1/(h_power.* 2000);

% Initialize the R
R = zeros(128, 128);


for i = 1:128
    for j = 1:128
        for k = 1:256
            if (0 < k+1-i) && (k+1-i < 129) && (0 < k+1-j) && (k+1-j < 129)
                if i == j
                    R(i,j) = R(i,j) + conj(h_down(k+1-i)) * h_down(k+1-j) + nvar;
                else
                    R(i,j) = R(i,j) + conj(h_down(k+1-i)) * h_down(k+1-j);
                end
            end
        end
    end
end


hd = zeros(128,1);

d = 68;

for i = 1:d
      hd(i) = h_down(d+1-i);
end




% The equalizer
c = R\conj(hd);

C = fft(c);

% Convolve the equalizer with the channel
convol = conv(c, h_down);

figure();
plot(abs(C));
title('Figure7: Magnitudes of the fft of c', fontsize=15);

figure();
plot(abs(ifft(C)));
title('Figure7b: IFFT of C', Fontsize=15);


% Calculate the SIR for 2b

max_2b = max(abs(convol));
sum_of_sqrs_2b = sum(abs(convol).^2);
SIR_2b = abs(max_2b / (sum_of_sqrs_2b -  max_2b));



delay = d;
w = zeros(128, 1); % Equalizer
w(d) = 1;
del = 2/(128.*h_power.*5);
tdata = zeros(128, 1);

% channel output array
u = zeros(128, 1);

% Initialize the squared-error
squared_error = 0;

num_of_iters = 1000000;

% choose point randomly

real_values = [-1./sqrt(2), 1./sqrt(2)];
ima_values = [-1i./sqrt(2), 1i./sqrt(2)];


for i = 1:num_of_iters
    rand_real = randi([1 2], 1, 1);
    rand_ima = randi([1 2], 1, 1);

    generated_data = real_values(rand_real) + ima_values(rand_ima);
    tdata = circshift(tdata, 1);
    tdata(1) = generated_data;
    
    % Compute the channel output as the sum of h(i) .* data(i)
    % Shift u down first
    u = circshift(u, 1);
    u(1) = h_down * tdata + normrnd(0, sqrt(nvar));
    
    

    equalizer_out = sum(w.*u);
    E = tdata(delay) - equalizer_out;

    if (i >= (num_of_iters - num_of_iters .* 0.05))
        squared_error = squared_error + abs(E).^2;
    end

    w = w + E .* del .* conj(u);

end

SINR = 0.05 .* num_of_iters ./ squared_error;




% Plot Figure.8

figure();
title('Figure8: magnitude of the fft of w against the magnitudes of the fft of c',Fontsize=15);
hold on;
plot(abs(fft(w)));
plot(abs(fft(c)));
hold off;


% 3b

delay = d;
w = zeros(128, 1); % Equalizer
w(d) = 1;
del = 2/(128.*h_power.*5);
tdata = zeros(128, 1);

% channel output array
u = zeros(128, 1);

% Initialize the squared-error
squared_error = 0;

num_of_iters = 1000000;

% choose point randomly

real_values = [-1./sqrt(2), 1./sqrt(2)];
ima_values = [-1i./sqrt(2), 1i./sqrt(2)];


for i = 1:num_of_iters
    rand_real = randi([1 2], 1, 1);
    rand_ima = randi([1 2], 1, 1);

    generated_data = real_values(rand_real) + ima_values(rand_ima);
    tdata = circshift(tdata, 1);
    tdata(1) = generated_data;
    
    % Compute the channel output as the sum of h(i) .* data(i)
    % Shift u down first
    u = circshift(u, 1);
    u(1) = h_down * tdata + normrnd(0, sqrt(nvar)) + sqrt(10) .* exp(1i .* pi .* i./3);
    
    

    equalizer_out = sum(w.*u);
    E = tdata(delay) - equalizer_out;

    if (i >= (num_of_iters - num_of_iters .* 0.05))
        squared_error = squared_error + abs(E).^2;
    end

    w = w + E .* del .* conj(u);

end

SINR = 0.05 .* num_of_iters ./ squared_error;




% Plot Figure.9

figure();
title('Figure9: fft of w adn compared with the noise-only case', Fontsize = 15);
hold on;
plot(abs(fft(w)));
plot(abs(fft(c)));
hold off;




function h=sqrtcos_v2( beta,span,sps )
% produces a vector of raised cosine filter coefficients 
% with roll-off factor beta, covering span symbols (odd), sampled at
% sps samples per symbol.
% This is a direct implementation of the equations for a filter 
% centered at the midpoint of the vector.
length=(span-1)*sps+1;
h=zeros(length,1);
b=zeros(length,1);
ts=sps; % this replaces t/sps*ts by i/ts; where Ts appears alone, have set
% the value to 1 for proper amplitude scaling.
mid=(length+1)/2;
h(mid)=(1-beta+4*beta/pi);
for i=1:mid-1
    j=mid-i;
    if (abs(i-ts/(4*beta))<=0.00001) % l'hopital rule region
        h(j)=beta*((1+2/pi)*sin(pi/(4*beta))+(1-2/pi)*cos(pi/(4*beta)));
        h(j)=h(j)/sqrt(2);
        h(mid+i)=h(j);
    else
        topj=sin(pi*-i/ts*(1-beta))+4*beta*-i/ts*cos(pi*-i/ts*(1+beta));
        h(j)=topj/(pi*-i/ts*(1-(4*beta*-i/ts)^2));
        topi=sin(pi*i/ts*(1-beta))+4*beta*i/ts*cos(pi*i/ts*(1+beta));
        h(mid+i)=topi/(pi*i/ts*(1-(4*beta*i/ts)^2));
    end
end
end
