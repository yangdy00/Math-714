function w = Derivativefft(v, range)
v_hat = fft(v);
N = length(v);
if mod(N,2) == 0
    w_hat = 1i*[0:N/2-1 0 -N/2+1:-1]' .* v_hat;
else
    w_hat = 1i*[0:(N-1)/2 (1-N)/2:-1]' .* v_hat;
end
w = real(ifft(w_hat))*2*pi/range;
