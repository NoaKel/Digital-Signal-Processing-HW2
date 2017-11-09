function [ MDFT ] = fft_windows( x ,N )

% adds zeros if necessary
modulo = mod(length(x),N);
if (modulo ~= 0)
    x = x(:);
    x = [x; zeros(length(x)+N-modulo,1)];
end


% creats MDFT matrix
MDFT = zeros(length(x)/N,N);

for i=1:(length(x)/N)
    x_segment = x(((i-1)*N+1:i*N));
    MDFT(i,:) = fftshift(fft(x_segment));
end

end

