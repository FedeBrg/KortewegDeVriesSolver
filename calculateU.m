function U=calculateU(order, delta_t, k, currentU)
    % delta_t 1 y 2   - delta_t / 2  3 y 4 - delta_t / 4  . 5 y 6
    U = currentU;
        h =  delta_t / (ceil(order / 2));

%     
%     if (order == 1)
% %         h = delta_t;
%         U = U.*exp(1i*h*(k.^3));    
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%     elseif (order == 2)
% %         h = delta_t;
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%         U = U.*exp(1i*h*(k.^3));    
%     elseif (order == 3)
% %         h = delta_t / 2;
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%         U = U.*exp(1i*h*(k.^3));
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%         U = U.*exp(1i*k.^3*h);
%     elseif(order == 4)
% %         h = delta_t / 2;
%         U = U.*exp(1i*h*(k.^3));    
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%         U = U.*exp(1i*h*(k.^3));    
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%     elseif(order == 5)
% %         h = delta_t / 4;
%         U = U.*exp(1i*h*(k.^3));    
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%         U = U.*exp(1i*h*(k.^3));    
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%         U = U.*exp(1i*h*(k.^3));    
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%     elseif(order == 6)
% %         h = delta_t / 4;
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%         U = U.*exp(1i*h*(k.^3));
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%         U = U.*exp(1i*h*(k.^3));
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%         U = U.*exp(1i*h*(k.^3));
%     elseif(order == 7)
% %         h = delta_t / 8;
%         U = U.*exp(1i*h*(k.^3));    
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%         U = U.*exp(1i*h*(k.^3));    
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%         U = U.*exp(1i*h*(k.^3));    
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%         U = U.*exp(1i*h*(k.^3));    
%         U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
%     end

    % delta_t 1 y 2   - delta_t / 2  3 y 4 - delta_t / 4  . 5 y 6
    for j = 1:(ceil(order / 2))
        if mod(order, 2) == 1
           U = U.*exp(1i*h*(k.^3));    
           U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
        else 
            U = U  - (3i*k*h).*fft((real(ifft(U))).^2);
            U = U.*exp(1i*h*(k.^3));
        end
    end
end