function kdv_yoshida()
% parpool('local',1);
    tic


%disp(orderKDV);
set(gca,'FontSize',18)
set(gca,'LineWidth',2)

Error = [];

N = 256;
x = linspace(-10,10,N);
delta_x = x(2) - x(1);
delta_k = 2*pi/(N*delta_x);

k = [0:delta_k:(N/2-1)*delta_k,0,-(N/2-1)*delta_k:delta_k:-delta_k];
c_1 = 13;

u = 1/2*c_1*(sech(sqrt(c_1)*(mod(x+3,20)-10)/2)).^2;

sol = @(x,t) (1/2*c_1*(sech(sqrt(c_1)*(mod(x+3-c_1*(t), 20)-10)/2)).^2);

w1 = 1/(2-2^(1/3));
w0 = 1-2*w1;

name = 'two_soliton.gif';
eval(['delete ',name])

delta_t = 0.4/N^2;
t=0;
% plot(x,u,'LineWidth',2)
% axis([-10 10 0 10])
% xlabel('x')
% ylabel('u')
% text(6,9,['t = ',num2str(t,'%1.2f')],'FontSize',18)
% drawnow
% gif_add_frame(gcf,name,2);
% drawnow

tmax = 0.5; nplt = floor((tmax/100)/delta_t); nmax = round(tmax/delta_t);
udata = u.'; tdata = 0;

U=fft(u);

% spmd(1)
time = 1;
for n = 1:nmax
    t = n*delta_t;
    
    U = U.*exp(1i*k.^3*delta_t*w1/2);
    U = U  - (3i*k*delta_t*w1).*fft((real(ifft(U))).^2);
    U = U.*exp(1i*k.^3*delta_t*(w1+w0)/2);
    U = U  - (3i*k*delta_t*w0).*fft((real(ifft(U))).^2);
    U = U.*exp(1i*k.^3*delta_t*(w1+w0)/2);
    U = U  - (3i*k*delta_t*w1).*fft((real(ifft(U))).^2);
    U = U.*exp(1i*k.^3*delta_t*w1/2);
    
    
    if mod(n,nplt) == 0
        u = real(ifft(U));
        for i =1:N
            u2(i) = sol(x(i),t);
        end
        udata = [udata u.']; tdata = [tdata t];
        
        Error = [Error abs(u-u2)'];
        
        if mod(n,4*nplt) == 0
           
%             plot(x,u,'LineWidth',2)
            plot(x,u,x,u2,'LineWidth',1)
            legend('Solucion Yoshida', 'Solucion Analitica', 'Location', 'southoutside');
            axis([-10 10 0 10])
            xlabel('x')
            ylabel('u')
            text(6,9,['t = ',num2str(t,'%1.2f')],'FontSize',18)
            drawnow
            %gif_add_frame(gcf,name,2);
            
            
        end
    end
end

figure

waterfall(x,tdata(1:4:end),udata(:,1:4:end)')
colormap(1e-6*[1 1 1]); view(-20,25)
xlabel x, ylabel t, axis([-10 10 0 tmax 0 10]), grid off
zlabel('u')
set(gca,'ztick',[0 10]), pbaspect([1 1 .13])
print -djpeg two_soliton

% end
toc

% lo del error 
figure
    plot(max(Error)), hold on
    plot(mean(Error)), 
    legend('Error global', 'Media error', 'Location', 'southoutside'),
    xlabel(['Dt = ', num2str(delta_t, '%1.5g')]);
end
