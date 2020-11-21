function kdv_strang()
% parpool('local',1);
    tic


%disp(orderKDV);
set(gca,'FontSize',18)
set(gca,'LineWidth',2)

Error =[];

N = 256;
x = linspace(-10,10,N);
delta_x = x(2) - x(1);
delta_k = 2*pi/(N*delta_x);

k = [0:delta_k:(N/2-1)*delta_k,0,-(N/2-1)*delta_k:delta_k:-delta_k];
c_1 = 13;

v = 1/2*c_1*(sech(sqrt(c_1)*(mod(x+3,20)-10)/2)).^2;
sol = @(x,t) (1/2*c_1*(sech(sqrt(c_1)*(mod(x+3-c_1*(t), 20)-10)/2)).^2);

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

tmax = 1; nplt = floor((tmax/100)/delta_t); nmax = round(tmax/delta_t);
vdata = v.'; tdata = 0;

V=fft(v);


% spmd(1)
time = 1;
for n = 1:nmax
    t = n*delta_t;
    
    V = V.*exp(1i*k.^3*delta_t/2);
    V = V  - (3i*k*delta_t).*fft((real(ifft(V))).^2);
    V = V.*exp(1i*k.^3*delta_t/2);     
    
    if mod(n,nplt) == 0
        v = real(ifft(V));
        for i =1:N
            u2(i) = sol(x(i),t);
        end
        vdata = [vdata v.']; tdata = [tdata t];
       
        Error = [Error abs(v-u2)'];
        
        if mod(n,4*nplt) == 0
           
%             plot(x,u,'LineWidth',2)
            plot(x,v,x,u2,x, Error(:, round(n/nplt))','LineWidth',1)
            legend('Solucion Strang', 'Solucion Analitica', 'Location', 'southoutside');
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

waterfall(x,tdata(1:4:end),vdata(:,1:4:end)')
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
