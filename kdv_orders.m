function kdv_orders(order1,order2)
% parpool('local',1);
    tic


%disp(orderKDV);
set(gca,'FontSize',18)
set(gca,'LineWidth',2)

N = 256;
x = linspace(-10,10,N);
delta_x = x(2) - x(1);
delta_k = 2*pi/(N*delta_x);

k = [0:delta_k:(N/2-1)*delta_k,0,-(N/2-1)*delta_k:delta_k:-delta_k];
c_1 = 16;
c_2 = 4;

u = 1/2*c_1*(sech(sqrt(c_1)*(x+8)/2)).^2 + 1/2*c_2*(sech(sqrt(c_2)*(x+1)/2)).^2;

name = 'two_soliton.gif';
eval(['delete ',name])

delta_t = 0.4/N^2;
t=0;
plot(x,u,'LineWidth',2)
axis([-10 10 0 10])
xlabel('x')
ylabel('u')
text(6,9,['t = ',num2str(t,'%1.2f')],'FontSize',18)
drawnow
gif_add_frame(gcf,name,2);
drawnow

tmax = 1; nplt = floor((tmax/100)/delta_t); nmax = round(tmax/delta_t);
udata = u.'; tdata = 0;

time = 1;

for i = 1:1:order1
    Us1{i} = fft(u);
end

for i = 1:1:order2
    Us2{i} = fft(u);
end
% spmd(1)

for n = 1:nmax
    t = n*delta_t;
    %ORDER 1
    for i = 1:order1
        Us1{i} = calculateU(i, delta_t, k, Us1{i});
    end
    
    gamma1 = 2*getGamma(order1);
    U1 = 0;
    for i = 1:order1
        U1 = U1 + gamma1(i)*Us1{i};
    end
    
    %ORDER 2
    for i = 1:order2
        Us2{i} = calculateU(i, delta_t, k, Us2{i});
    end
    
    gamma2 = 2*getGamma(order2);
    U2 = 0;
    for i = 1:order2
        U2 = U2 + gamma2(i)*Us2{i};
    end
    %%%%%%%
    if mod(n,nplt) == 0
        u = real(ifft(U1));
        u2 = real(ifft(U2));

        udata = [udata u.']; tdata = [tdata t];
        if mod(n,4*nplt) == 0
            subplot(1,2,1)
            t2(time) = t;
            error(time) = mean(abs(u2 - u));
            plot(x,u,'LineWidth',2)
            axis([-10 10 0 10])
            xlabel('x')
            ylabel('u')
            text(6,9,['t = ',num2str(t,'%1.2f')],'FontSize',10)
            subplot(1,2,2)
            plot(t2,error);
            xlabel('time[s]')
            ylabel('Mean square Error')
            time = time + 1;
            drawnow
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
end