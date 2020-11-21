function kdv_solve_parallel(orderKDV)
%parpool(orderKDV);

delete(gcp('nocreate'));

myCluster = parcluster('local');
numOfCores = orderKDV;

myCluster.NumWorkers = numOfCores;  % 'Modified' property now TRUE
saveProfile(myCluster);    % 'local' profile now updated,
disp(numOfCores);                           % 'Modified' property now FALSE
parpool('local', numOfCores);


tic
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

U = fft(u);

spmd(orderKDV)
    for n = 1:nmax
        if (labindex == 1)
            t = n*delta_t;

            for i = 2:orderKDV
                labSend(U,i);
            end
            U = calculateU(labindex,delta_t,k,U);
            gamma = 2*getGamma(orderKDV);
            U = gamma(1)*U;
            for i = 2:orderKDV
                U = U + gamma(i)*labReceive(i);
            end
            
            if mod(n,nplt) == 0
                u = real(ifft(U));
                udata = [udata u.']; tdata = [tdata t];
                if mod(n,4*nplt) == 0
                    plot(x,u,'LineWidth',2)
                    axis([-10 10 0 10])
                    xlabel('x')
                    ylabel('u')
                    text(6,9,['t = ',num2str(t,'%1.2f')],'FontSize',18)
                    drawnow
                    gif_add_frame(gcf,name,2);
                end
            end
        end
        
        if (labindex ~= 1)
            U = labReceive(1);
            U = calculateU(labindex,delta_t,k,U);
            labSend(U,1);
        end
        
    end
if(labindex == 1)
    figure
    waterfall(x,tdata(1:4:end),udata(:,1:4:end)')
    colormap(1e-6*[1 1 1]); view(-20,25)
    xlabel x, ylabel t, axis([-10 10 0 tmax 0 10]), grid off
    zlabel('u')
    set(gca,'ztick',[0 10]), pbaspect([1 1 .13])
    print -djpeg two_soliton
end

end

toc
end