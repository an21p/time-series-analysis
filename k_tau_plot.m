function [] = k_tau_plot(name, price, n)
% GHE K(tau) % from Xi Ji
    % q=1
    K_UP=zeros(n,1);
    K_DOWN=zeros(n,1);
    K1=zeros(n,1);
        for tao=1:n
            for t=1:size(price,1)-tao
            K_up=abs((price(t+tao,:))-(price(t,:)));
            K_UP(tao,:)=K_UP(tao,:)+K_up;
            K_down=abs(price(t,:));
            K_DOWN(tao,:)=K_DOWN(tao,:)+K_down;
            end
            K1(tao,:)=K_UP(tao,:)./K_DOWN(tao,:);
        end
      log_tao=log(1:n);
      log_K1=log(K1);

    % q=2
    K_UP=zeros(n,1);
    K_DOWN=zeros(n,1);
    K2=zeros(n,1);
        for tao=1:n
            for t=1:size(price,1)-tao
            K_up=abs(price(t+tao,:)-price(t,:)).^2;
            K_UP(tao,:)=K_UP(tao,:)+K_up;
            K_down=abs(price(t,:)).^2;
            K_DOWN(tao,:)=K_DOWN(tao,:)+K_down;
            end
            K2(tao,:)=K_UP(tao,:)./K_DOWN(tao,:);
        end
      log_tao=log(1:n);
      log_K2=log(K2);

    % q=3
    K_UP=zeros(n,1);
    K_DOWN=zeros(n,1);
    K3 =zeros(n,1);
        for tao=1:n
            for t=1:size(price,1)-tao
            K_up=abs(price(t+tao,:)-price(t,:)).^3;
            K_UP(tao,:)=K_UP(tao,:)+K_up;
            K_down=abs(price(t,:)).^3;
            K_DOWN(tao,:)=K_DOWN(tao,:)+K_down;
            end
            K3(tao,:)=K_UP(tao,:)./K_DOWN(tao,:);
        end
      log_tao=log(1:n);
      log_K3=log(K3);

    % GHE plot: log(K(\tau)) vs log(\tau)
    figure
    plot(log_tao,log_K1,'-or',log_tao,log_K2,'-ob',log_tao,log_K3,'-ok')
    legend('q=1','q=2','q=3', 'Location' , 'best')
    title(['Hurst exponent ', name])
    xlabel('log (\tau)')
    ylabel('log (K(\tau))')
end
