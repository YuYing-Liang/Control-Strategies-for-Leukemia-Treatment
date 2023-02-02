function [f] = calc_p_star(A,b,C,x_init,rho,mn,u,k)
    f=0;
    w=[1,1,1];
    x_k = x_init;
    for i=1:k
        x_k = A*x_k + b*u(i);
        y_k = C*x_k;
        y_obs = rho(:,:,i) + mn(:,:,i);
        for j=1:length(y_k)
            f = f + ((y_k(j,1) - y_obs(j,1))/w(j))^2;
        end
    end
    f = (1/k)*f;
end