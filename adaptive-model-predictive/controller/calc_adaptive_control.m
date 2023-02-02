function [f] = calc_adaptive_control(A, b, C, rho, x_k, u, q, r, Hp, Hu, k)
    f=0;
    for i=1:Hp
        x_k = A*x_k + b*u(i);
        z = C*x_k;
        if k+i <= length(rho(1,1,:))
            y_err = z - rho(:,:,k+i);
            f = f + q*norm(y_err);
        end
    end
    
    for i=1:Hu
        f = f + r*u(i)^2;
    end
end