function [f] = calc_adaptive_control(A, b, C, rho, x, u, q, r, Hp, Hu, k)
    f=0;
    x_k = A*x + b*u;
    for i=1:Hp
        z = C*x_k;
        if k+i > length(rho(1,1,:))
            f = f + q*norm(z - rho(:,:,length(rho(1,1,:))));
        else
            f = f + q*norm(z - rho(:,:,k+i));
        end
        x_k = A*x_k + b*u;
    end
    
    for i=1:Hu
        f = f + r*u^2;
    end
    disp(f);
end