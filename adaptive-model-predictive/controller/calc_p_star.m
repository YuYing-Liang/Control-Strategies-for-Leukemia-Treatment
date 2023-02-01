function [f] = calc_p_star(A,b,C,x,y,mn,u,k)
    f=0;
    for i=1:k
        x_star = A*x(:,:,i) + b*u(:,i);
        y_star = C*x_star;
        y_obs = y(:,:,i) + mn(:,:,i);
        f = f + norm(y_star - y_obs);
%         for j=1:length(y_star)
%             f = f + (y_star(j,1) - y_obs(j,1))^2;
%         end
    end
    f = (1/k)*f;
end