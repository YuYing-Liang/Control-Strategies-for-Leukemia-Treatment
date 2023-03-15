function f = f_nonlin(x_sym,u_sym,params,t_s)
    g_grav = 9.81;
    x = num2cell(x_sym);
    p = num2cell(params);
    [x1, x2, x3, x4] = x{:};
    [m_cart, m_pend, l_pend] = p{:};
    
    f = [...
        x2;
        (l_pend*(x4^2)*sin(x3) - g_grav*cos(x3)*sin(x3) + u_sym/m_pend) / (m_cart/m_pend + (sin(x3))^2);
        x4;
        ((m_cart+m_pend)/(m_pend*l_pend)*g_grav*sin(x3) - (x4^2)*cos(x3)*sin(x3) - cos(x3)*u_sym/(m_pend*l_pend)) / (m_cart/m_pend + (sin(x3))^2);
    ];
    f = [x1;x2;x3;x4] + t_s * f;
end