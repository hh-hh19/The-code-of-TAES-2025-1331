function dY = uavWithObserverAndControllerDynamics(t, Y, params, neighbors_in)

    N = size(neighbors_in, 1);
    K = 21;
    dY = zeros(N * K, 1);
    l_1   = params.l_1;    l_2   = params.l_2;
    zeta_1= params.zeta_1; zeta_2= params.zeta_2;
    zeta_3= params.zeta_3; zeta_4= params.zeta_4; zeta_5= params.zeta_5;
    r_s   = params.r_s;    r_c   = params.r_c;    ita   = params.ita;
    if ~isfield(params,'traj_type'), params.traj_type = 'curve'; end
    if ~isfield(params,'target_v'),  params.target_v  = [45;45;5]; end
    vtx = params.target_v(1); vty = params.target_v(2); vtz = params.target_v(3);
    neighbors = zeros(N, N);
    for i = 1:N
        for j = i+1:N
            bi = (i-1) * K; bj = (j-1) * K;
            xi = Y(bi+1); yi = Y(bi+2); zi = Y(bi+3);
            xj = Y(bj+1); yj = Y(bj+2); zj = Y(bj+3);
            distance = norm([xi-xj, yi-yj, zi-zj]);
            if distance <= r_c
                neighbors(i,j) = 1; neighbors(j,i) = 1;
            end
        end
    end
    neighbors = convert_adjacency_matrix(neighbors); 

    for i = 1:N
        base = (i-1)*K;
        Yi = Y(base+1:base+K);
        x=Yi(1); y=Yi(2); z=Yi(3);
        Vx=Yi(4); Vy=Yi(5); Vz=Yi(6);

        px_hat=Yi(7);  py_hat=Yi(8);  pz_hat=Yi(9);
        iotax=Yi(10);  iotay=Yi(11);  iotaz=Yi(12);
        nux  =Yi(13);  nuy  =Yi(14);  nuz  =Yi(15);
        Thetax_hat=Yi(16); Thetay_hat=Yi(17); Thetaz_hat=Yi(18);
        px_0=Yi(19);  py_0=Yi(20);  pz_0=Yi(21);

        sumx_pi_pj_hat=0; sumy_pi_pj_hat=0; sumz_pi_pj_hat=0;
        sumx_Thetai_Thetaj_hat=0; sumy_Thetai_Thetaj_hat=0; sumz_Thetai_Thetaj_hat=0;
        sumx_w=0; sumy_w=0; sumz_w=0;

        for j = neighbors(i,:)
            if j==0, continue; end
            bj = (j-1)*K; Yj = Y(bj+1:bj+K);
            xj=Yj(1); yj=Yj(2); zj=Yj(3);
            pxj_hat=Yj(7); pyj_hat=Yj(8); pzj_hat=Yj(9);
            Thetaxj_hat=Yj(16); Thetayj_hat=Yj(17); Thetazj_hat=Yj(18);

            sumx_pi_pj_hat = sumx_pi_pj_hat + (px_hat - pxj_hat);
            sumy_pi_pj_hat = sumy_pi_pj_hat + (py_hat - pyj_hat);
            sumz_pi_pj_hat = sumz_pi_pj_hat + (pz_hat - pzj_hat);

            sumx_Thetai_Thetaj_hat = sumx_Thetai_Thetaj_hat + (Thetax_hat - Thetaxj_hat);
            sumy_Thetai_Thetaj_hat = sumy_Thetai_Thetaj_hat + (Thetay_hat - Thetayj_hat);
            sumz_Thetai_Thetaj_hat = sumz_Thetai_Thetaj_hat + (Thetaz_hat - Thetazj_hat);

            norm_pij = norm([x, y, z] - [xj, yj, zj]);
            w_norm_pij = ita^2 / (norm_pij - r_s);
            if (norm_pij - r_s) <= 0, w_norm_pij = 15000; end

            sumx_w = sumx_w + (w_norm_pij) * (x - xj) / norm_pij;
            sumy_w = sumy_w + (w_norm_pij) * (y - yj) / norm_pij;
            sumz_w = sumz_w + (w_norm_pij) * (z - zj) / norm_pij;
        end

        dpx_hat = iotax; dpy_hat = iotay; dpz_hat = iotaz;
        diotax = -l_2*(iotax - nux) - l_1*(l_2^2+1)/l_2^2*(sumx_pi_pj_hat-(px_0-px_hat));
        diotay = -l_2*(iotay - nuy) - l_1*(l_2^2+1)/l_2^2*(sumy_pi_pj_hat-(py_0-py_hat));
        diotaz = -l_2*(iotaz - nuz) - l_1*(l_2^2+1)/l_2^2*(sumz_pi_pj_hat-(pz_0-pz_hat));
        dnux   = -(l_1/l_2^2) * (sumx_pi_pj_hat-(px_0-px_hat));
        dnuy   = -(l_1/l_2^2) * (sumy_pi_pj_hat-(py_0-py_hat));
        dnuz   = -(l_1/l_2^2) * (sumz_pi_pj_hat-(pz_0-pz_hat));

        ux = -zeta_1*(x-px_hat) + sumx_w + Thetax_hat - zeta_2*sumx_Thetai_Thetaj_hat;
        uy = -zeta_1*(y-py_hat) + sumy_w + Thetay_hat - zeta_2*sumy_Thetai_Thetaj_hat;
        uz = -zeta_1*(z-pz_hat) + sumz_w + Thetaz_hat - zeta_2*sumz_Thetai_Thetaj_hat;

        dThetax_hat = -zeta_3*sumx_Thetai_Thetaj_hat - zeta_4*Thetax_hat + ((zeta_1*zeta_4 - zeta_5)/zeta_1) * (-sumx_w + zeta_1*(x - px_hat));
        dThetay_hat = -zeta_3*sumy_Thetai_Thetaj_hat - zeta_4*Thetay_hat + ((zeta_1*zeta_4 - zeta_5)/zeta_1) * (-sumy_w + zeta_1*(y - py_hat));
        dThetaz_hat = -zeta_3*sumz_Thetai_Thetaj_hat - zeta_4*Thetaz_hat + ((zeta_1*zeta_4 - zeta_5)/zeta_1) * (-sumz_w + zeta_1*(z - pz_hat));

        dx = Vx; dy = Vy; dz = Vz;
        dVx = ux; dVy = uy; dVz = uz;

        if ~isfield(params,'traj_type'), params.traj_type = 'line_const'; end
        if ~isfield(params,'vxy'),  params.vxy = 45; end
        if ~isfield(params,'vz'),   params.vz  = 5;  end
        if ~isfield(params,'psi0'), params.psi0= 0;  end
        vxy  = params.vxy;  vz = params.vz;  psi0 = params.psi0;

        switch lower(params.traj_type)
            case {'line','line_const'} 
                dpx_0 = vxy * cos(psi0);    dpy_0 = vxy * sin(psi0);    dpz_0 = vz;
            case 'line_circ'
                dpx_0 = 45 * cos(t/5);
                dpy_0 = 45 * sin(t/5);
                dpz_0 = vz;
            otherwise 
                if ~isfield(params,'dpsi'), params.dpsi = deg2rad(10); end
                if ~isfield(params,'t0'),   params.t0   = 12.5;       end
                if ~isfield(params,'tau'),  params.tau  = 5.0;        end
                u = (t - params.t0) / max(params.tau, eps);
                if     u <= 0, s = 0;
                elseif u >= 1, s = 1;
                else,          s = 3*u^2 - 2*u^3;
                end
                psi   = psi0 + params.dpsi * s;
                dpx_0 = vxy * (cos(psi)+0.2);
                dpy_0 = vxy * (sin(psi)+0.2);
                dpz_0 = vz;
        end

        dY(base+1:base+K) = [dx; dy; dz; dVx; dVy; dVz; ...
                             dpx_hat; dpy_hat; dpz_hat; ...
                             diotax; diotay; diotaz; ...
                             dnux; dnuy; dnuz; ...
                             dThetax_hat; dThetay_hat; dThetaz_hat; ...
                             dpx_0; dpy_0; dpz_0];
    end
end
