clear; clc; close all;
warning('off','MATLAB:ode45:IntegrationTolNotMet');

NUM_LINE_CONST = 600;
NUM_LINE_CIRC  = 600;
NUM_CURVE      = 1200;


tspan   = [0 25];
dt_samp = 0.5;
t_fixed = tspan(1):dt_samp:tspan(2);
T = numel(t_fixed);
K = 21; 
Ntot = 12;

output_line_const = zeros(NUM_LINE_CONST, T, Ntot, K);
output_line_circ  = zeros(NUM_LINE_CIRC,  T, Ntot, K);
output_curve      = zeros(NUM_CURVE,      T, Ntot, K);
meta_group_sizes_line_const = zeros(NUM_LINE_CONST,2);
meta_group_sizes_line_circ  = zeros(NUM_LINE_CIRC, 2);
meta_group_sizes_curve      = zeros(NUM_CURVE,     2);
cnt_line_const = 0; cnt_line_circ = 0; cnt_curve = 0;

paramsA = struct('m',20,'g',9.80665,'rho',1.225,'S_w',1.37, ...
    'C_D0',0.02,'k_d',0.1,'k_n',1,'l_1',1.055,'l_2',2.35, ...
    'zeta_1',1.5,'zeta_2',0.08,'zeta_3',-0.2581333,'zeta_4',3.3, ...
    'zeta_5',0.11,'r_s',20,'r_c',100,'r_o',300,'ita',25);

paramsB = paramsA; paramsB.r_s = 20; paramsB.r_c = 100;

odeopt = odeset('InitialStep',0.05, ...
                'RelTol',1e-3, 'AbsTol',1e-6, ...
                'MaxStep',0.1, ... 
                'OutputFcn', @(t,y,flag) timeoutOutputFcn(t,y,flag, 3));
skipped_cnt = 0;  


G1_CHOICES = [3,4,5,6,7,8,9];
space_box = [0 1000;0 1000;0 400];


TOTAL_RUNS = NUM_LINE_CONST + NUM_LINE_CIRC + NUM_CURVE;

useWaitbar = usejava('desktop') && feature('ShowFigureWindows');
tStart = tic;
lastN = 0;
lastUpdate = 0;

if useWaitbar
    hwb = waitbar(0, 'Initializing...', 'Name', 'Generating UAV dataset');
else
    fprintf('Generating %d samples...\n', TOTAL_RUNS);
end


while (cnt_line_const < NUM_LINE_CONST) || ...
      (cnt_line_circ  < NUM_LINE_CIRC)  || ...
      (cnt_curve      < NUM_CURVE)

    if cnt_line_const < NUM_LINE_CONST
        mode = 'line_const';
    elseif cnt_line_circ < NUM_LINE_CIRC
        mode = 'line_circ';
    else
        mode = 'curve';
    end

    N1 = G1_CHOICES(randi(numel(G1_CHOICES)));
    N2 = Ntot - N1;
    target0 = [ ...
        randInRange(space_box(1,1), space_box(1,2)); ...
        randInRange(space_box(2,1), space_box(2,2)); ...
        randInRange(space_box(3,1), space_box(3,2))  ...
    ];
    params1 = paramsA; params2 = paramsB;

    switch mode
        case 'line_const'
            params1.traj_type = 'line_const';
            params2.traj_type = 'line_const';
            psi0 = 2*pi*rand; vxy = 40+10*rand; vz = 3+5*rand;
            params1.vxy=vxy; params2.vxy=vxy;
            params1.vz=vz;   params2.vz=vz;
            params1.psi0=psi0; params2.psi0=psi0;

        case 'line_circ'
            params1.traj_type = 'line_circ';
            params2.traj_type = 'line_circ';
            params1.vz = 5; params2.vz = 5;

        case 'curve'
            params1.traj_type = 'curve';
            params2.traj_type = 'curve';
            vxy  = 40+10*rand; vz = 3+5*rand; psi0 = 2*pi*rand;
            dpsi = deg2rad(30 + 50*rand); 
            t0_turn   = 6 + 13*rand;            
            tau  = 4 + 2*rand;            
            params1.vxy=vxy; params2.vxy=vxy; params1.vz=vz; params2.vz=vz;
            params1.psi0=psi0; params2.psi0=psi0;
            params1.dpsi=dpsi; params2.dpsi=dpsi;
            params1.t0=t0_turn;     params2.t0=t0_turn;
            params1.tau=tau;   params2.tau=tau;
    end

    center1 = target0 + sampleOnShell(params1.r_o/2, params1.r_o);
    center2 = target0 + sampleOnShell(params2.r_o/2, params2.r_o);

    Y0_1 = zeros(N1*K,1); 
    Y0_2 = zeros(N2*K,1);

    for i = 1:N1
        b=(i-1)*K;
        pos = center1 + sampleInBall(params1.r_c/2);
        vel = randSpeedVec(30,40);
        Y0_1(b+1:b+3)  = pos;                 
        Y0_1(b+4:b+6)  = vel;                  
        Y0_1(b+19:b+21)= target0;         
    end

    for i = 1:N2
        b=(i-1)*K;
        pos = center2 + sampleInBall(params2.r_c/2);
        vel = randSpeedVec(10,15);
        Y0_2(b+1:b+3)  = pos;
        Y0_2(b+4:b+6)  = vel;
        Y0_2(b+19:b+21)= target0;
    end

    neighbors_1 = zeros(N1,N1); 
    neighbors_2 = zeros(N2,N2);


    ok = true;
    try
        [T1,Y1] = ode45(@(t,Y) uavWithObserverAndControllerDynamics(t,Y,params1,neighbors_1), tspan, Y0_1, odeopt);
        [T2,Y2] = ode45(@(t,Y) uavWithObserverAndControllerDynamics(t,Y,params2,neighbors_2), tspan, Y0_2, odeopt);
    
        ok = ok & isSolutionOK(T1,Y1,t_fixed, N1*K);
        ok = ok & isSolutionOK(T2,Y2,t_fixed, N2*K);
    
    catch ME
        ok = false;
    end
    
    if ~ok
        continue;
    end

    Y1_fix = interp1(T1, Y1, t_fixed, 'linear');
    Y2_fix = interp1(T2, Y2, t_fixed, 'linear');
    Y1_fix_3d = permute(reshape(Y1_fix.',[K,N1,T]), [3,2,1]);
    Y2_fix_3d = permute(reshape(Y2_fix.',[K,N2,T]), [3,2,1]);
    Y12_full  = cat(2, Y1_fix_3d, Y2_fix_3d); 

    switch mode
        case 'line_const'
            cnt_line_const = cnt_line_const + 1;
            output_line_const(cnt_line_const,:,:,:) = Y12_full;
            meta_group_sizes_line_const(cnt_line_const,:) = [N1,N2];

        case 'line_circ'
            cnt_line_circ = cnt_line_circ + 1;
            output_line_circ(cnt_line_circ,:,:,:) = Y12_full;
            meta_group_sizes_line_circ(cnt_line_circ,:) = [N1,N2];

        case 'curve'
            cnt_curve = cnt_curve + 1;
            output_curve(cnt_curve,:,:,:) = Y12_full;
            meta_group_sizes_curve(cnt_curve,:) = [N1,N2];
    end
    
    done = cnt_line_const + cnt_line_circ + cnt_curve;

    if done > lastN && (done == TOTAL_RUNS || done - lastN >= 5 || (toc(tStart) - lastUpdate) > 0.5)
        elapsed = toc(tStart);
        rate = done / max(elapsed, eps);
        eta  = (TOTAL_RUNS - done) / max(rate, eps);
        msg  = sprintf('%.1f%% | %d/%d | elapsed %s | ETA %s', ...
                       100*done/TOTAL_RUNS, done, TOTAL_RUNS, sec2hms(elapsed), sec2hms(eta));
        if useWaitbar
            if isvalid(hwb)
                waitbar(done / TOTAL_RUNS, hwb, msg);
            end
        else

            fprintf('\r%s', msg);
        end
        lastN = done;
        lastUpdate = toc(tStart);
    end


end

output = cat(1, output_line_const, output_line_circ, output_curve); 
meta_group_sizes = cat(1, meta_group_sizes_line_const, meta_group_sizes_line_circ, meta_group_sizes_curve);

meta_traj_code = [ones(NUM_LINE_CONST,1); 2*ones(NUM_LINE_CIRC,1); zeros(NUM_CURVE,1)];

save('matrix12_2400_fullseq.mat', 'output', 'K', 'T', 'dt_samp', ...
     'meta_group_sizes', 'meta_traj_code', '-v7.3');

fprintf('Done. line_const=%d, line_circ=%d, curve=%d, total=%d, T=%d\n', ...
        cnt_line_const, cnt_line_circ, cnt_curve, size(output,1), T);


function x = randInRange(a,b)
    x = a + (b-a) * rand;
end

function v = sampleOnShell(rmin, rmax)
    r = rmin + (rmax - rmin) * rand;    
    u = randn(3,1); u = u / max(norm(u), eps);  
    v = r * u;
end

function v = sampleInBall(R)
    u = randn(3,1); u = u / max(norm(u), eps);  
    rho = R * (rand)^(1/3);
    v = rho * u;
end

function v = randSpeedVec(vmin, vmax)
    speed = vmin + (vmax - vmin) * rand;
    u = randn(3,1); u = u / max(norm(u), eps);
    v = speed * u;
end

function status = timeoutOutputFcn(~, ~, flag, maxTimeSec)

    persistent tStart
    switch flag
        case 'init'
            tStart = tic;  status = 0;
        case 'done'
            status = 0;
        otherwise
            if toc(tStart) > maxTimeSec
                status = 1;
            else
                status = 0;
            end
    end
end

function ok = isSolutionOK(T, Y, t_fixed, K)

    ok = true;
    if isempty(T) || size(Y,1) ~= numel(T) || size(Y,2) ~= K
        ok = false; return;
    end

    if T(1) > t_fixed(1) || T(end) < t_fixed(end)
        ok = false; return;
    end
    if any(~isfinite(Y(:))) || any(isnan(Y(:)))
        ok = false; return;
    end
    if any(abs(Y(:)) > 1e6)
        ok = false; return;
    end
end


function s = sec2hms(x)

    if ~isscalar(x) || ~isfinite(x) || x < 0
        s = '--:--:--'; return;
    end
    x = floor(x);
    h = floor(x/3600);
    m = floor(mod(x,3600)/60);
    ssec = mod(x,60);
    s = sprintf('%02d:%02d:%02d', h, m, ssec);
end



