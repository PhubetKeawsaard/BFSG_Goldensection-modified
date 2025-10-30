function [xmin,fmin,Xk,Fk,Gk,Lk,nF,nG,nH,CHN,IFLAG,k] = BFGS(Fcn,x0,epsilon,mu,eta,itmax)
% Fcn(x,1)->f(x), Fcn(x,2)->g(x)
    x = x0(:);
    IFLAG = 1;
    % histories / counters
    Xk = x.'; Fk = []; Gk = []; CHN = [];Lk=[];
    nF = 0; nG = 0; nH = 0;
    % initial f,g
    f = Fcn(x,1); nF = nF+1;
    g = Fcn(x,2); nG = nG+1;
    n = numel(x);
    H = eye(n); I = eye(n);
    % line-search params
    for k = 1:itmax
        % --- direction ---
        sdir = -H * g;
        if g' * sdir >= 0
            sdir = -g; CHN = [CHN; 1];   % fallback
        else
            CHN = [CHN; 0];
        end
        % --- line search ---
        fund1d = @(a) Fcn(x+a*sdir,1) ;
        g1d =@(a) Fcn(x+a*sdir,2) ;
        [alpha, ~] = Goldensection(0, 1,fund1d,g1d, sdir, 1e-8,60,mu, eta);
        if ~(alpha == 1)
            nF=60+nF ;
        end
        Lk = [Lk;alpha] ;
        % --- update ---
        x_new = x + alpha * sdir;
        f_new = Fcn(x_new,1); nF = nF+1;
        g_new = Fcn(x_new,2); nG = nG+1;
        % histories
        Xk = [Xk; x_new.'];
        Fk = [Fk; f_new];
        Gk = [Gk; (g_new).'];
        % --- stopping criteria ---
        dx_norm    = norm(x_new - x);
        x_scale    = max(1, norm(x));
        step_small = dx_norm <= mu * x_scale;
        df         = abs(f_new - f);
        f_scale    = max(1, abs(f));
        f_small    = df <= eta * f_scale;
        g_small    = norm(g_new) <= epsilon;
        if g_small || (step_small && f_small)
            xmin = x_new; fmin = f_new; IFLAG = 0; return;
        end
        % --- BFGS (inverse) update ---
        s = x_new - x;
        y = g_new - g;
        sy = s' * y;
        if sy > 1e-12 * norm(s) * norm(y)
            rho = 1 / sy;
            V = I - rho * (s * y.');
            H = V * H * V.' + rho * (s * s.');
        end
        x = x_new; f = f_new; g = g_new;
    end
    % max iters
    xmin = x; fmin = f;
end
