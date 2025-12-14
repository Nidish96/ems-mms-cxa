using juliajim.HARMONIC

# * Setup System
function rocfun!(du, u, p, t)
    (; zt, w0, al, f, Om) = p;
    
    du[1] = u[2];
    du[2] = -w0^2*u[1]-2zt*w0*u[2] - al*u[1]^3 + f*cos(Om*t);
    return du;
end
  
# * Slow Flow
function sflow!(du, u, p, t=0.0)
    (; zt, w0, kt, fs, f, Om, typ) = p;

    taust = (kt*norm(u)<fs) ? π : acos(1-2fs/kt/norm(u));
    Fnl = kt/2π*[2taust-sin(2taust) 1-cos(2taust);
                 -(1-cos(2taust)) 2taust-sin(2taust)]*u;

    if typ==:MMS
        du[:] = [-zt*w0 -(Om-w0);(Om-w0) -zt*w0]*u +
                1/2w0*[0 1;-1 0]*Fnl + [0;f/2w0];
    elseif typ==:EMS
        du[:] = [-zt*w0 -(Om^2-w0^2)/2Om;(Om^2-w0^2)/2Om -zt*w0]*u +
                1/2Om*[0 1;-1 0]*Fnl + [0;f/2Om];
    elseif typ==:EMS0 
        wa = sqrt(w0^2+kt);
        du[:] = [-zt*w0 -(wa^2-w0^2)/2wa-(Om-wa);(wa^2-w0^2)/2wa+(Om-wa) -zt*w0]*u +
                1/2wa*[0 1;-1 0]*Fnl + [0;f/2wa];
    end

    return du;
end

function sflow_harms(u, p)
    (; zt, w0, kt, fs, f, Om, typ) = p;

    taust = (kt*norm(u)<fs) ? π : acos(1-2fs/kt/norm(u));
    Fnl3 = -(im*norm(u)*kt*exp(-3im*atan(-u[2],u[1]))*(exp(im*2taust)-1)*(2exp(2im*taust)+1))/12π;

    if typ==:MMS
        A3 = Fnl3/8w0^2;
        A5 = 0im;
    elseif typ==:EMS
        A3 = Fnl3/8Om^2;
        A5 = 0;
    elseif typ==:EMS0
        wa = sqrt(w0^2+kt);
        A3 = Fnl3/8wa^2;
        A5 = 0;        
    end
    return A3, A5
end


# * Harmonic Balance

function hbresfun!(du, u, p, h, N)
    (; zt, w0, kt, fs, f, Om) = p;
    
    E, _ = HARMONICSTIFFNESS(1.0, 2zt*w0, w0^2, Om, h);

    ut = AFT(u, h,N, :f2t);
    f0 = 0;
    ft = kt*ut;
    it = 0;
    while abs((ft[end]-f0)/fs)>eps()^(4//5) || it<1
        f0 = ft[end];
        slip = false;
        for ti in 1:N
            tim1 = mod(ti-1-1,N)+1;
            fsp = kt*(ut[ti]-ut[tim1]) + ft[tim1];
            ft[ti] = (abs(fsp)<fs) ? fsp : fs*sign(fsp);
            slip = (abs(fsp)<fs) ? slip : true;
        end
        if !slip
            ft = ft.-sum(ft)/N;
        end
        it += 1;
    end
    
    Fnl = AFT(ft, h,N, :t2f);

    _, _, _, rinds, _ = HINDS(1, h);

    du[:] = E*u+Fnl;
    du[rinds[1]] -= f;

    return du;
end
