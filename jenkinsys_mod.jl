using juliajim.HARMONIC
using LinearAlgebra

# * Harmonic Coefficients
function Fnl_n(u, p)
    (; zt, w0, kt, fs, f, Om, typ, order) = p;

    d0sig = (typ==:MMS) ? w0 : (typ==:EMS) ? Om : (typ==:EMS0) ? sqrt(w0^2+kt) : 0;

    taust = (kt*norm(u)<fs) ? π : acos(1-2fs/kt/norm(u));
    dtaust_da0 = (kt*norm(u)<fs) ? 0 : -sqrt(fs/(norm(u)*kt-fs))/norm(u);

    if order == 1
        Fnl = kt/2π*[2taust-sin(2taust) 1-cos(2taust);
                     -(1-cos(2taust)) 2taust-sin(2taust)]*u;
        dFnl = kt/2π*[2taust-sin(2taust) 1-cos(2taust);
                      -(1-cos(2taust)) 2taust-sin(2taust)] +
               dtaust_da0*(kt/π*[1-cos(2taust) sin(2taust);
                                 -sin(2taust) 1-cos(2taust)]*u).*u'/norm(u);
    elseif order == 2
        Fnl = (kt/π/d0sig/24/√30)^2*[-(17*cos(6*taust)+308*cos(4*taust)-185*cos(2*taust)-140) 5*(sin(6*taust)+26*sin(4*taust)-55*sin(2*taust));
	                             -(5*(sin(6*taust)+26*sin(4*taust)-55*sin(2*taust))) -(17*cos(6*taust)+308*cos(4*taust)-185*cos(2*taust)-140)]*u;
        dFnl = (kt/π/d0sig/24/√30)^2*([-(17*cos(6*taust)+308*cos(4*taust)-185*cos(2*taust)-140) 5*(sin(6*taust)+26*sin(4*taust)-55*sin(2*taust));
	                               -(5*(sin(6*taust)+26*sin(4*taust)-55*sin(2*taust))) -(17*cos(6*taust)+308*cos(4*taust)-185*cos(2*taust)-140)] +
		                      dtaust_da0*([2*(51*sin(6*taust)+616*sin(4*taust)-185*sin(2*taust)) 10*(3*cos(6*taust)+52*cos(4*taust)-55*cos(2*taust));
	                                           -(10*(3*cos(6*taust)+52*cos(4*taust)-55*cos(2*taust))) 2*(51*sin(6*taust)+616*sin(4*taust)-185*sin(2*taust))]*u).*u'/norm(u)
        );
    else
        error("Order > 2 not implemented.")
    end

    return Fnl, dFnl;
end

# * Slow Flow
function sflow!(du, u, p, t=0.0)
    (; zt, w0, kt, fs, f, Om, typ, order) = p
    d0sig = (typ==:MMS) ? w0 : (typ==:EMS) ? Om : (typ==:EMS0) ? sqrt(w0^2+kt) : 0;
    
    if typ == :MMS
        d1sig = Om-w0;
        if order == 1
            Fnl, dFnl = Fnl_n(u, p);
            du[:] = [-zt*w0 -d1sig; d1sig -zt*w0]*u + 1/2w0*[0 1;-1 0]*Fnl +
                    f/2w0*[0;1];
            J = [-zt*w0 -d1sig; d1sig -zt*w0] + 1/2w0*[0 1;-1 0]*dFnl;
        elseif order == 2
            d1u_du = sflow!(du, u, (;p..., order=1), t);
            d1_u = copy(du);
            d1sq_u = d1u_du*d1_u;
            Fnl, dFnl = Fnl_n(u, p);
            
            d2_u = 1/2w0*[0 1;-1 0]*d1sq_u + 1/w0*[-d1sig zt*w0;-zt*w0 -d1sig]*d1_u +
                   1/2w0*[-2zt*w0*d1sig 0;0 -2zt*w0*d1sig]*u +
                   1/2w0*[0 1;-1 0]*Fnl;
            # Incomplete Gradient (use sparingly)
            d2u_du = 1/w0*[-d1sig zt*w0;-zt*w0 -d1sig]*d1u_du +
                     1/2w0*[-2zt*w0*d1sig 0;0 -2zt*w0*d1sig] +
                     1/2w0*[0 1;-1 0]*dFnl;

            du[:] = d1_u + d2_u;
            J = d1u_du + d2u_du;
        else
            error("Order > 2 not implemented");
        end
    elseif typ == :EMS
        if order == 1
            Fnl, dFnl = Fnl_n(u, p);
            du[:] = [-zt*w0 -(Om^2-w0^2)/2Om;(Om^2-w0^2)/2Om -zt*w0]*u +
                    1/2Om*[0 1;-1 0]*Fnl + f/2Om*[0;1];
            J = [-zt*w0 -(Om^2-w0^2)/2Om;(Om^2-w0^2)/2Om -zt*w0] +
                1/2Om*[0 1;-1 0]*dFnl;
        elseif order == 2
            d1u_du = sflow!(du, u, (;p..., order=1), t);
            d1_u = copy(du);
            d1sq_u = d1u_du*d1_u;
            Fnl, dFnl = Fnl_n(u, p);

            d2_u = 1/2Om*[0 1;-1 0]*d1sq_u + 1/Om*[0 zt*w0;-zt*w0 0]*d1_u +
                   1/2Om*[0 1;-1 0]*Fnl;
            # Incomplete Gradient (use sparingly)
            d2u_du = 1/Om*[0 zt*w0;-zt*w0 0]*d1u_du + 1/2Om*[0 1;-1 0]*dFnl;

            du[:] = d1_u + d2_u;
            J = d1u_du + d2u_du;
        else
            error("Order > 2 not implemented");
        end
    elseif typ == :EMS0
        omb = sqrt(w0^2+kt);
        d1sig = Om-omb;
        if order == 1
            Fnl, dFnl = Fnl_n(u, p);
            du[:] = [-zt*w0 -(omb^2-w0^2)/2omb-d1sig;(omb^2-w0^2)/2omb+d1sig -zt*w0]*u +
                    1/2omb*[0 1;-1 0]*Fnl + f/2omb*[0;1];
            J = [-zt*w0 -(omb^2-w0^2)/2omb-d1sig;(omb^2-w0^2)/2omb+d1sig -zt*w0] +
                1/2omb*[0 1;-1 0]*dFnl;
        elseif order == 2
            d1u_du = sflow!(du, u, (;p..., order=1), t);
            d1_u = copy(du);
            d1sq_u = d1u_du*d1_u;
            Fnl, dFnl = Fnl_n(u, p);

            d2_u = 1/omb*[0 1;-1 0]*d1sq_u + 1/omb*[-d1sig zt*w0;-zt*w0 -d1sig]*d1_u +
                   1/2omb*[-2zt*w0*d1sig -d1sig^2;d1sig^2 -2zt*w0*d1sig]*u +
                   1/2omb*[0 1;-1 0]*Fnl;
            # Incomplete Jacobian (use sparingly)
            d2u_du = 1/omb*[-d1sig zt*w0;-zt*w0 -d1sig]*d1u_du +
                     1/2omb*[-2zt*w0*d1sig -d1sig^2;d1sig^2 -2zt*w0*d1sig] +
                     1/2omb*[0 1;-1 0]*dFnl;

            du[:] = d1_u + d2_u;
            J = d1u_du + d2u_du;
        else
            error("Order > 2 not implemented");
        end
    else
        error("Unrecognized method: ", typ);
    end
    return J;
end

function sflow_harms(u, p)
    (; zt, w0, kt, fs, f, Om, typ, order) = p

    d0sig = (typ==:MMS) ? w0 : (typ==:EMS) ? Om : (typ==:EMS0) ? sqrt(w0^2+kt) : 0;

    taust = (kt*norm(u)<fs) ? π : acos(1-2fs/kt/norm(u));
    Fnl3 = im*kt/12π*(1-(2exp(im*2taust)-1)*exp(im*2taust))*norm(u)*exp(-3im*atan(-u[2],u[1]));
    Fnl5 = im*kt/60π*(1-(3exp(im*2taust)-2)*exp(im*4taust))*norm(u)*exp(-5im*atan(-u[2],u[1]));

    if kt*norm(u)<fs
        Fnl3 *= 0;
        Fnl5 *= 0;
    end


    A3 = Fnl3/8d0sig^2;
    A5 = Fnl5/24d0sig^2;

    return A3, A5
end

# * HB Residue
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

    if du !== nothing
        du[:] = E*u+Fnl;
        du[rinds[1]] -= f;        
    end

    return Fnl;
end
