using juliajim.HARMONIC
using juliajim.MDOFUTILS

# * Setup System
function rocfun!(du, u, p, t)
    (; c, w0, al, f, Om) = p;
    
    du[1] = u[2];
    du[2] = -w0^2*u[1]+c*u[2] - al*u[1]^2*u[2] + f*cos(Om*t);
    return du;
end
  
# * Slow Flow
function sflow!(du, u, p, t=0.0)
    (; c, w0, al, f, Om, typ, order) = p;

    if typ==:MMS
        if order==1
            du[1] = -(u[2]*(Om-w0))-(u[1]*u[2]^2*al)/8-(u[1]^3*al)/8+(u[1]*c)/2;
            du[2] = u[1]*(Om-w0)+f/2w0-(u[2]^3*al)/8-(u[1]^2*u[2]*al)/8+(u[2]*c)/2;
        elseif order==2
            sflow!(du, u, (;p..., order=1), t);
            du[1] += -(7u[2]^5*al^2)/256w0-(7u[1]^2*u[2]^3*al^2)/128w0-(7u[1]^4*u[2]*al^2)/256w0+(u[2]^3*c*al)/8w0+(u[1]^2*u[2]*c*al)/8w0-(u[2]*c^2)/8w0+(3*u[2]^2*f*al)/32w0^2+(u[1]^2*f*al)/32w0^2-(f*c)/8w0^2;
            du[2] += -(f*(Om-w0))/4w0^2+(7u[1]*u[2]^4*al^2)/256w0+(7u[1]^3*u[2]^2*al^2)/128w0+(7u[1]^5*al^2)/256w0-(u[1]*u[2]^2*c*al)/8w0-(u[1]^3*c*al)/8w0+(u[1]*c^2)/8w0-(u[1]*u[2]*f*al)/16w0^2;
        elseif order==3
            sflow!(du, u, (;p..., order=2), t);
            du[1] += (19u[1]*u[2]^6*al^3)/8192w0^2+(57u[1]^3*u[2]^4*al^3)/8192w0^2+(57u[1]^5*u[2]^2*al^3)/8192w0^2+(19u[1]^7*al^3)/8192w0^2-(9u[1]*u[2]^4*c*al^2)/1024w0^2-(9u[1]^3*u[2]^2*c*al^2)/512w0^2-(9u[1]^5*c*al^2)/1024w0^2-(u[1]*u[2]^3*f*al^2)/128w0^3-(u[1]^3*u[2]*f*al^2)/128w0^3;
            du[2] += (19u[2]^7*al^3)/8192w0^2+(57u[1]^2*u[2]^5*al^3)/8192w0^2+(57u[1]^4*u[2]^3*al^3)/8192w0^2+(19u[1]^6*u[2]*al^3)/8192w0^2-(9u[2]^5*c*al^2)/1024w0^2-(9u[1]^2*u[2]^3*c*al^2)/512w0^2-(9u[1]^4*u[2]*c*al^2)/1024w0^2-(15u[2]^4*f*al^2)/2048w0^3-(7u[1]^2*u[2]^2*f*al^2)/1024w0^3+(u[1]^4*f*al^2)/2048w0^3;
        end
    elseif typ==:EMS
        if order==1
            du[1] = (u[2]*w0^2)/2Om-(u[1]*u[2]^2*al)/8-(u[1]^3*al)/8-(u[2]*Om)/2+(u[1]*c)/2;
            du[2] = -(u[1]*w0^2)/2Om-(u[2]^3*al)/8-(u[1]^2*u[2]*al)/8+(u[1]*Om)/2+f/2Om+(u[2]*c)/2;
        elseif order==2
            sflow!(du, u, (;p..., order=1), t);
            du[1] += -(u[2]*w0^4)/8Om^3+(u[2]*w0^2)/4Om-(7u[2]^5*al^2)/256Om-(7u[1]^2*u[2]^3*al^2)/128Om-(7u[1]^4*u[2]*al^2)/256Om+(u[2]^3*c*al)/8Om+(u[1]^2*u[2]*c*al)/8Om+(3*u[2]^2*f*al)/32Om^2+(u[1]^2*f*al)/32Om^2-(u[2]*Om)/8-(u[2]*c^2)/8Om-(f*c)/8Om^2;
            du[2] += (u[1]*w0^4)/8Om^3-(u[1]*w0^2)/4Om-(f*w0^2)/8Om^3+(7u[1]*u[2]^4*al^2)/256Om+(7u[1]^3*u[2]^2*al^2)/128Om+(7u[1]^5*al^2)/256Om-(u[1]*u[2]^2*c*al)/8Om-(u[1]^3*c*al)/8Om-(u[1]*u[2]*f*al)/16Om^2+(u[1]*Om)/8+(u[1]*c^2)/8Om+f/8Om;
        elseif order==3
            sflow!(du, u, (;p..., order=2), t);

            du[1] += (u[2]*w0^6)/16Om^5-(3u[2]*w0^4)/16Om^3+(7u[2]^5*al^2*w0^2)/512Om^3+(7u[1]^2*u[2]^3*al^2*w0^2)/256Om^3+(7u[1]^4*u[2]*al^2*w0^2)/512Om^3-(u[2]^3*c*al*w0^2)/16Om^3-(u[1]^2*u[2]*c*al*w0^2)/16Om^3-(3u[2]^2*f*al*w0^2)/64Om^4-(u[1]^2*f*al*w0^2)/64Om^4+(3u[2]*w0^2)/16Om+(u[2]*c^2*w0^2)/16Om^3+(f*c*w0^2)/16Om^4+(19u[1]*u[2]^6*al^3)/8192Om^2+(57u[1]^3*u[2]^4*al^3)/8192Om^2+(57u[1]^5*u[2]^2*al^3)/8192Om^2+(19u[1]^7*al^3)/8192Om^2-(7u[2]^5*al^2)/512Om-(7u[1]^2*u[2]^3*al^2)/256Om-(7u[1]^4*u[2]*al^2)/512Om-(9u[1]*u[2]^4*c*al^2)/1024Om^2-(9u[1]^3*u[2]^2*c*al^2)/512Om^2-(9u[1]^5*c*al^2)/1024Om^2-(u[1]*u[2]^3*f*al^2)/128Om^3-(u[1]^3*u[2]*f*al^2)/128Om^3+(u[2]^3*c*al)/16Om+(u[1]^2*u[2]*c*al)/16Om+(3u[2]^2*f*al)/64Om^2+(u[1]^2*f*al)/64Om^2-(u[2]*Om)/16-(u[2]*c^2)/16Om-(f*c)/16Om^2;
            du[2] += -((u[1]*w0^6)/16Om^5)+(3u[1]*w0^4)/16Om^3+(f*w0^4)/16Om^5-(7u[1]*u[2]^4*al^2*w0^2)/512Om^3-(7u[1]^3*u[2]^2*al^2*w0^2)/256Om^3-(7u[1]^5*al^2*w0^2)/512Om^3+(u[1]*u[2]^2*c*al*w0^2)/16Om^3+(u[1]^3*c*al*w0^2)/16Om^3+(u[1]*u[2]*f*al*w0^2)/32Om^4-(3u[1]*w0^2)/16Om-(u[1]*c^2*w0^2)/16Om^3-(f*w0^2)/8Om^3+(19u[2]^7*al^3)/8192Om^2+(57u[1]^2*u[2]^5*al^3)/8192Om^2+(57u[1]^4*u[2]^3*al^3)/8192Om^2+(19u[1]^6*u[2]*al^3)/8192Om^2+(7u[1]*u[2]^4*al^2)/512Om+(7u[1]^3*u[2]^2*al^2)/256Om+(7u[1]^5*al^2)/512Om-(9u[2]^5*c*al^2)/1024Om^2-(9u[1]^2*u[2]^3*c*al^2)/512Om^2-(9u[1]^4*u[2]*c*al^2)/1024Om^2-(15u[2]^4*f*al^2)/2048Om^3-(7u[1]^2*u[2]^2*f*al^2)/1024Om^3+(u[1]^4*f*al^2)/2048Om^3-(u[1]*u[2]^2*c*al)/16Om-(u[1]^3*c*al)/16Om-(u[1]*u[2]*f*al)/32Om^2+(u[1]*Om)/16+(u[1]*c^2)/16Om+f/16Om;
        end
    end

    return du;
end

function sflow_harms(u, p)
    (; c, w0, al, f, Om, typ, order) = p;

    if typ==:MMS
        if order==1
            A3 = -(u[2]^3-3im*u[1]*u[2]^2-3u[1]^2*u[2]+im*u[1]^3)*al/32w0;
            A5 = 0im;
        elseif order==2
            A3, A5 = sflow_harms(u, (;p..., order=1));
            A3 += (im*u[2]^5*al^2)/1024w0^2+(3u[1]*u[2]^4*al^2)/1024w0^2-(im*u[1]^2*u[2]^3*al^2)/512w0^2+(u[1]^3*u[2]^2*al^2)/512w0^2-(3im*u[1]^4*u[2]*al^2)/1024w0^2-(u[1]^5*al^2)/1024w0^2+(im*u[2]^3*c*al)/128w0^2+(3u[1]*u[2]^2*c*al)/128w0^2-(3im*u[1]^2*u[2]*c*al)/128w0^2-(u[1]^3*c*al)/128w0^2+(5im*u[2]^2*f*al)/256w0^3+(5u[1]*u[2]*f*al)/128w0^3-(5im*u[1]^2*f*al)/256w0^3
            A5 += -5al^2*(im*u[2]^5+5u[1]*u[2]^4-10im*u[1]^2*u[2]^3-10u[1]^3*u[2]^2+5im*u[1]^4*u[2]+u[1]^5)/3072w0^2;
        else
            @warn "Higher harmonics not typed up for order>2."
            return sflow_harms(u, (;p..., order=2));
        end
    elseif typ==:EMS
        if order==1
            A3 = -(u[2]^3-3im*u[1]*u[2]^2-3u[1]^2*u[2]+im*u[1]^3)*al/32Om;
            A5 = 0im;
        elseif order==2
            A3, A5 = sflow_harms(u, (;p..., order=1));
            A3 += (u[2]^3*al*w0^2)/64Om^3-(3im*u[1]*u[2]^2*al*w0^2)/64Om^3-(3u[1]^2*u[2]*al*w0^2)/64Om^3+(im*u[1]^3*al*w0^2)/64Om^3+(im*u[2]^5*al^2)/1024Om^2+(3u[1]*u[2]^4*al^2)/1024Om^2-(im*u[1]^2*u[2]^3*al^2)/512Om^2+(u[1]^3*u[2]^2*al^2)/512Om^2-(3im*u[1]^4*u[2]*al^2)/1024Om^2-(u[1]^5*al^2)/1024Om^2-(u[2]^3*al)/64Om+(3im*u[1]*u[2]^2*al)/64Om+(3u[1]^2*u[2]*al)/64Om-(im*u[1]^3*al)/64Om+(im*u[2]^3*c*al)/128Om^2+(3u[1]*u[2]^2*c*al)/128Om^2-(3im*u[1]^2*u[2]*c*al)/128Om^2-(u[1]^3*c*al)/128Om^2+(5im*u[2]^2*f*al)/256Om^3+(5u[1]*u[2]*f*al)/128Om^3-(5im*u[1]^2*f*al)/256Om^3;
            A5 += -5al^2*(im*u[2]^5+5u[1]*u[2]^4-10im*u[1]^2*u[2]^3-10u[1]^3*u[2]^2+5im*u[1]^4*u[2]+u[1]^5)/3072Om^2;
        else
            @warn "Higher harmonics not typed up for order>2."
            return sflow_harms(u, (;p..., order=2));
        end
    end
    
    return A3, A5
end


# * Harmonic Balance

function hbresfun!(du, u, p, h, N)
    (; c, w0, al, f, Om) = p;
    
    E, _ = HARMONICSTIFFNESS(1.0, -c, w0^2, Om, h);
    D1, _ = HARMONICSTIFFNESS(0, 1.0, 0.0, Om, h);
    
    ut = AFT(u, h,N, :f2t);
    udt = AFT(D1*u, h,N, :f2t);    
    ft = al*ut.^2 .*udt;
    Fnl = AFT(ft, h,N, :t2f);

    _, _, _, rinds, _ = HINDS(1, h);

    du[:] = E*u+Fnl;
    du[rinds[1]] -= f;

    return du;
end

# * HB Residue of Slow Flow
function hbslow!(R, Uw, p, h, N; U0=nothing)    
    Nhc = NHC(h);

    Ut = AFT(reshape(Uw[1:end-1], 2,Nhc)', h,N, :f2t);  # N, 2
    Rt = zeros(typeof(Ut[1,1]*p.Om), N,2);
    for i in 1:N
        Rt[i, :] = sflow!(Rt[i,:], Ut[i,:], p);
    end

    D1, _ = HARMONICSTIFFNESS(zeros(2,2), I(2), zeros(2,2), Uw[end], h);
    R[:] = D1*Uw[1:end-1]-AFT(Rt, h,N, :t2f)'[:];

    if !(U0 === nothing)
        if eltype(U0) <: Vector
            R[:] /= sum(norm.([Uw].-U0).^(-2));
        else
            R[:] /= norm(Uw-U0)^2;
        end
    end

    return R;
end
