# * Analytical Forced Response Solution
ansol(p) = p.f/(p.w0^2-p.Om^2+2im*p.zt*p.w0*p.Om);

# * Setup Linear System
function rocfun!(du, u, p, t)
    (; zt, w0, f, Om) = p;
    
    du[1] = u[2];
    du[2] = -w0^2*u[1]-2zt*w0*u[2] + f*cos(Om*t);
    return du;
end
  
# * Slow Flow
function sflow!(du, u, p, t=0.0)
    (; zt, w0, f, Om, typ, order) = p;

    if typ==:MMS
        if order==1
            du[1] = -zt*w0*u[1] - (Om-w0)u[2];
            du[2] = -zt*w0*u[2] + (Om-w0)u[1] + f/2w0;
        elseif order>=2  # Nothing changes for higher orders when ω2, ω3, … = 0.
            sflow!(du, u, (;p..., order=1), t);
            du[1] += -(zt*(2u[2]*zt*w0^2-f))/4w0;
            du[2] += -(f*(Om-w0)-2u[1]*zt^2*w0^3)/4w0^2;
        elseif order<=-2
            sflow!(du, u, (;p..., order=1), t);
            du[1] += -((u[2]*(Om-w0)^2)/2w0)-u[1]*zt*(Om-w0);
            du[2] += (u[1]*(Om-w0)^2)/2w0-u[2]*zt*(Om-w0);
        end
    elseif typ==:EMS
        if order==1
            du[1] = -zt*w0*u[1] - (Om^2-w0^2)/2Om*u[2];
            du[2] = -zt*w0*u[2] + (Om^2-w0^2)/2Om*u[1] + f/2Om;
        elseif order==2
            sflow!(du, u, (;p..., order=1), t);
            du[1] += f*zt*w0/4Om^2 - (w0^4+(2zt^2-1)*2Om^2*w0^2+Om^4)/8Om^3*u[2];
            du[2] += f*(Om^2-w0^2)/8Om^3 + (w0^4+(2zt^2-1)*2Om^2*w0^2+Om^4)/8Om^3*u[1];
        elseif order==3
            sflow!(du, u, (;p..., order=2), t);
            du[1] += (u[2]*w0^6+Om^2*(4u[2]*zt^2-3u[2])*w0^4-2f*Om*zt*w0^3+Om^4*(3u[2]-4u[2]*zt^2)*w0^2+2f*Om^3zt*w0-u[2]*Om^6)/16Om^5;
            du[2] += (-u[1]*w0^6+Om^2*((3u[1]-4u[1]*zt^2)*w0^4-2f*w0^2)+f*w0^4+Om^4*((4u[1]*zt^2-3u[1])*w0^2+f)+u[1]*Om^6)/16Om^5;
        elseif order==4
            sflow!(du, u, (;p..., order=3), t);
            du[1] += f*zt*w0*(Om^2-w0^2)/8Om^4 + (w0^6+(4Om^2*zt^2-3Om^2)*w0^4+(-4Om^4*zt^2+3Om^4)*w0^2-Om^6)/16Om^5*u[2];
            du[2] += f*(Om^2-w0^2)^2/16Om^5 - (w0^6+(4Om^2*zt^2-3Om^2)*w0^4+(-4Om^4*zt^2+3Om^4)*w0^2-Om^6)/16Om^5*u[1];
        end
    end

    return du;
end
