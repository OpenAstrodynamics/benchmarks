function [v, v0] = lambert(k, r0, r, tof, short, numiter, rtol)
if short == true
    t_m = 1;
else
    t_m = -1;
end

norm_r0 = norm(r0);
norm_r = norm(r);
cos_dnu = dot(r0, r) / (norm_r0 * norm_r);

A = t_m * sqrt(norm_r * norm_r0 * (1 + cos_dnu));

if A == 0
    error('Cannot compute orbit, phase angle is 180 degrees');
end

psi = 0;
psi_low = -4 * pi;
psi_up = 4 * pi;

count = 0;
converged = false;
while count < numiter
    y = norm_r0 + norm_r + A * (psi * c3(psi) - 1) / sqrt(c2(psi));
    if A > 0 && y < 0
        while y < 0
            psi_low = psi;
            psi = (0.8 * (1.0 / c3(psi)) * ...
                (1.0 - (norm_r0 + norm_r) * sqrt(c2(psi)) / A));
            y = norm_r0 + norm_r + A * (psi * c3(psi) - 1) / sqrt(c2(psi));
        end
    end
    
    xi = sqrt(y / c2(psi));
    tof_new = (xi^3 * c3(psi) + A * sqrt(y)) / sqrt(k);
    
    if abs((tof_new - tof) / tof) < rtol
        converged = true;
        break
    else
        count = count + 1;
        if tof_new <= tof
            psi_low = psi;
        else
            psi_up = psi;
        end
        psi = (psi_up + psi_low) / 2;
    end
end

if ~converged
    error('Maximum number of iterations reached');
end

f = 1 - y / norm_r0;
g = A * sqrt(y / k);
gdot = 1 - y / norm_r;

v0 = (r - f * r0) / g;
v = (gdot * r - r0) / g;