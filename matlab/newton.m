function p = newton(x0, func, deriv, varargin)

switch nargin
    case 3
        maxiter = 50;
        tol = 1e-8;
    case 4
        maxiter = varargin{1};
        tol = 1e-8;
    case 5
        maxiter = varargin{1};
        tol = varargin{2};
end

p0 = x0;
for ii = 1:maxiter
    p = p0 - func(p0)/deriv(p0);
    if abs(p - p0) < tol
        return
    end
    p0 = p;
end
error('Not converged.');
