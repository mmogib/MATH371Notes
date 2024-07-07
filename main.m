f = @(x) x^2 - 2;
a = 0;
b = 2;
TOL = 1e-7;
N0 = 1000;

[p, TT, FLAG] = bisect(f, a, b, TOL, N0);
if FLAG
    disp('Approximate root:');
    disp(p);
    disp('Iteration matrix:');
    disp(TT);
else
    disp('Method failed after maximum number of iterations');
end