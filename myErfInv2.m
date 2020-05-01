

% // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% // this routine does the bookkeeping for an infection event
% //Quick & dirty, tolerance under +-6e-3.
% // Work based on "A handy approximation for the error function and its inverse" by Sergei Winitzki.
% // https:% //stackoverflow.com/questions/631629/erfx-and-math-h
function res = myErfInv2(x)
if x < 0
    sgn =  -1.0;
else
    sgn = 1.0;
end
x = -x*x + 1.0;        % // x = 1 - x*x;
lnx = log(x);

tt1 = 2.0/(M_PI*0.147) + 0.5 * lnx;
tt2 = lnx/(0.147) ;

res = (sgn*sqrt(-tt1 + sqrt(tt1*tt1 - tt2)));
end
