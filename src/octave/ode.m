1;
global me = 'ode.m';

function err(fmt, varargin)
  global me
  fmt = ["[%s] " fmt];
  arg = {me, varargin{:}};
  error(fmt, arg{:});
end

function write_vector(file, r, d)
  X = 1; Y = 2; Z = 3;
  if nargin != 3; err("write_vector: nargin = %d", nargin); end
  if size(r, 1) != size(d, 1); err("write_vector: r and d missmatch"); end
  f = fopen(file, "w");
  if f == -1; err("fail to open '%s'", file); end
  n = size(r, 1)
  for i = 1:n
    r0 = r(i, :); r1 = r0 + d(i, :);
    fprintf(f, "%g %g %g\n", r0(X), r0(Y), r0(Z));
    fprintf(f, "%g %g %g\n", r1(X), r1(Y), r1(Z));
    fprintf(f, "\n\n");
  end
  fclose(f);
end

function dq = f(q, param)
  if nargin != 2; err("f: nargin = %d", nargin); end
  if numel(q) != 3; err("f: numel(q)=%d", numel(q)); end
  if numel(param) != 3; err("f: numel(param)=%d", numel(param)); end

  b = param(1); sg = param(2); r = param(3);
  x = q(1); y = q(2); z = q(3);

  dx = sg*(y - x);
  dy = -x*z + r*x - y;
  dz = x*y - b*z;

  dq = zeros(size(q));
  dq(1) = dx; dq(2) = dy; dq(3) = dz;
end

function dp = df(p, q, param)
  if nargin != 3; err("df: nargin = %d", nargin); end
  if numel(p) != 3; err("df: numel(p)=%d", numel(p)); end
  if numel(q) != 3; err("df: numel(q)=%d", numel(q)); end
  if numel(param) != 3; err("fun: numel(param)=%d", numel(param)); end

  b = param(1); sg = param(2); r = param(3);
  xp = p(1); yp = p(2); zp = p(3);
  x  = q(1);  y = q(2);  z = q(3);

  dx = sg*(yp - xp);
  dy = (r - z)*xp - yp - x*zp;
  dz = y*xp + x*yp - b*zp;

  dp = zeros(size(p));
  dp(1) = dx; dp(2) = dy; dp(3) = dz;
endfunction

function dQ = dF(Q, param)
  if nargin != 2; err(": nargin = %d", nargin); end
  if numel(param) != 3; err("fun: numel(param)=%d", numel(param)); end
  if numel(Q) != 6; err("fun: numel(Q)=%d", numel(Q)); end

  q = Q(1:3); p = Q(4:6);

  dq = f(q, param);
  dp = df(p, q, param);

  dQ = zeros(size(Q));
  dQ(1:3) = dq; dQ(4:6) = dp;
endfunction

function x = norm0(x)
  eps = 1e-12;
  d = sumsq(x);
  if abs(d) > eps; x /= sqrt(d); end
endfunction
function Q = norm(Q)
  Q(4:6) = norm0(Q(4:6));
endfunction

function [to, x] = chain(df, upd, x0,   t0, t, nchunk)
  nt = numel(t);
  no = nchunk*nt - nchunk + 1;
  x =  zeros(no, numel(x0));
  to = zeros(no, 1);
  lo = 1;
  for i = 1:nchunk
    hi = lo + nt - 1;
    x0 = upd(x0);
    x(lo:hi, :) = lsode(df, x0, t0 + t);
    to(lo:hi)   = t0 + t;
    t0 += t(end);
    x0  = x(hi, :); lo = hi;
  end
endfunction

param = [4.0, 16.0, 45.92];

t0 = 0;
t  = (0:0.001:0.001)';
p0 = [1, 0, 0];
q0 = [1, 1, 1];
Q0 = [q0, p0];
Q  = lsode(@(x, t) dF(x, param), Q0, t);

nchunk = 10000;
[t, Q] = chain(@(x, t) dF(x, param),
	       @norm, Q0, t0, t, nchunk);
dlmwrite("q.txt", [t Q], ' ');

q = Q(:, 1:3);
p = Q(:, 4:6);
eps = 2.0;
write_vector("vec", q(1:10:end, :), eps * p(1:10:end, :));
