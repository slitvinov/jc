function farazmand()
  global me;
  me = 'farazmand.m';

  X = 1; Y = 2; Z = 3;
  t0 = 0;
  t  = (0:0.1:1200)';
  param  = [0.01, 2*pi, 0.1, 0.1];
  q0 = [0 0.01 0.01];
				% p0 = [1 1 0]/sqrt(2);
  p0 = f(q0, param); p0 = p0/sqrt(sum(p0.^2, 2));
  Q0 = [q0, p0];

		       % lsode_options('integration method', 'adams');
		       % Q = lsode(@(x, t) dF(x, param), Q0, t);
  tspan = [t(1) t(end)];

  o = odeset(odeset(), 'RelTol', 1e-8);
  [t, Q] = ode45(@(t, x) dF(x, param), tspan, Q0, o);

  q = Q(:, 1:3);
  p = Q(:, 4:6);

  v = f(q, param);
  v = v ./ sqrt(sum(v.^2, 2));

  dlmwrite('q.txt', [t q], ' ');
  dlmwrite('p.txt', [t p], ' ');
  dlmwrite('v.txt', [t v], ' ');
endfunction

function err(fmt, varargin)
  global me
  fmt = ['[%s] ' fmt];
  arg = {me, varargin{:}};
  error(fmt, arg{:});
end

function dq = f(q, param)
  if nargin ~= 2; err('f: nargin = %d', nargin); end
  if size(q, 2) ~= 3; q = q'; end
  if size(q, 2) ~= 3; err('f: size(q)=[%d %d]', size(q)); end
  if numel(param) ~= 4; err('f: numel(param)=%d', numel(param)); end

  al = param(1); om = param(2); la = param(3); be = param(4);
  x = q(:, 1); y = q(:, 2); z = q(:, 3);
  
  dx = al*x + om*y + al*x.^2 + 2*om*x.*y + z.^2;
  dy = -om*x + al.*y - om*x.^2 + 2*al*x.*y;
  dz = -la*z - (la + be)*x.*z;

  dq = zeros(size(q));
  dq(:, 1) = dx; dq(:, 2) = dy; dq(:, 3) = dz;
end

function dp = df(p, q, param)
  if nargin ~= 3; err('df: nargin = %d', nargin); end
  if numel(p) ~= 3; err('df: numel(p)=%d', numel(p)); end
  if numel(q) ~= 3; err('df: numel(q)=%d', numel(q)); end
  if numel(param) ~= 4; err('df: numel(param)=%d', numel(param)); end

  al = param(1); om = param(2); la = param(3); be = param(4);
  x  = q(1);  y = q(2);  z = q(3);
  xp = p(1); yp = p(2); zp = p(3);

  L = [
       2*om*y+2*al*x+al   2*om*x+om               2*z;
       2*al*y-2*om*x-om   2*al*x+al                 0;
       -(la+be)*z                 0   (-(la+be)*x)-la];

  dp = L*p  - (p' * L * p) * p;
end

function dQ = dF(Q, param)
  if nargin ~= 2; err('dF: nargin = %d', nargin); end
  if numel(param) ~= 4; err('dF: numel(param)=%d', numel(param)); end
  if size(Q, 1) ~= 6; err('dF: size(Q)=[%d %d]', size(Q)); end

  q = Q(1:3); p = Q(4:6);

  dq = f(q, param);
  dp = df(p, q, param);

  dQ = zeros(size(Q));
  dQ(1:3) = dq; dQ(4:6) = dp;
end
