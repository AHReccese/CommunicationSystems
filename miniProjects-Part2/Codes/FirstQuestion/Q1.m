syms t x(t) c(t) r(t) v(t) w
x(t) = exp(-t^2);
fc   = 10;
xfourier(w) = fourier(x(t));
c(t) = cos(2*pi*fc*t);
cfourier(w) = fourier(c(t));
v(t) = x(t)*c(t);
vfourier(w) = fourier(v(t));
r(t) = v(t)*c(t);
syms t
helper(t) = (10/pi)*sinc((10/pi)*t);
syms w
HLP(w) = fourier(helper(t));
syms t
lowpassFilter(t) = HLP(t);


% fourier series
syms w
X(w) = fourier(x);
C(w) = fourier(c);

V(w) = fourier(v);
R(w) = fourier(r);
D(w) = R(w)*HLP(w);
syms w t
a(w,t) = D(w)*exp(1i*w*t);
d(t) = 1/(2*pi)*int(a,w,-inf,inf);

subplot(3,2,1)
fplot(x(t))
title("inputSignal-x(t)")
xlabel('t');
ylabel('x(t)')
grid on

subplot(3,2,2)
fplot(c(t))
title("localOcil-c(t)")
xlabel('t');
ylabel('c(t)')
grid on

subplot(3,2,3)
fplot(v(t))
title("modulated-v(t)")
xlabel('t');
ylabel('v(t)')
grid on

subplot(3,2,4)
fplot(r(t))
title("response-r(t)")
xlabel('t');
ylabel('r(t)')
grid on


subplot(3,2,5)
fplot(helper(t))
title("lowPassFilter-lowPass(t)")
xlabel('t');
ylabel('hlp(t)')
grid on
xlim([-15 15])

subplot(3,2,6)
fplot(d(t))
title("demodulate-d(t)")
xlabel('t');
ylabel('d(t)')
grid on

figure

subplot(3,2,1)
fplot(X(w))
title("inputSignal-X(w)")
xlabel('w');
ylabel('X(w)')
grid on

subplot(3,2,2)
fplot(C(w))
title("localOcil-C(w)")
xlabel('w');
ylabel('C(w)')
grid on
subplot(3,2,3)
fplot(V(w))
title("modulated-v(w)")
xlim([-20 20])
ylim([0 1])
xlabel('w');
ylabel('V(w)')
grid on
subplot(3,2,4)
fplot(R(w))
title("response-R(w)")
xlabel('w');
xlim([-40 40])
ylabel('R(w)')
grid on



subplot(3,2,5)
fplot(HLP(w))
title("lowPassFilter-HLP(w)")
xlim([-20 20])
ylim([0 2])
xlabel('w');
ylabel('HLP(w)')
grid on

subplot(3,2,6)
fplot(D(w))
title("demodulate-D(w)")
xlabel('w');
ylabel('D(w)')
grid on


