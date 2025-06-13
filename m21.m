% ex4-2. speckle pattern (툴박스 없이)
clear, close all
addpath C:\Users\kdo77\Desktop\MUST

%% parameters
probe = 'L11-5v';
must = getparam(probe);
must.attenuation = 0.5; % attenuation [dB/cm/MHz]
must.fs = 100e6; % sampling freq [Hz]
must.c = 1540; % sound speed [m/s]
must.fc = 7.5e6; % 중심 주파수 (필요시 수동 설정)
rng(0)
figure(1), viewxdcr(must)
figure(2), [pulse,t] = getpulse(must);
plot(t/1e-6,pulse)

%% phantom
N = 1000;
W = 5e-3; H = 5e-3; L = 5e-3;
z0 = 10e-3;
xs = rand(N,1)*W - W/2;
ys = rand(N,1)*L - L/2;
zs = rand(N,1)*H - H/2 + z0;
rc = ones(N,1);
figure(3), scatter3(xs*1000,ys*1000,zs*1000,1)
xlabel x(mm), ylabel y(mm), zlabel z(mm)

%% simulation
txdel = txdelay(must,0);
must.fnumber = 2;
px = 0.05e-3;
DR = 30;
x = -2e-3:px:2e-3;
z = 8e-3:px:12e-3;
[zz,xx] = ndgrid(z,x);

tic
B = cell(10,1);
for i = 1:10
    RF = simus(xs, ys, zs, rc, txdel, must);
    SIG = das(RF,xx,zz,txdel,must);
    
    % IQ 신호 생성 (hilbert 수동 구현 + FIR 저역통과 필터)
    IQ = hilbert_manual(SIG);  % analytic signal
    IQ = fir_lowpass(IQ, must.fs, must.fc);  % 대체 필터링

    B{i} = bmode(IQ,DR);
    xs = xs + 2*px;
end
toc

%% B-mode 시각화
figure(4)
montage(B,'size',[1 4])

%% 상관계수 계산 (corrcoef 사용)
corrs = zeros(9,1);
for i = 1:9
    b1 = double(B{i}(:));
    b2 = double(B{i+1}(:));
    r = corrcoef(b1, b2);
    corrs(i) = r(1,2);
end

figure(5)
plot(1:9, corrs, '-o', 'LineWidth', 1.5)
xlabel('프레임 번호'), ylabel('상관계수')
title('인접 B-mode 프레임 간 상관계수')
grid on

%% 함수: Hilbert 수동 구현
function x_analytic = hilbert_manual(x)
    N = size(x,1);
    X = fft(x);
    h = zeros(N,1);
    if mod(N,2)==0
        h([1 N/2+1]) = 1;
        h(2:N/2) = 2;
    else
        h(1) = 1;
        h(2:(N+1)/2) = 2;
    end
    X = bsxfun(@times, X, h);
    x_analytic = ifft(X);
end

%% 함수: FIR 저역통과 필터
function y = fir_lowpass(x, fs, fc)
    n = 50;
    Wn = fc / (fs/2);
    b = fir1(n, Wn, hamming(n+1));
    y = filter(b, 1, x);  % 한쪽 방향 필터 (대신 filtfilt는 사용 불가)
end
