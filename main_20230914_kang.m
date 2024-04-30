clear;
% profile on

% per = parpool(16); % Create parallel pool on cluster
% 
ps = 0.1:0.01:1; % turgor pressure P1

delta = 0.5; %細胞壁の厚み P2
nu = 0.5;%  P3
% T = 200000;
T = 50000;

%------画像の出力設定--------------------------
n1 = T/20;
% ic = 0;
c = jet(20+1);
% c = summer(10+1);
%--------------------------------

%-------初期値設定-------------------------
A = 1;
r=5.0;
norm=1.0;
ani=1.2;
%--------------------------------

dt=0.001;%tstep

% 時間発展
N=30;
L1 = length(ps);

% parfor ip = 1:L1
for ip = 11
    P = ps(ip);
    bs = 2.5:0.05:10; % parameter of growth area
    L2 = length(bs);
%     for ib = 1:L2
    for ib = 51
        B = bs(ib);
        ic = 0;
        v = mainfunc(P,A,B,delta,nu,T,n1,ic,c,r,norm,ani,dt,N);
    end
end

% delete(per); % shut down the parallel pool.

%--------------------------------------------------------------------------
function v = mainfunc(P,A,B,delta,nu,T,n1,ic,c,r,norm,ani,dt,N)
%Description:
%P；膨圧
%nu：Poisson's ratioみたい,
%delta：細胞壁の厚み
%n1:　グラフ出力用
%ic: 色設定
%c: colormap
%dt: t-step
%N: 初期値の点の数
%A,B,norm,ani,r初期値のパラメータ

%-----保存先を指定する------------------------------------------------------
ph = '/home/kang/Desktop/data3/'; % savepath

fd = ['p=' num2str(P,'%.2f') '_lg=' num2str(B,'%.2f')];%フォルダを作る
if exist([ph fd],'dir')==0
    mkdir([ph fd])
end

%輪郭座標保存
sxy = [ph fd '/xy/'];
if exist([ph fd '/xy'],'dir')==0
    mkdir([ph fd '/xy'])
end

%曲率情報保存
k = [ph fd '/kappa/'];
if exist([ph fd '/kappa'],'dir')==0
    mkdir([ph fd '/kappa'])
end

%表面張力保存
stress = [ph fd '/sigma/'];
if exist([ph fd '/sigma'],'dir')==0
    mkdir([ph fd '/sigma'])
end

%生長速度保存
sv = [ph fd '/velocity/'];
if exist([ph fd '/velocity'],'dir')==0
    mkdir([ph fd '/velocity'])
end

%細胞壁伸展性保存
sGphi = [ph fd '/Gphi/'];
if exist([ph fd '/Gphi'],'dir')==0
    mkdir([ph fd '/Gphi'])
end

%法線方向と中心線の角度
sphi = [ph fd '/phi/'];
if exist([ph fd '/phi'],'dir')==0
    mkdir([ph fd '/phi'])
end
%--------------------------------------------------------------------------

%初期値を作る---------------------------------------------------------------
theta = zeros(N,1);
x = zeros(N,1);
y = zeros(N,1);

for i = 1:N
    theta(i) = pi*(i-1)/(2.0*(N-1));
    x(i) = norm*r*cos(theta(i));
    y(i) = ani*r*sin(theta(i));
end
%-------------------------------------------


%------------Main-----------------------------
v = zeros(T,1); %生長速度保存
for t=1:T

    epss = zeros(N,1);%接線方向ひずみ率
    epst = zeros(N,1);%法線方向ひずみ率
    f = zeros(N,1);%式(15)(16)の積分の部分を計算するため
    int1 = zeros(N,1);%式(15)(16)の積分の部分を保存して，ｖtとvn計算するため
    vt = zeros(N,1); %接線方向の生長速度
    vn = zeros(N,1);%法線方向の生長速度
    ds = zeros(N,1);%s(i)-s(i-1)
    
    %-曲率の計算-----------------------------------------------------------
    %まず3点を使って接線方向の曲率kappasを計算します。r_s = 1/kappas
    [~,tmpr,tmpk] = curvature([[-x;flip(x)],[y;flip(y)]]);
    kappas = flip(1./tmpr(N+1:end));%子午線方向の曲率（接線方向）

    %phiを計算します。
    tmpalpha = flip(atan(tmpk(N+1:end,2)./tmpk(N+1:end,1)));
    phi = pi/2-tmpalpha;
    
    %r_theta = sin(phi)/x,  xは各店のｘ-axisの座標の値
    %法線方向の曲率kappatは1/r_theta
    kappat = sin(phi)./x;
    kappat(N) = kappas(N);%円環方向の曲率
    %---------------------------------------------------------------------


    %曲率kappasとkappatに基づいて表面張力の計算　式（９）と（１０）
    sigmas = P./(2.*delta.*kappat);
    sigmat = P.*(2-kappas./kappat)./(2.*delta.*kappat);
    %---------------------------------------------------------------------

    %         頂点からの距離の計算
    s1 = 0;
    s2 = zeros(N,1);
    for i = 2:N
        ds(i) = sqrt((x(i)-x(i-1))^2+(y(i)-y(i-1))^2);
        s2(i) = s2(i-1) + ds(i);
        s1 = s1+ s2(i);
    end
    s1 = s2(N,1);
    GPhi = FGPhi(A,B,s2,s1,2);

    integrals1 = 0;
    for i = 2:N-1
        beta=2*nu*nu-2*nu+2; % equation (5)
        K=sqrt(beta*sigmas(i)*sigmas(i)+beta*sigmat(i)*sigmat(i)+(beta-6*nu)*sigmas(i)*sigmat(i)); % equation (5)
        sigmae=sqrt(0.5*(sigmas(i)-sigmat(i))*(sigmas(i)-sigmat(i))+0.5*(sigmat(i))*(sigmat(i))+0.5*(sigmas(i))*(sigmas(i))); %equation (1)
        epss(i)=GPhi(i)*0.5*sigmae*(sigmas(i)-nu*sigmat(i))/K; % equation (5)
        epst(i)=GPhi(i)*0.5*sigmae*(sigmat(i)-nu*sigmas(i))/K; % equation (5)

        % equation (16)
        if i==2
            f(i)=((epss(i-1)-kappas(i-1)*epst(i-1)/kappat(i-1))/sin(pi/2.0)+(epss(i)-kappas(i)*epst(i)/kappat(i))/sin(phi(i)))*(ds(i)/2.0);
            integrals1=integrals1+f(i);
        elseif i>2
            f(i)=((epss(i-1)-kappas(i-1)*epst(i-1)/kappat(i-1))/sin(phi(i-1))+(epss(i)-kappas(i)*epst(i)/kappat(i))/sin(phi(i)))*(ds(i)/2.0);
            integrals1=integrals1+f(i);
        end
        int1(i)=integrals1;
        vt(i)=sin(phi(i))*int1(i);
        vn(i)=epst(i)/kappat(i)-cos(phi(i))*int1(i);
    end

    i = i+1;
    if i==N
        beta = 2*nu*nu-2*nu+2; %equation (5)
        K = sqrt(beta*sigmas(i)*sigmas(i)+beta*sigmat(i)*sigmat(i)+(beta-6*nu)*sigmas(i)*sigmat(i)); %equation (5)
        sigmae = sqrt(0.5*(sigmas(i)-sigmat(i))*(sigmas(i)-sigmat(i))+0.5*(sigmat(i))*(sigmat(i))+0.5*(sigmas(i))*(sigmas(i))); %equation (1)
        epss(i) = GPhi(i)*0.5*sigmae*(sigmas(i)-nu*sigmat(i))/K; % equation (5)
        epst(i) = GPhi(i)*0.5*sigmae*(sigmat(i)-nu*sigmas(i))/K; % equation (5)

        % equation (16)
        integrals1 = integrals1+((epss(i-1)-kappas(i-1)*epst(i-1)/kappat(i-1))/sin(phi(i-1))+0.0)*(ds(i)/2.0); 
        int1(i) = integrals1;
        vt(i) = 0.0;
        vn(i) = epst(i)/kappat(i)-int1(i);
    end

    %   形状変化
    x(1)=x(1);
    y(1)=y(1);    
    for i = 2:N
        x(i)=x(i)+vt(i)*cos(phi(i))*dt+vn(i)*sin(phi(i))*dt;
        y(i)=y(i)+vt(i)*sin(phi(i))*dt+vn(i)*cos(phi(i))*dt;
    end

    if mod(t-1,n1)==0 % save and plot
        ic = ic + 1;
        plot(x,y,'.','color',c(ic,:),LineWidth=1.5,MarkerSize=10)
        hold on
        plot(-x,y,'.','color',c(ic,:),LineWidth=1.5,MarkerSize=10)
        axis equal;
        pause(0.01)
%         parsave([sxy 't=' num2str(ic,'%.3d') '.txt'],x(cc:end),y(cc:end),...
%             [k 't=' num2str(ic,'%.3d') '.txt'],kappas(cc:end),kappat(cc:end),...
%             [stress 't=' num2str(ic,'%.3d') '.txt'],sigmas(cc:end),sigmat(cc:end),...
%             [sv 't=' num2str(ic,'%.3d') '.txt'],vt(cc:end),vn(cc:end),...
%             [sphi 't=' num2str(ic,'%.3d') '.txt'],phi(cc:end),...
%             [sGphi 't=' num2str(ic,'%.3d') '.txt'],GPhi(cc:end));
    end

    %ds=0.35を指定します。dsの爆発と特異点による爆発を回避するため
    [x,y,N] = interpolation(x,y,ds,N,GPhi,0.35,t); %interpolation

    v(t) = vn(end);
end
axis off

figure(2) % plot growth rate
plot(v(200:end),'k',LineWidth=2)
xticks([])
yticks([])
xlabel('Time')
ylabel('Growth rate')
end

function [r,theta] = polarCoordinate(x,y,N)
r = zeros(N,1);
theta = zeros(N,1);
for i = 1:N
    r(i) = sqrt(x(i)^2+y(i)^2);
    theta(i) = atan(y(i)/x(i));
end
end

function [x,y] = xyCoordinate(r,theta,N)
x = zeros(N,1);
y = zeros(N,1);
for i = 1:N
    x(i) = r(i)*cos(theta(i));
    y(i) = r(i)*sin(theta(i));
end
end

function y = FGPhi(A,lg,s,S,n) %細胞壁伸展性の分布
% n = 1;
y = zeros(length(s),1);
for i = 1:length(s)
    tmps = S-s(i);
    if tmps<=lg
        y(i,1) = A*cos((pi*tmps)/(2*lg))^n;
%         y(i) = 1;
%         y(i) = (1-tmps/lg)^n;
    else
        y(i,1) = 0;
    end
end
end

function [x_new,y_new,N_new] = interpolation(x,y,ds,N,Gphi,threshold,t)
s = zeros(N,1);
for i=2:N
    s(i,1) = s(i-1)+ds(i);
end
[r,theta] = polarCoordinate(x,y,N);

s_new = s(N):-threshold:s(1);
N_new = length(s_new);
s_new(N_new) = s(1);
s_new = flip(s_new);

theta_new = interp1(s,abs(theta),s_new,'spline');
r_new = interp1(s,r,s_new,'spline');
[x_new,y_new] = xyCoordinate(r_new,theta_new,N_new);
end






















