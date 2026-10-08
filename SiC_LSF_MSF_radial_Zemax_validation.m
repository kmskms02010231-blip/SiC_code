%% SiC_LSF_MSF_radial_Zemax_validation.m
% =========================================================================
%  LSF + MSF 랜덤 표면 + radial ripple로 MTF를 세 가지 방법으로 계산해 서로 비교하는 스크립트
%
%   [A] 연구 모델 : PSD -> ACV(한켈) -> 구조함수 -> STF -> MTF   (검증 대상)
%   [B] MATLAB    : PSD -> 실제 2D 표면 -> 반사 위상 -> pupil 자기상관 OTF
%                   (ACV/STF 공식을 쓰지 않음, 여러 장 평균)
%   [C] Zemax     : [B]와 같은 표면 -> 위상맵(.DAT) -> OpticStudio Grid Phase
%                   -> Zemax FFT MTF   (ZOS-API로 MATLAB에서 자동 실행)
%
%  PSD 조건은 연구 코드 1~7절과 동일
%   - LSF : Lorentzian 1/(1+(f/rho_hp)^p), rho_hp = 0.005, p = 2.60, 15 nm
%   - MSF : LSF 끝값에서 시작하는 멱함수 S_L_edge*(f/f_LM)^(-p_M), 13 nm
%           (p_M은 fzero로 13 nm가 되도록 결정)
%   - Radial ripple : 0.21 cycles/mm, 5 nm RMS (연구 코드 4-1절 ripple_radial 첫 줄)
%                     sigma_M = 13 nm 안에 ripple 포함 -> MSF 배경 = sqrt(13^2-5^2) = 12 nm
%                     배경 PSD는 연구 코드 7절의 alpha_M(cosine 전환) 보정 그대로
%                     ACV 모델: C_r(r) = sigma_r^2 * J0(2*pi*fr*r)  (연구 코드 10-3절)
%   - Raster ripple, HSF 제외
%
%  판정 기준
%   - A vs B : 모든 주파수에서 |B-A| < 2*표준오차 이면 연구 모델 통과
%   - B vs C : 같은 표면끼리 RMSE ~1e-3 이하이면 MATLAB 계산이 Zemax와 일치
%
%  필요한 것 : MATLAB R2019b 이상, Ansys Zemax OpticStudio (ZOS-API 라이선스)
%  실행 시간 : 표면 150장 약 10분 + Zemax 10장 약 8분 (RAM 16 GB 이상 권장)
% =========================================================================
clc; clearvars; close all;

%% 1. 광학계 파라미터
D          = 2000;                                % 유효 구경 지름 [mm]
R          = D/2;                                 % 구경 반지름 [mm]
lambda_nm  = 650;                                 % 파장 [nm]
lambda     = lambda_nm*1e-6;                      % 파장 [mm]
Fnum       = 13.91049;                            % F-number
f_eff      = Fnum*D;                              % 유효 초점거리 [mm] (= 27,820.98 mm)
theta_deg  = 0;                                   % 입사각 [deg]
pixel_um   = 7;                                   % 검출기 픽셀 피치 [um]
nu_nyq     = 1/(2*pixel_um*1e-3);                 % 영상 Nyquist 주파수 [cycles/mm] (= 71.43)
nu_c       = 1/(lambda*Fnum);                     % 회절 차단 주파수 [cycles/mm] (= 110.6)
q          = 4*pi*cosd(theta_deg)/lambda_nm;      % 표면 높이 -> 반사 위상 변환 상수 [rad/nm]
                                                  % (반사면은 광경로가 2h라서 4*pi*h/lambda)

%% 2. 대역 경계 및 PSD 형태 파라미터 (연구 코드 1절과 동일)
rho_hp = 0.005;                                   % Lorentzian knee 주파수 [cycles/mm]
p      = 2.60;                                    % Lorentzian 고주파 감쇠 지수 [-]

f_min  = 1/D;                                     % LSF 시작 (구경당 1주기) [cycles/mm]
f_LM   = 7/D;                                     % LSF-MSF 경계 (구경당 7주기) [cycles/mm]
f_MH   = sqrt(10/(lambda*1378.014));              % MSF-HSF 경계 [cycles/mm] (= 3.341)

sigma_L_nm = 15;                                  % LSF 목표 RMS [nm]
sigma_M_nm = 13;                                  % MSF 목표 RMS [nm] (ripple 포함 총량)

%% 2-1. Radial ripple 파라미터 (연구 코드 4-1절)
fr_radial  = 0.21;                                % radial ripple 주파수 [cycles/mm] (주기 4.76 mm)
sig_radial = 5;                                   % radial ripple RMS [nm]
A_radial   = sqrt(2)*sig_radial;                  % 사인파 진폭(peak) [nm] (RMS = A/sqrt(2))
sigma_M_bg_nm = sqrt(sigma_M_nm^2 - sig_radial^2);% MSF 배경(연속 PSD)이 채울 RMS [nm] (= 12 nm)

%% 3. LSF PSD 정규화 계수 h_L (연구 코드 2절과 동일)
% 40 m x 40 m 가상 도메인, 2048 x 2048 격자의 LSF 대역 합으로 정규화
widthXY = 40000;                                  % [공간] 가상 도메인 폭 [mm]
ss      = 2048;                                   % 격자 한 변 점 수 [-]
A_area  = widthXY^2;                              % [공간] 전체 면적 [mm^2]

u   = ((0:ss-1)-ss/2)/widthXY;                    % [주파수] 1D 축 [cycles/mm]
[U, V] = meshgrid(u, u);
rho = sqrt(U.^2 + V.^2);                          % [주파수] 반경 주파수 [cycles/mm]

base_lorentz = 1./(1+(rho/rho_hp).^p);
mask_L = (rho >= f_min) & (rho <= f_LM);          % LSF 대역 마스크
h_L    = sum(base_lorentz(mask_L));               % LSF 정규화 계수

n_bins = round((f_LM - f_min)/(1/widthXY));       % 연구 코드와 같은 bin 수
clear U V rho base_lorentz mask_L

%% 4. LSF / MSF PSD 함수 정의 (연구 코드 5~7절과 동일)
% LSF PSD [nm^2 mm^2]
S_L_at   = @(f) (sigma_L_nm^2*A_area/h_L) * (1./(1+(f/rho_hp).^p));
S_L_edge = S_L_at(f_LM);                          % LSF PSD의 f_LM 값 -> MSF 시작값

% MSF PSD: S_L_edge*(f/f_LM)^(-p_M), 분산이 sigma_M^2가 되도록 p_M 결정
f_edges_M  = logspace(log10(f_LM), log10(f_MH), n_bins+1);
f_axis_M   = [f_LM, sqrt(f_edges_M(1:end-1).*f_edges_M(2:end)), f_MH];
var_M_of_p = @(pp) 2*pi*trapz(f_axis_M, S_L_edge*(f_axis_M/f_LM).^(-pp).*f_axis_M);
p_M        = fzero(@(pp) var_M_of_p(pp) - sigma_M_nm^2, p);
S_M_nat    = @(f) S_L_edge*(f/f_LM).^(-p_M);      % 13 nm 기준 MSF PSD (모양 결정용)

% 연구 코드 7절: f_LM에서 LSF와 연속 -> f_trans까지 cosine 전환 -> 이후 alpha_M배로 평탄하게 낮춤
% 배경 분산이 sigma_M_bg^2가 되도록 alpha_M 결정
trans_frac = 0.15;
f_trans    = f_LM*(f_MH/f_LM)^trans_frac;
t_of_f     = @(f) min(max(log(f/f_LM)/log(f_trans/f_LM), 0), 1);
g_of_f     = @(f, a) a + (1-a)*0.5*(1+cos(pi*t_of_f(f)));
var_M_of_a = @(a) 2*pi*trapz(f_axis_M, g_of_f(f_axis_M,a).*S_M_nat(f_axis_M).*f_axis_M);
alpha_M    = fzero(@(a) var_M_of_a(a) - sigma_M_bg_nm^2, sigma_M_bg_nm^2/sigma_M_nm^2);
S_M_fun    = @(f) g_of_f(f,alpha_M) .* S_M_nat(f);   % MSF 배경 PSD [nm^2 mm^2]

fprintf('[PSD]  h_L = %.4g,  p_M = %.4f,  alpha_M = %.4f\n', h_L, p_M, alpha_M);
fprintf('       LSF 분산 = %.2f nm^2 (목표 %.0f)\n', ...
    2*pi*integral(@(f) S_L_at(f).*f, f_min, f_LM), sigma_L_nm^2);
fprintf('       MSF 배경 분산 = %.2f nm^2 (목표 %.0f) + radial ripple %.0f nm^2\n', ...
    2*pi*integral(@(f) S_M_fun(f).*f, f_LM, f_MH), sigma_M_bg_nm^2, sig_radial^2);

%% 5. 시뮬레이션 설정
n_grid    = 2049;                                 % 2 m 표면 격자 한 변 점 수 [-]
N_ens     = 150;                                  % 평균에 쓸 랜덤 표면 장 수 [-]
                                                  % (1장 표준편차 ~0.0186 -> 표준오차 0.0015 목표: N ~ 150)
N_zemax   = 10;                                   % 그중 Zemax로 1:1 비교할 장 수 [-]
zemax_fft = 2048;                                 % Zemax FFT MTF pupil 샘플링 (256~4096)
rng_seed  = 20261007;                             % 난수 시드 (같은 값이면 같은 표면)
run_zemax = true;                                 % false면 Zemax 단계 생략

out_dir = fullfile(pwd, 'out_LSF_MSF_radial_Zemax');     % 결과 저장 폴더
dat_dir = fullfile(out_dir, 'grid_phase_maps');   % Zemax 입력 위상맵 폴더
if ~exist(dat_dir,'dir'), mkdir(dat_dir); end

%% 6. 표면 격자, pupil, lag 축
dx       = D/(n_grid-1);                          % 표면 격자 간격 [mm] (= 0.977 mm)
f_map_ny = 1/(2*dx);                              % 표면 격자 Nyquist [cycles/mm] (= 0.512)
f_M_max  = min(f_MH, 0.96*f_map_ny);              % 격자로 표현할 MSF 상한 [cycles/mm]
                                                  % (p_M이 커서 MSF 분산의 99.999%가 이 아래)
fprintf('[격자]  dx = %.4f mm,  MSF 사용 상한 = %.4f /mm,  MSF 분산 포함률 = %.4f %%\n', ...
    dx, f_M_max, 100*integral(@(f) S_M_fun(f).*f, f_LM, f_M_max)/integral(@(f) S_M_fun(f).*f, f_LM, f_MH));

if fr_radial > 0.8*f_map_ny
    error('radial ripple %.3f /mm가 격자로 표현 불가 (격자 Nyquist %.3f /mm) -> n_grid를 키우세요', fr_radial, f_map_ny);
end
fprintf('[격자]  radial ripple 한 주기당 %.1f 점\n', 1/(fr_radial*dx));

x_mm  = linspace(-R, R, n_grid);                  % [공간] x 축 [mm]
y_mm  = fliplr(x_mm);                             % [공간] y 축 [mm] (DAT 첫 행 = +y)
pupil = hypot(x_mm, y_mm.') <= R;                 % 원형 pupil 마스크

% 영상 주파수 nu <-> 표면 lag r 관계: r = lambda*f_eff*nu
r_nyq = lambda*f_eff*nu_nyq;                      % Nyquist에 해당하는 lag [mm] (= 1,292 mm)
n_lag = floor(r_nyq/dx) + 3;                      % Nyquist보다 조금 넘게 계산 (내삽용)
r_lag = (0:n_lag-1)*dx;                           % lag 축 [mm]
nu    = r_lag/(lambda*f_eff);                     % 영상 주파수 축 [cycles/mm]
in    = r_lag <= r_nyq;                           % 0 ~ Nyquist 구간
L_fft = 2^nextpow2(n_grid + n_lag);               % 자기상관 FFT 길이 (원형 wrap 방지)

% 대역별 PSD (격자 범위로 자름). f=0에서 Inf*0 = NaN이 되지 않도록 max 처리
W_L = @(f) S_L_at(f) .* (f>=f_min & f<=f_LM);
W_M = @(f) S_M_fun(max(f,f_LM)) .* (f>=f_LM & f<=f_M_max);

%% 7. [A] 연구 모델: PSD -> ACV -> 구조함수 -> STF -> MTF
% 한켈 변환 C(r) = 2*pi*∫ S(f) J0(2*pi*f*r) f df  (음수 ACV 그대로 = raw)
C_L = hankel_acv(W_L, f_min, f_LM,  r_lag);       % LSF ACV [nm^2]
C_M = hankel_acv(W_M, f_LM,  f_M_max, r_lag);     % MSF ACV [nm^2]

C_R = sig_radial^2 * besselj(0, 2*pi*fr_radial*r_lag);   % radial ripple ACV (연구 코드 10-3절)

% raw ACV
Dh_raw  = 2*((C_L(1)+C_M(1)+C_R(1)) - (C_L+C_M+C_R));     % 구조함수 [nm^2]
STF_raw = exp(-0.5*q^2*Dh_raw);                   % STF

% trust ACV (연구 코드 11절 방식: 처음 0 이하가 되는 지점부터 0) - 비교용
C_Lt = C_L; k0 = find(C_Lt<=0,1,'first'); if ~isempty(k0), C_Lt(k0:end) = 0; end
C_Mt = C_M; k0 = find(C_Mt<=0,1,'first'); if ~isempty(k0), C_Mt(k0:end) = 0; end
Dh_tr   = 2*((C_Lt(1)+C_Mt(1)+C_R(1)) - (C_Lt+C_Mt+C_R));   % ripple은 연구 코드처럼 trust 미적용
STF_tr  = exp(-0.5*q^2*Dh_tr);

% 회절 MTF: 해석식 / 이산 pupil (B와 같은 격자)
nn       = min(nu/nu_c, 1);
MTF_diff = (2/pi)*(acos(nn) - nn.*sqrt(1-nn.^2));             % 해석식
[c0, ~]  = otf_cuts(zeros(n_grid), pupil, q, n_lag, L_fft);
MTF_core = real(c0);                                          % 이산 pupil 회절 MTF

MTF_A_raw = MTF_core .* STF_raw;                  % 연구 모델 (raw ACV)
MTF_A_tr  = MTF_core .* STF_tr;                   % 연구 모델 (trust ACV)
fprintf('[A]  ACV 최솟값: LSF %.2f nm^2, MSF %.2f nm^2 (음수 보존)\n', min(C_L), min(C_M));

%% 8. [B] MATLAB: 랜덤 표면 N_ens장 생성 -> 각 장의 OTF -> 평균
% 표면은 3 부분으로 만들어 더함 (경계 반복 효과를 없애려고 2 m보다 넓게 생성 후 잘라냄)
%  - LSF          : 간격 7.8 mm, 16 m 영역에서 생성 -> 2049 격자로 spline 보간
%  - MSF <= 0.05  : 간격 1.95 mm, 8 m 영역에서 생성 -> 보간
%  - MSF  > 0.05  : 간격 dx, 약 4 m 영역에서 바로 생성
f_split = 0.05;                                   % MSF 저/고주파 나누는 주파수 [cycles/mm]
dx_L = D/256;  M_L = 2048;                        % LSF 생성 격자
dx_m = D/1024; M_m = 4096;                        % 저주파 MSF 생성 격자
M_h  = 2^nextpow2(n_grid);                        % 고주파 MSF 생성 격자
xL = ((0:M_L-1)-M_L/2)*dx_L;  iL = find(abs(xL) <= R+4*dx_L);   % 2 m 부근만 잘라낼 인덱스
xm = ((0:M_m-1)-M_m/2)*dx_m;  im = find(abs(xm) <= R+4*dx_m);
W_M_lo = @(f) W_M(f).*(f<=f_split);
W_M_hi = @(f) W_M(f).*(f> f_split);

rng(rng_seed, 'twister');
OTF_x = complex(zeros(N_ens, n_lag));            % 각 장의 x방향 OTF
OTF_y = complex(zeros(N_ens, n_lag));            % 각 장의 y방향 OTF
RMS   = zeros(N_ens, 4);                          % [LSF, MSF배경, ripple, 합] pupil RMS [nm]
r_grid = hypot(x_mm, y_mm.');                     % 각 격자점의 중심 거리 [mm] (radial ripple용)

for k = 1:N_ens
    % (1) LSF 표면
    hL    = synth_surface(W_L, dx_L, M_L);
    G     = griddedInterpolant({xL(iL), xL(iL)}, hL(iL,iL), 'spline');
    h_lsf = G({x_mm, x_mm});

    % (2) MSF 표면 = 저주파 + 고주파
    hm    = synth_surface(W_M_lo, dx_m, M_m);
    G     = griddedInterpolant({xm(im), xm(im)}, hm(im,im), 'spline');
    h_msf = G({x_mm, x_mm});
    hh    = synth_surface(W_M_hi, dx, M_h);
    h_msf = h_msf + hh(1:n_grid, 1:n_grid);

    % (3) radial ripple: 거울 중심 기준 동심원 사인파, 위상은 장마다 무작위
    %     (연구 코드 ACV 모델 sigma^2*J0는 위상 앙상블 평균에 해당)
    phi   = 2*pi*rand;
    h_rip = A_radial*sin(2*pi*fr_radial*r_grid + phi);

    % (4) 합친 표면 (pupil 평균 제거; pupil 밖 값도 유지 -> Zemax 경계 불연속 방지)
    h = h_lsf + h_msf + h_rip;
    h = h - mean(h(pupil));

    % (4) 이 표면의 OTF (반사 위상 exp(i*q*h)의 pupil 자기상관)
    [OTF_x(k,:), OTF_y(k,:)] = otf_cuts(h, pupil, q, n_lag, L_fft);
    RMS(k,:) = [std(h_lsf(pupil),1), std(h_msf(pupil),1), std(h_rip(pupil),1), std(h(pupil),1)];

    % (5) 앞 N_zemax장은 Zemax 입력 위상맵으로 저장: phase = 4*pi*h/lambda [rad]
    if run_zemax && k <= N_zemax
        write_grid_phase_dat(fullfile(dat_dir, sprintf('combined_%03d.dat',k)), q*h, dx);
    end
    if k == 1
        save(fullfile(out_dir,'surface_001.mat'), 'x_mm','y_mm','pupil','h_lsf','h_msf','h_rip','h','-v7.3');
    end
    fprintf('  표면 %3d/%d   RMS LSF/MSF배경/ripple/합 = %.2f / %.2f / %.2f / %.2f nm\n', k, N_ens, RMS(k,:));
end
clear hL hm hh h_lsf h_msf h_rip h G

% 앙상블 평균: 복소 OTF를 먼저 평균하고 절댓값 (통계 모델 A와 같은 정의)
MTF_B_x = abs(mean(OTF_x,1));
MTF_B_y = abs(mean(OTF_y,1));
SE_x    = std(OTF_x,0,1)/sqrt(N_ens);             % 평균의 표준오차
SE_y    = std(OTF_y,0,1)/sqrt(N_ens);

%% 9. A vs B 비교
idx   = find(in); idx = idx(2:end);               % 0 제외, Nyquist까지
at_ny = @(c) interp1(nu, c, nu_nyq);              % Nyquist 값 (내삽)
z_raw = max(abs([(MTF_B_x(idx)-MTF_A_raw(idx))./SE_x(idx), (MTF_B_y(idx)-MTF_A_raw(idx))./SE_y(idx)]));
z_tr  = max(abs([(MTF_B_x(idx)-MTF_A_tr(idx))./SE_x(idx),  (MTF_B_y(idx)-MTF_A_tr(idx))./SE_y(idx)]));

fprintf('\n[A vs B]  표면 %d장, 평균 RMS LSF/MSF배경/ripple/합 = %.2f / %.2f / %.2f / %.2f nm\n', N_ens, mean(RMS));
fprintf('  raw ACV   : RMSE = %.5f,  최대 |z| = %.2f  (2 이하면 통과)\n', ...
    sqrt(mean((MTF_B_x(idx)-MTF_A_raw(idx)).^2)), z_raw);
fprintf('  trust ACV : RMSE = %.5f,  최대 |z| = %.2f\n', ...
    sqrt(mean((MTF_B_x(idx)-MTF_A_tr(idx)).^2)), z_tr);
fprintf('  Nyquist MTF: 회절 %.4f | A(raw) %.4f | A(trust) %.4f | B %.4f / %.4f (SE %.4f)\n', ...
    at_ny(MTF_diff), at_ny(MTF_A_raw), at_ny(MTF_A_tr), at_ny(MTF_B_x), at_ny(MTF_B_y), at_ny(SE_x));

%% 10. [C] Zemax: 위상맵 -> Grid Phase -> FFT MTF  (ZOS-API)
% 광학계: 무한물체 -> 1면 Grid Phase + Stop (EPD 2000 mm) -> 2면 이상렌즈 f = f_eff -> 상면
% 주의: Zemax Tangential = y방향 주파수, Sagittal = x방향 주파수
%       (x방향 줄무늬 표면에서 Sagittal만 떨어지는 것으로 확인)
if run_zemax
    write_grid_phase_dat(fullfile(dat_dir,'flat.dat'), zeros(n_grid), dx);        % 평면(기준) 위상맵

    % (1) OpticStudio 실행
    zroot = winqueryreg('HKEY_CURRENT_USER','Software\Zemax','ZemaxRoot');
    NET.addAssembly(fullfile(zroot,'ZOS-API','Libraries','ZOSAPI_NetHelper.dll'));
    if ZOSAPI_NetHelper.ZOSAPI_Initializer.Initialize() ~= 1, error('OpticStudio 초기화 실패'); end
    zdir = char(ZOSAPI_NetHelper.ZOSAPI_Initializer.GetZemaxDirectory());
    NET.addAssembly(fullfile(zdir,'ZOSAPI_Interfaces.dll'));
    NET.addAssembly(fullfile(zdir,'ZOSAPI.dll'));
    conn = ZOSAPI.ZOSAPI_Connection();
    app  = conn.CreateNewApplication();
    if isempty(app) || ~app.IsValidLicenseForAPI, error('ZOS-API 라이선스 없음'); end
    sys      = app.PrimarySystem;
    data_dir = char(app.ZemaxDataDir);            % Grid Phase는 이 폴더 기준으로 불러옴

    case_name = [arrayfun(@(k) sprintf('combined_%03d',k), 1:N_zemax, 'UniformOutput', false), {'flat'}];
    Z = cell(numel(case_name), 1);                % 각 case의 [freq, T, S]

    for c = 1:numel(case_name)
        tic;
        % (2) 새 시스템: 구경, 파장, 축상 시야
        sys.New(false);
        sys.SystemData.Aperture.ApertureValue = D;
        sys.SystemData.Wavelengths.GetWavelength(1).Wavelength = lambda_nm/1000;
        sys.LDE.InsertNewSurfaceAt(2);

        % (3) 1면 = Grid Phase + Stop, 위상맵 DAT 불러오기
        s1 = sys.LDE.GetSurfaceAt(1);
        imp = fullfile(data_dir, ['SIC_' case_name{c} '.DAT']);
        copyfile(fullfile(dat_dir,[case_name{c} '.dat']), imp, 'f');
        s1.ChangeType(s1.GetSurfaceTypeSettings(ZOSAPI.Editors.LDE.SurfaceType.GridPhase));
        s1.ImportData.ImportDataFile(System.String(imp));
        s1.IsStop = true;  s1.Thickness = 0;  s1.Comment = case_name{c};

        % (4) 2면 = 이상렌즈(Paraxial), 초점거리 f_eff, 뒤 두께 f_eff
        s2 = sys.LDE.GetSurfaceAt(2);
        s2.ChangeType(s2.GetSurfaceTypeSettings(ZOSAPI.Editors.LDE.SurfaceType.Paraxial));
        s2.GetSurfaceCell(ZOSAPI.Editors.LDE.SurfaceColumn.Par1).DoubleValue = f_eff;
        s2.Thickness = f_eff;

        if c == 1 || c == numel(case_name)       % 1번 표면과 flat만 .zos 저장 (용량 큼)
            sys.SaveAs(System.String(fullfile(out_dir,[case_name{c} '.zos'])));
        end

        % (5) FFT MTF (0 ~ Nyquist)
        mtf = sys.Analyses.New_FftMtf();
        st  = mtf.GetSettings();
        st.MaximumFrequency = nu_nyq;
        st.SampleSize = fft_sample_enum(zemax_fft);
        mtf.ApplyAndWaitForCompletion();
        rr = mtf.GetResults();
        ds = rr.DataSeries(1);
        fr = ds.XData.Data.double;  yy = ds.YData.Data.double;
        Z{c} = [fr(:), yy];                       % [freq, Tangential, Sagittal]
        mtf.Close();
        delete(imp);
        fprintf('  Zemax %-13s  %.0f s   Nyquist T/S = %.5f / %.5f\n', case_name{c}, toc, Z{c}(end,2), Z{c}(end,3));
    end
    app.CloseApplication();

    %% 11. B vs C 비교 (같은 표면끼리 1:1)
    fz = Z{end}(:,1);                             % Zemax 주파수 축
    Zflat_S = Z{end}(:,3);
    err = zeros(N_zemax, 4);                      % [RMSE S-X, RMSE T-Y, Nyq S-X, Nyq T-Y]
    ZS = zeros(N_zemax, numel(fz)); ZT = ZS; MX = ZS; MY = ZS;
    for k = 1:N_zemax
        ZS(k,:) = interp1(Z{k}(:,1), Z{k}(:,3), fz);   % Zemax Sagittal
        ZT(k,:) = interp1(Z{k}(:,1), Z{k}(:,2), fz);   % Zemax Tangential
        MX(k,:) = interp1(nu, abs(OTF_x(k,:)), fz);    % 같은 표면의 MATLAB X
        MY(k,:) = interp1(nu, abs(OTF_y(k,:)), fz);    % 같은 표면의 MATLAB Y
        err(k,:) = [sqrt(mean((ZS(k,:)-MX(k,:)).^2)), sqrt(mean((ZT(k,:)-MY(k,:)).^2)), ...
                    ZS(k,end)-MX(k,end), ZT(k,end)-MY(k,end)];
    end
    fprintf('\n[B vs C]  같은 표면 %d장\n', N_zemax);
    fprintf('  RMSE (S-X) = %.5f ~ %.5f,  RMSE (T-Y) = %.5f ~ %.5f\n', min(err(:,1)), max(err(:,1)), min(err(:,2)), max(err(:,2)));
    fprintf('  Nyquist 차 (Zemax-MATLAB) = %.5f ~ %.5f\n', min(err(:,3:4),[],'all'), max(err(:,3:4),[],'all'));
    fprintf('  Zemax flat Nyquist = %.5f (해석값 %.5f, 샘플링 한계로 약간 낮음)\n', Zflat_S(end), at_ny(MTF_diff));
end

%% 12. 그래프
% (1) Zemax에 넣은 1번 표면
S1 = load(fullfile(out_dir,'surface_001.mat'));
figure(1); set(gcf,'Position',[50 50 1800 420]);
H = {S1.h_lsf, S1.h_msf, S1.h_rip, S1.h};  T = {'LSF','MSF 배경','Radial ripple','합 (Zemax 입력)'};
for k = 1:4
    hk = H{k} - mean(H{k}(S1.pupil));  hk(~S1.pupil) = NaN;
    subplot(1,4,k); imagesc(S1.x_mm, S1.y_mm, hk); axis image xy; colorbar;
    title(sprintf('%s  RMS %.2f nm', T{k}, std(hk(S1.pupil),1))); xlabel('x [mm]'); ylabel('y [mm]');
end

% (2) A vs B vs C
figure(2); set(gcf,'Position',[80 80 1400 480]);
subplot(1,2,1);
plot(nu(in), MTF_diff(in), '-', 'Color', [.6 .6 .6]); hold on;
plot(nu(in), MTF_A_raw(in), 'k-', 'LineWidth', 2);
plot(nu(in), MTF_A_tr(in),  'm:', 'LineWidth', 2);
plot(nu(in), MTF_B_x(in),   'b--', 'LineWidth', 1.3);
leg = {'회절한계','A: raw ACV','A: trust ACV','B: MATLAB 평균'};
if run_zemax, plot(fz, mean(ZS,1), 'r-.', 'LineWidth', 1.3); leg{end+1} = 'C: Zemax 평균'; end
xline(nu_nyq, ':r', 'HandleVisibility','off'); grid on;
xlabel('Image frequency \nu [cycles/mm]'); ylabel('MTF'); legend(leg, 'Location','northeast');
title(sprintf('LSF %g + MSF %g (배경 %g + radial %g nm @ %.2f/mm)', sigma_L_nm, sigma_M_nm, sigma_M_bg_nm, sig_radial, fr_radial));
subplot(1,2,2);
fill([nu(in) fliplr(nu(in))], [2*SE_x(in) fliplr(-2*SE_x(in))], [.88 .88 .88], 'EdgeColor','none'); hold on;
plot(nu(in), MTF_B_x(in)-MTF_A_raw(in), 'b', nu(in), MTF_B_x(in)-MTF_A_tr(in), 'm', 'LineWidth', 1.3);
grid on; xlabel('\nu [cycles/mm]'); ylabel('B - A');
legend({'±2 표준오차','B - A(raw)','B - A(trust)'}, 'Location','southeast');
title('회색 띠 안에 있으면 통과');

% (3) 같은 표면 1장: MATLAB vs Zemax (오른쪽은 회절한계로 나눈 STF)
if run_zemax
    figure(3); set(gcf,'Position',[110 110 1400 480]);
    core_z = interp1(nu, MTF_core, fz);           % 열벡터
    subplot(1,2,1);
    plot(fz, core_z, '-', 'Color', [.6 .6 .6]); hold on;
    plot(fz, MX(1,:), 'k-', fz, ZS(1,:), 'b--', fz, MY(1,:), '-', 'Color', [.4 .4 .4]);
    plot(fz, ZT(1,:), 'r--'); grid on;
    legend({'회절한계','MATLAB X','Zemax Sagittal','MATLAB Y','Zemax Tangential'});
    xlabel('\nu [cycles/mm]'); ylabel('MTF'); title('표면 #1: MATLAB vs Zemax');
    subplot(1,2,2);
    plot(fz, MX(1,:)./core_z.', 'k-', fz, ZS(1,:)./Zflat_S.', 'b--', fz, MY(1,:)./core_z.', '-', 'Color', [.4 .4 .4]);
    hold on; plot(fz, ZT(1,:)./Zflat_S.', 'r--'); grid on; ylim([0.8 1.0]);
    legend({'MATLAB X / 회절','Zemax S / Zemax flat','MATLAB Y / 회절','Zemax T / Zemax flat'});
    xlabel('\nu [cycles/mm]'); ylabel('STF = MTF / 회절 MTF (항상 1 이하)'); title('표면에 의한 저하만 비교');
end

%% 13. 결과 저장
save(fullfile(out_dir,'results.mat'), 'nu','in','MTF_diff','MTF_core','MTF_A_raw','MTF_A_tr', ...
    'MTF_B_x','MTF_B_y','SE_x','SE_y','C_L','C_M','C_R','STF_raw','STF_tr','RMS','p_M','h_L','alpha_M');
writetable(table(nu(in).', MTF_diff(in).', MTF_A_raw(in).', MTF_A_tr(in).', MTF_B_x(in).', MTF_B_y(in).', SE_x(in).', ...
    'VariableNames', {'nu_cyc_per_mm','MTF_diff','MTF_A_raw','MTF_A_trust','MTF_B_X','MTF_B_Y','StdErr'}), ...
    fullfile(out_dir,'A_vs_B.csv'));
if run_zemax
    writetable(table(fz, core_z, Zflat_S, MX(1,:).', ZS(1,:).', MY(1,:).', ZT(1,:).', ...
        'VariableNames', {'nu_cyc_per_mm','MTF_diff','Zemax_flat','MATLAB_X','Zemax_S','MATLAB_Y','Zemax_T'}), ...
        fullfile(out_dir,'surface001_MATLAB_vs_Zemax.csv'));
    save(fullfile(out_dir,'results.mat'), 'Z','ZS','ZT','MX','MY','err','-append');
end
fprintf('\n결과 저장: %s\n', out_dir);


%% ========================== 로컬 함수 ==========================

function h = synth_surface(Wfun, dx, M)
% 랜덤 표면 생성: 백색잡음 FFT x sqrt(PSD) -> 역FFT
%  Wfun : PSD 함수 [nm^2 mm^2],  dx : 격자 간격 [mm],  M : 격자 점 수
%  분산 = Σ W(f) * df^2  (가우스 랜덤 과정 그대로, RMS 강제 재조정 없음)
H  = fft2(randn(M));
f1 = [0:M/2-1, -M/2:-1].'/(M*dx);                 % FFT 순서의 주파수 축 [cycles/mm]
df = 1/(M*dx);
for c = 1:512:M                                   % 메모리 절약을 위해 열 묶음으로 처리
    id = c:min(c+511, M);
    FR = hypot(f1, f1(id).');                     % 반경 주파수
    H(:,id) = H(:,id) .* (M*df*sqrt(Wfun(FR)));
end
h = ifft2(H, 'symmetric');                        % 실수 표면 [nm]
end

function [ox, oy] = otf_cuts(h_nm, pupil, q, n_lag, L)
% 표면 h의 OTF (x방향, y방향 단면)
%  U = P*exp(i*q*h),  OTF(d) = Σ U(x) conj(U(x+d)) / Σ|P|^2
%  행(열)마다 1D FFT 자기상관을 구해 더함 (2D 전체 FFT보다 메모리 적게 씀)
U  = double(pupil) .* exp(1i*q*h_nm);
np = nnz(pupil);  n = size(U,1);
Px = zeros(1,L);  Py = zeros(L,1);
for r = 1:256:n
    id = r:min(r+255, n);
    Px = Px + sum(abs(fft(U(id,:), L, 2)).^2, 1);   % x방향
    Py = Py + sum(abs(fft(U(:,id), L, 1)).^2, 2);   % y방향
end
cx = conj(ifft(Px));  cy = conj(ifft(Py)).';
ox = cx(1:n_lag)/np;  oy = cy(1:n_lag)/np;
end

function C = hankel_acv(Wfun, fa, fb, r)
% ACV C(r) = 2*pi*∫ W(f) J0(2*pi*f*r) f df   (음수 그대로)
%  J0가 빨리 진동하므로 큰 r에서도 위상 간격 0.25 rad 이하가 되도록 주파수 축 구성
n_lin = min(150000, ceil((fb-fa)*2*pi*max(r)/0.25) + 1);
f = unique([logspace(log10(fa), log10(fb), 20000), linspace(fa, fb, max(n_lin,2))]).';
w = Wfun(f).*f;
C = zeros(size(r));
blk = max(1, floor(4e6/numel(f)));
for k = 1:blk:numel(r)
    id = k:min(k+blk-1, numel(r));
    C(id) = 2*pi*trapz(f, w.*besselj(0, 2*pi*f*r(id)), 1);
end
end

function write_grid_phase_dat(path, phase_rad, dx)
% Zemax Grid Phase DAT 형식으로 저장
%  헤더 : nx ny dx dy 0 0 0
%  본문 : 점마다 "phase dphase/dx dphase/dy d2phase/dxdy 0" -> 미분은 0 (Zemax가 자동 계산)
%  가장자리에 한 칸씩 여유를 둠 (pupil 경계 광선이 격자 밖으로 빠지지 않게)
phase_rad = [phase_rad(:,1), phase_rad, phase_rad(:,end)];
phase_rad = [phase_rad(1,:); phase_rad; phase_rad(end,:)];
[ny, nx] = size(phase_rad);
fid = fopen(path, 'w');
fprintf(fid, '! Reflection phase map: phase = 4*pi*h/lambda [rad]\n');
fprintf(fid, '%d %d %.15g %.15g 0 0 0\n', nx, ny, dx, dx);
d = phase_rad.';
fprintf(fid, '%.10g 0 0 0 0\n', d(:));
fclose(fid);
end

function e = fft_sample_enum(n)
% Zemax FFT MTF 샘플링 값 -> ZOS-API enum
switch n
    case 256,  e = ZOSAPI.Analysis.SampleSizes.S_256x256;
    case 512,  e = ZOSAPI.Analysis.SampleSizes.S_512x512;
    case 1024, e = ZOSAPI.Analysis.SampleSizes.S_1024x1024;
    case 2048, e = ZOSAPI.Analysis.SampleSizes.S_2048x2048;
    case 4096, e = ZOSAPI.Analysis.SampleSizes.S_4096x4096;
    otherwise, error('zemax_fft는 256/512/1024/2048/4096 중 하나');
end
end
