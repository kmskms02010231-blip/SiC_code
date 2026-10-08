

clc; clearvars; close all;


%% 1. 대역 경계 및 PSD 형태 파라미터
D          = 2000;                                   % 유효 구경 지름 D [mm]
lambda_nm  = 650;                                     % 파장 [nm]
lambda     = lambda_nm*1e-6;                          % 파장 [mm]

rho_hp   = 0.005;                                     % Lorentzian PSD knee 주파수 (MASTER_real의 rho_hp) [cycles/mm]
p        = 2.60;                                      % Lorentzian PSD 고주파 감쇠 지수 (MASTER_real의 p) [-]

f_min = 1/D;                                          % 유한 구경이 표현하는 최저 표면주파수 [cycles/mm]
f_LM  = 7/D;                                          % LSF-MSF 경계 [cycles/mm]
f_MH  = sqrt(10/(lambda*1378.014));                   % MSF-HSF 경계 [cycles/mm]
f_max = 1/(2*0.325e-3);                               % HSF 최고주파수 [cycles/mm]

%% 2. LSF 2D 격자 및 정규화 계수 h_L
widthXY = 40000;                                      % [공간 영역] 격자 가상 도메인 폭 [mm]
ss      = 2048;                                       % 격자 한 변 점 개수 [-]
A_area  = widthXY^2;                                  % [공간 영역] 전체 면적 [mm^2]

u = ((0:ss-1)-ss/2)/widthXY;                          % [주파수 영역] 1D 축 [cycles/mm]
[U, V] = meshgrid(u, u);

rho = sqrt(U.^2 + V.^2);                              % [주파수 영역] 2D 등방 반경주파수 [cycles/mm]

base_lorentz = 1./(1+(rho/rho_hp).^p);
mask_L = (rho >= f_min) & (rho <= f_LM);               % LSF 대역 마스크

h_L = sum(base_lorentz(mask_L));                       % LSF 대역만의 정규화 계수 (식 22의 h)

df          = 1/widthXY;                              % 격자 주파수 해상도 [cycles/mm]
n_bins      = round((f_LM - f_min)/df);               % LSF 대역의 반경 해상도 한계 — 더 키우면 빈 bin(NaN) 생김

f_edges_M = logspace(log10(f_LM), log10(f_MH), n_bins+1);
f_axis_M  = [f_LM, sqrt(f_edges_M(1:end-1).*f_edges_M(2:end)), f_MH];   % 경계점 + 로그축 bin 중심

f_edges_H = logspace(log10(f_MH), log10(f_max), n_bins+1);
f_axis_H  = [f_MH, sqrt(f_edges_H(1:end-1).*f_edges_H(2:end)), f_max];

%% 3. [HSF 목표 RMS별 σ_L-σ_M 가능 조합 경계

sigma_H_targets = [1, 1.5, 2];                         % 고정해서 반복할 HSF RMS 값들 [nm]

sigma_L_grid = linspace(1, 30, 40);                    % 탐색할 σ_L 범위 [nm]
sigma_M_grid = linspace(1, 30, 40);                    % 탐색할 σ_M 범위 [nm]
[SIGMA_L, SIGMA_M] = meshgrid(sigma_L_grid, sigma_M_grid);

p_cap = 1e4;                                            % p_M 탐색 상한 — 이산화된 f_axis_M의 첫 구간
                                                        % 때문에 매우 큰 p에서도 분산이 0으로 가지 않고
                                                        % 바닥값에 수렴하므로, 그 이상은 사실상 무한대로 취급
sigma_H_max = zeros(size(SIGMA_L));
for k = 1:numel(SIGMA_L)
    S_L_edge_k = (SIGMA_L(k)^2*A_area/h_L) * (1/(1+(f_LM/rho_hp)^p));

    S_M_of_p_k   = @(pp) S_L_edge_k*(f_axis_M/f_LM).^(-pp);
    var_M_of_p_k = @(pp) 2*pi*trapz(f_axis_M, S_M_of_p_k(pp).*f_axis_M);
    if var_M_of_p_k(0) <= SIGMA_M(k)^2
        p_M_k = 0;                                     % σ_M 목표 자체가 이미 감소 조건에서 불가능(너무 큼)
    elseif var_M_of_p_k(p_cap) >= SIGMA_M(k)^2
        p_M_k = p_cap;                                 % σ_M 목표가 이산화 바닥값보다 작음(너무 작음) -> 최대 감쇠로 근사
    else
        p_M_k = fzero(@(pp) var_M_of_p_k(pp) - SIGMA_M(k)^2, [0, p_cap]);
    end
    S_M_edge_k = S_L_edge_k*(f_MH/f_LM)^(-p_M_k);

    sigma_H_max(k) = sqrt(2*pi*trapz(f_axis_H, S_M_edge_k*f_axis_H));   % p_H=0(평탄)일 때 최대 HSF RMS
end

figure(1);
contour(SIGMA_L, SIGMA_M, sigma_H_max, sigma_H_targets, 'ShowText', 'on', 'LineWidth', 2.0);
set(gca, 'FontSize', 13, 'LineWidth', 1.0, 'Box', 'on');
grid on;
xlabel('\sigma_L [nm]');
ylabel('\sigma_M [nm]');
title('HSF RMS별 실현 가능한 \sigma_L-\sigma_M 경계 (선 위쪽=\sigma_M/\sigma_L 큰 쪽이 가능 영역)', 'FontWeight', 'normal');

%% 4. 관심 변수 — 대역별 목표 RMS
sigma_L_nm = 15;                                       % LSF 대역 목표 RMS [nm]
sigma_M_nm = 1;                                       % MSF 대역 목표 RMS [nm] (ripple 포함 총량)
sigma_H_nm = 1;                                        % HSF 대역 목표 RMS [nm]

%% 4-1. Ripple 파라미터 및 MSF 배경 RMS 보정

% [주파수 fr(cycles/mm), RMS sigma_r(nm)] — radial(circular) ripple 목록
ripple_radial = [ ...
    0.21, 5; ...
    0.33, 0; ...
    0.40, 0];

% [주파수 fr(cycles/mm), RMS sigma_r(nm), 방향각 angle_deg(x축 기준)] — raster(x+y 혼합) ripple 목록
% angle_deg = 0  -> 순수 x-raster, angle_deg = 90 -> 순수 y-raster(투영 시 x축 STF에 기여 없음)
ripple_raster = [ ...
    1.01, 5,  0];

sigma_ripple_radial2 = sum(ripple_radial(:,2).^2);      % radial ripple 총 분산 [nm^2]
sigma_ripple_raster2 = sum(ripple_raster(:,2).^2);      % raster ripple 총 분산 [nm^2]
sigma_ripple2_total  = sigma_ripple_radial2 + sigma_ripple_raster2;

sigma_M_bg_nm = sqrt(sigma_M_nm^2 - sigma_ripple2_total);   % MSF 배경(연속 PSD)이 채울 RMS

fprintf('\n[MSF RMS 배분]\n');
fprintf('  지정 sigma_M   = %7.3f nm\n', sigma_M_nm);
fprintf('  배경(연속 PSD) = %7.3f nm\n', sigma_M_bg_nm);
fprintf('  Ripple         = %7.3f nm   (분산 %.4f nm^2)\n', sqrt(sigma_ripple2_total), sigma_ripple2_total);

%% 5. LSF PSD 정규화
sPSD_2D_L = (sigma_L_nm^2*A_area/h_L) * base_lorentz .* mask_L;    % LSF 2D PSD [nm^2 mm^2]

%% 6. LSF 2D PSD -> 1D PSD
f_edges_L   = linspace(f_min, f_LM, n_bins+1);
f_axis_L    = (f_edges_L(1:end-1)+f_edges_L(2:end))/2;

% 방사형 데이터를 구간별 평균 -> 1D로 변환
S_L = zeros(1, n_bins);
for k = 1:n_bins
    ring_idx = (rho >= f_edges_L(k)) & (rho < f_edges_L(k+1));
    S_L(k) = mean(sPSD_2D_L(ring_idx));
end

S_L_at   = @(f) (sigma_L_nm^2*A_area/h_L) * (1./(1+(f/rho_hp).^p));
S_L_min  = S_L_at(f_min);
S_L_edge = S_L_at(f_LM);                              % LSF PSD를 f_LM에서 평가한 값 -> MSF 시작값
f_axis_L = [f_min, f_axis_L, f_LM];
S_L      = [S_L_min, S_L, S_L_edge];

sig2_L_nm2_check = 2*pi*trapz(f_axis_L, S_L.*f_axis_L);

fprintf('\n[PSD 정규화 검증]  (분산/목표^2)\n');
fprintf('  LSF       : %.2f\n', sig2_L_nm2_check/sigma_L_nm^2);

%% 7. MSF 그래프 생성

S_M_of_p    = @(p) S_L_edge*(f_axis_M/f_LM).^(-p);
var_M_of_p  = @(p) 2*pi*trapz(f_axis_M, S_M_of_p(p).*f_axis_M);
p_M         = fzero(@(p) var_M_of_p(p) - sigma_M_nm^2, p);   % 전체 sigma_M_nm 기준 -> 모양은 ripple과 무관
S_M_natural = S_M_of_p(p_M);

trans_frac = 0.15;                                       % 전환 구간 폭 = MSF 대역 log-range의 15%(표시/전환용)
f_trans    = f_LM*(f_MH/f_LM)^trans_frac;                 % 전환이 끝나는 주파수(그 이후는 alpha_M로 평탄)
t_of_f     = @(f) min(max(log(f/f_LM)/log(f_trans/f_LM), 0), 1);
g_of_f     = @(f, a) a + (1-a)*0.5*(1+cos(pi*t_of_f(f)));  % f_LM에서 g=1(연속), f_trans부터 g=a(평탄)
                                                            % (지수형은 감쇠가 대역 전체로 퍼져 slope 재적합과
                                                            %  똑같이 alpha가 극단으로 작아지며 역전이 재발함 -> 기각)

var_M_of_a = @(a) 2*pi*trapz(f_axis_M, g_of_f(f_axis_M,a).*S_M_natural.*f_axis_M);
alpha_M0   = sigma_M_bg_nm^2/sigma_M_nm^2;                 % 초기값: 상수배로 근사했을 때의 비율
alpha_M    = fzero(@(a) var_M_of_a(a) - sigma_M_bg_nm^2, alpha_M0);

S_M      = g_of_f(f_axis_M, alpha_M) .* S_M_natural;       % f_LM에서 LSF와 정확히 연속, 이후 완만히 감쇠
S_M_edge = S_M(end);                                        % 보정된 MSF PSD의 f_MH 값 -> HSF 시작값

sig2_M_nm2_check = 2*pi*trapz(f_axis_M, S_M.*f_axis_M);
fprintf('  MSF(배경) : %.2f   (기울기 p_M = %.2f, 평탄구간 진폭계수 alpha_M = %.3f)\n', ...
    sig2_M_nm2_check/sigma_M_bg_nm^2, p_M, alpha_M);

%% 8. HSF 그래프 생성
S_H_of_p   = @(p) S_M_edge*(f_axis_H/f_MH).^(-p);
var_H_of_p = @(p) 2*pi*trapz(f_axis_H, S_H_of_p(p).*f_axis_H);
p_H = fzero(@(p) var_H_of_p(p) - sigma_H_nm^2, p);
S_H = S_H_of_p(p_H);

sig2_H_nm2_check = 2*pi*trapz(f_axis_H, S_H.*f_axis_H);
fprintf('  HSF       : %.2f   (기울기 p_H = %.2f)\n', sig2_H_nm2_check/sigma_H_nm^2, p_H);

%% 9. PSD 그래프

c_L = [0.00 0.45 0.70];
c_M = [0.10 0.60 0.30];
c_H = [0.85 0.33 0.10];
c_R = [0.55 0.00 0.55];                                % ripple 표기 색

figure(2);
loglog(f_axis_L, S_L, '-', 'Color', c_L, 'LineWidth', 2.0); hold on;
loglog(f_axis_M, S_M, '-', 'Color', c_M, 'LineWidth', 2.0);
loglog(f_axis_H, S_H, '-', 'Color', c_H, 'LineWidth', 2.0);

% 임의의 f에서 배경(연속 PSD) 값을 돌려주는 함수 — 범프가 어느 대역 위에 있든
% 그 대역의 배경 곡선에서 솟아났다가 다시 내려오도록 기준선으로 사용한다.
bg_at = @(f) (f<=f_LM).*S_L_at(f) + ...
             (f>f_LM & f<=f_MH).*(g_of_f(f,alpha_M).*S_L_edge.*(f/f_LM).^(-p_M)) + ...
             (f>f_MH).*(S_M_edge*(f/f_MH).^(-p_H));

rel_width      = 0.10;                                   % 범프 폭 = fr의 10%(표시용 값, 계산에는 무관)
ripple_fr_all  = [ripple_radial(:,1); ripple_raster(:,1)];
ripple_sig_all = [ripple_radial(:,2); ripple_raster(:,2)];
ripple_fr_nz   = ripple_fr_all(ripple_sig_all > 0);
ripple_sig_nz  = ripple_sig_all(ripple_sig_all > 0);

for kk = 1:numel(ripple_fr_nz)
    fr_k   = ripple_fr_nz(kk);
    sig2_k = ripple_sig_nz(kk)^2;
    df_k   = rel_width * fr_k;
    f_lo   = fr_k - df_k/2;
    f_hi   = fr_k + df_k/2;

    % 삼각형 범프의 넓이(2*pi*fr_k 가중, 밑변 df_k, 높이 h_peak) = sig2_k 가 되도록 역산:
    %   2*pi*fr_k * (1/2 * df_k * h_peak) = sig2_k  ->  h_peak = sig2_k / (pi*fr_k*df_k)
    h_peak = sig2_k / (pi*fr_k*df_k);
    f_bump = [f_lo, fr_k, f_hi];
    S_bump = [bg_at(f_lo), bg_at(fr_k)+h_peak, bg_at(f_hi)];
    if kk == 1
        plot(f_bump, S_bump, '-', 'Color', c_R, 'LineWidth', 2.0);
    else
        plot(f_bump, S_bump, '-', 'Color', c_R, 'LineWidth', 2.0, 'HandleVisibility', 'off');
    end
end

xline(f_LM, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.1, 'HandleVisibility', 'off');
xline(f_MH, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.1, 'HandleVisibility', 'off');
set(gca, 'FontSize', 13, 'LineWidth', 1.0, 'Box', 'on');
grid on;
xlabel('Surface frequency [cycles/mm]');
ylabel('PSD [nm^2\cdotmm^2]');
legend({'LSF','MSF(배경)','HSF','Ripple'}, 'Location', 'southwest', 'Box', 'off', 'FontSize', 12);
title('Surface PSD (연속 성분 + ripple)', 'FontWeight', 'normal');

% PSD 함수화 (임의의 f에서 값을 뽑을 수 있게) — 10-1/10-2에서 공통으로 사용
S_L_fun = S_L_at;                                                     % LSF PSD 함수
S_M_fun = @(f) g_of_f(f,alpha_M) .* S_L_edge.*(f/f_LM).^(-p_M);       % MSF 배경 PSD 함수 (진폭 보정 포함)
S_H_fun = @(f) S_M_edge*(f/f_MH).^(-p_H);                             % HSF PSD 함수
opts    = {'RelTol',1e-4,'AbsTol',1e-4};

% r축: 로그축으로 r=0 근방까지 촘촘하게
% 최소값은 f_max 가 만드는 상관길이(~1/f_max)보다 충분히 작게 잡아야 함
N_r   = 600;
r_min = 1/(50*f_max);
r_max = D;
r     = [0, logspace(log10(r_min), log10(r_max), N_r)];

%% 10-1. PSD -> ACV, FT(2D FFT, Wiener-Khinchin) 방식 — LSF만 검증용으로 계산
% Wiener-Khinchin: ACV = IFFT(PSD).
% 2절에서 만든 sPSD_2D_L 이용 -> FFT로 ACV 계산
% MSF/HSF는 f_max(~1538 cycles/mm)까지 담으려면 격자 한 변이 수만~수십만 점이 되어야
% 해서(메모리상 불가능) 여기서는 생략

C_2D_L_FFT = real((ss/widthXY)^2 * fftshift(ifft2(ifftshift(sPSD_2D_L))));  % [nm^2], 2D ACV

x_fft = ((0:ss-1)-ss/2)*(widthXY/ss);                 % [공간 영역] 1D 축 [mm], 격자와 같은 컨벤션
[X_fft, Y_fft] = meshgrid(x_fft, x_fft);
rho_fft = sqrt(X_fft.^2 + Y_fft.^2);

N_r_fft     = 30;                                     % 격자 셀 크기(dx=widthXY/ss~20mm) 대비 적당한 bin 수
r_edges_fft = linspace(0, D, N_r_fft+1);
r_axis_fft  = (r_edges_fft(1:end-1)+r_edges_fft(2:end))/2;

for k = 1:N_r_fft
    ring_idx   = (rho_fft >= r_edges_fft(k)) & (rho_fft < r_edges_fft(k+1));
    C_L_FFT(k) = mean(C_2D_L_FFT(ring_idx));
end

%% 10-2. PSD -> ACV, 한켈 변환 방식 C(r) = 2*pi*∫ S(f)*J0(2*pi*f*r)*f df (대역별 선형중첩)

for k = 1:numel(r)
    C_L(k) = 2*pi*integral(@(f) S_L_fun(f).*besselj(0,2*pi*f*r(k)).*f, f_min, f_LM, opts{:});
    C_M(k) = 2*pi*integral(@(f) S_M_fun(f).*besselj(0,2*pi*f*r(k)).*f, f_LM,  f_MH, opts{:});
    C_H(k) = 2*pi*integral(@(f) S_H_fun(f).*besselj(0,2*pi*f*r(k)).*f, f_MH,  f_max, opts{:});
end

fprintf('\n[ACV 검증]  (r=0에서 분산과 일치해야 함)\n');
fprintf('  LSF       : %9.2f   (목표 %8.2f)\n', C_L(1), sigma_L_nm^2);
fprintf('  MSF(배경) : %9.2f   (목표 %8.2f)\n', C_M(1), sigma_M_bg_nm^2);
fprintf('  HSF       : %9.2f   (목표 %8.2f)\n', C_H(1), sigma_H_nm^2);

% FT(10-1) vs 한켈(10-2) 비교 — 같은 r_axis_fft 위치에서 한켈 결과를 보간해 대조.
% 값이 0을 넘나드는 구간이 있어 점별 상대오차 대신 전체 분산(C_L(1)) 대비 오차로 비교.
C_L_hankel_at_fft = interp1(r, C_L, r_axis_fft, 'pchip');
err_rms_pct = sqrt(mean((C_L_FFT - C_L_hankel_at_fft).^2)) / C_L(1) * 100;
fprintf('\n[FT vs 한켈 검증, LSF]  RMS 오차 / 전체분산 = %.3f%%\n', err_rms_pct);

%% 10-3. Ripple ACV — 해석적(closed-form) 계산, 수치 적분/변환 불필요
% radial: C(r)  = sigma_r^2 * J0(2*pi*fr*r)                (등방 앙상블 근사, 조리개 중심·다주기 가정)
% raster: C(r)  = sigma_r^2 * cos(2*pi*(fr*cosd(angle))*r) (x축으로 정확히 사영, 근사 아님)

fr_radial   = ripple_radial(:,1);                                 % [N_radial x 1]
sig2_radial = ripple_radial(:,2).^2;
C_ripple_radial = sum(sig2_radial .* besselj(0, 2*pi*fr_radial*r), 1);   % [nm^2], 1 x numel(r)

fr_raster_eff = ripple_raster(:,1) .* cosd(ripple_raster(:,3));   % x축 사영 유효 주파수 [N_raster x 1]
sig2_raster   = ripple_raster(:,2).^2;
C_ripple_raster = sum(sig2_raster .* cos(2*pi*fr_raster_eff*r), 1);      % [nm^2], 1 x numel(r)

C_ripple_total = C_ripple_radial + C_ripple_raster;

fprintf('\n[Ripple ACV 검증]  r=0에서 C_ripple = %.2f  (목표 %.2f)\n', C_ripple_total(1), sigma_ripple2_total);

% symlog 방식 x축 준비
r0    = r_min;                    % 선형<->로그 전환 스케일
x_sym = asinh(r/r0);              % r=0 -> 0, r>>r0에서는 로그축처럼 동작

tick_r     = [0, 1e-4, 1e-2, 1, 1e2, r_max];
tick_x     = asinh(tick_r/r0);
tick_label = arrayfun(@(v) sprintf('%g', v), tick_r, 'UniformOutput', false);

%% 11. 대역별 신뢰구간 판정 및 결합

idx_trust_L = find(C_L <= 0, 1, 'first') - 1;   % 각 대역이 처음 0 밑으로 내려가기 직전까지
idx_trust_M = find(C_M <= 0, 1, 'first') - 1;
idx_trust_H = find(C_H <= 0, 1, 'first') - 1;

C_L_trust = zeros(size(r)); C_L_trust(1:idx_trust_L) = C_L(1:idx_trust_L);
C_M_trust = zeros(size(r)); C_M_trust(1:idx_trust_M) = C_M(1:idx_trust_M);
C_H_trust = zeros(size(r)); C_H_trust(1:idx_trust_H) = C_H(1:idx_trust_H);

C_total_surf = C_L_trust + C_M_trust + C_H_trust;               % 표면(연속 PSD, ripple 제외) 전체 ACV
C_total_raw  = C_L + C_M + C_H + C_ripple_total;                 % 신뢰구간 처리 이전의 원본(검증용)
C_total      = C_total_surf + C_ripple_total;                    % 공식 채택 결과(ripple 포함) — 이후 STF 계산까지 이 변수를 사용

fprintf('\n[대역별 신뢰구간]  (한켈 적분을 신뢰하는 r 상한)\n');
fprintf('  LSF       : r <= %10.2f mm\n', r(idx_trust_L));
fprintf('  MSF(배경) : r <= %10.2f mm\n', r(idx_trust_M));
fprintf('  HSF       : r <= %10.6f mm\n', r(idx_trust_H));

%% 12. 대역별 ACV 그래프 (fig 3) — ripple 포함
figure(3);
plot(x_sym, C_L_trust, '-', 'Color', c_L, 'LineWidth', 2.0); hold on;
plot(x_sym, C_M_trust, '-', 'Color', c_M, 'LineWidth', 2.0);
plot(x_sym, C_H_trust, '-', 'Color', c_H, 'LineWidth', 2.0);
% plot(x_sym, C_ripple_total, '-', 'Color', c_R, 'LineWidth', 2.0);
set(gca, 'FontSize', 13, 'LineWidth', 1.0, 'Box', 'on', 'XTick', tick_x, 'XTickLabel', tick_label);
grid on;
xlabel('Lag r [mm] (symlog)');
ylabel('ACV [nm^2]');
% legend({'LSF','MSF(배경)','HSF','Ripple'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
legend({'LSF','MSF(배경)','HSF'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
title('Surface ACV by band (ripple 포함)', 'FontWeight', 'normal');

%% 13. 전체 ACV 그래프 (fig 4, 정규화된 전체 상관함수)
rho_ACV = C_total / C_total(1);

figure(5);
plot(x_sym, rho_ACV, '-', 'Color', [0.20 0.20 0.20], 'LineWidth', 2.0);
set(gca, 'FontSize', 13, 'LineWidth', 1.0, 'Box', 'on', 'XTick', tick_x, 'XTickLabel', tick_label);
grid on;
xlabel('Lag r [mm] (symlog)');
ylabel('Normalized ACV \rho(r)');
title('Surface ACV (ripple 포함)', 'FontWeight', 'normal');

%% 14. [검증용] 신뢰구간 처리 전/후 비교 그래프
figure(4);
plot(x_sym, C_total_raw, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.3); hold on;
plot(x_sym, C_total,     '-',  'Color', [0.20 0.20 0.20], 'LineWidth', 2.0);
xline(asinh(r(idx_trust_L)/r0), ':', 'Color', c_L, 'LineWidth', 1.1, 'HandleVisibility', 'off');
xline(asinh(r(idx_trust_M)/r0), ':', 'Color', c_M, 'LineWidth', 1.1, 'HandleVisibility', 'off');
xline(asinh(r(idx_trust_H)/r0), ':', 'Color', c_H, 'LineWidth', 1.1, 'HandleVisibility', 'off');
set(gca, 'FontSize', 13, 'LineWidth', 1.0, 'Box', 'on', 'XTick', tick_x, 'XTickLabel', tick_label);
grid on;
xlabel('Lag r [mm] (symlog)');
ylabel('ACV [nm^2]');
legend({'Before trust (raw)','After trust (adopted)'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
title('Surface ACV (before vs after trust, ripple 포함)', 'FontWeight', 'normal');

%% 15. 광학 파라미터 및 ACV -> STF 계산
% STF(nu) = exp(-0.5*q^2*Dh(r)), Dh(r) = 2*sigma_total^2*(1-rho_ACV(r))
% nu = r/(lambda*f_eff): 표면(pupil) lag r을 이미지 공간 주파수 nu로 재변환
% Dh_total = Dh_surf + Dh_ripple 이므로 STF_total = STF_surf.*STF_ripple (곱셈 분해, 상호 무상관 가정)

Fnum      = 13.91049;                          % F-number
theta_deg = 0;                                 % 입사각 [deg]
f_eff     = Fnum*D;                            % 유효 초점거리 [mm]
nu_c      = 1/(lambda*Fnum);                   % 회절 컷오프 주파수 [cycles/mm] (이미지 공간)
q_1pnm    = 4*pi*cosd(theta_deg)/lambda_nm;    % 표면오차->위상 변환 상수 [1/nm]

nu      = r/(lambda*f_eff);                    % r -> nu 매핑
nu_norm = nu/nu_c;                             % 정규화 주파수 nu/nu_c (0~1이 회절 통과 대역)

Dh_L_nm2      = 2*(C_L_trust(1)      - C_L_trust);        % 대역별 높이 구조함수 [nm^2]
Dh_M_nm2      = 2*(C_M_trust(1)      - C_M_trust);
Dh_H_nm2      = 2*(C_H_trust(1)      - C_H_trust);
Dh_surf_nm2   = 2*(C_total_surf(1)   - C_total_surf);      % 표면(ripple 제외) 전체 구조함수
Dh_ripple_nm2 = 2*(C_ripple_total(1) - C_ripple_total);    % ripple만의 구조함수
Dh_nm2        = 2*(C_total(1)        - C_total);           % 전체(표면+ripple) 구조함수

STF_L      = exp(-0.5*q_1pnm^2 .* Dh_L_nm2);       % 대역별 STF(nu)
STF_M      = exp(-0.5*q_1pnm^2 .* Dh_M_nm2);
STF_H      = exp(-0.5*q_1pnm^2 .* Dh_H_nm2);
STF_surf   = exp(-0.5*q_1pnm^2 .* Dh_surf_nm2);     % 표면만(ripple 제외) — 참고/검증용
STF_ripple = exp(-0.5*q_1pnm^2 .* Dh_ripple_nm2);   % ripple만 — 참고/검증용
STF        = exp(-0.5*q_1pnm^2 .* Dh_nm2);          % 최종 채택 STF(nu) — 표면+ripple 전체

%% 16. 대역별 STF 그래프 (fig 6)

figure(6);
plot(nu, STF_L, '-', 'Color', c_L, 'LineWidth', 2.0); hold on;
plot(nu, STF_M, '-', 'Color', c_M, 'LineWidth', 2.0);
plot(nu, STF_H, '-', 'Color', c_H, 'LineWidth', 2.0);
plot(nu, STF_ripple, '-', 'Color', c_R, 'LineWidth', 2.0);
set(gca, 'FontSize', 13, 'LineWidth', 1.0, 'Box', 'on');
grid on;
xlim([0 max(nu)]);
ylim([0 1.02]);
xlabel('Normalized frequency \nu/\nu_c');
ylabel('STF(\nu)');
legend({'LSF','MSF(배경)','HSF', 'Ripple'}, 'Location', 'southwest', 'Box', 'off', 'FontSize', 16);
title('Surface Transfer Function by band', 'FontWeight', 'normal');

%% 17. 전체 STF 그래프 (fig 7)

figure(7);
plot(nu, STF_surf, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); hold on;
plot(nu, STF,      '-',  'Color', [0.00 0.30 0.60], 'LineWidth', 2.0);
set(gca, 'FontSize', 13, 'LineWidth', 1.0, 'Box', 'on');
grid on;
ylim([0 1.05]);
xlim([0 max(nu)]);
xlabel('Image frequency \nu [cycles/mm]');
ylabel('STF(\nu)');
legend({'표면만(ripple 제외)','표면+ripple'}, 'Location', 'southwest', 'Box', 'off', 'FontSize', 12);
title('Surface Transfer Function', 'FontWeight', 'normal');

%% 18. 회절 MTF 및 표면오차(STF) 반영 시스템 MTF (fig 8)

nu_norm_c = min(nu_norm, 1);
MTF_diff  = (2/pi)*(acos(nu_norm_c) - nu_norm_c.*sqrt(1-nu_norm_c.^2));   % 회절한계 MTF
MTF_surf  = MTF_diff .* STF_surf;                                         % 표면오차(ripple 제외) 반영 MTF
MTF_total = MTF_diff .* STF;                                              % 표면오차+ripple 반영 시스템 MTF

nu_nyq = 71.43;                                % 픽셀 피치 7um 기준 Nyquist 주파수 [cycles/mm]

figure(8);
plot(nu, MTF_diff,  '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5); hold on;
xline(nu_nyq, '--r')
plot(nu, MTF_surf,  '-',  'Color', [0.30 0.30 0.30], 'LineWidth', 1.5);
plot(nu, MTF_total, '-',  'Color', [0.70 0.10 0.10], 'LineWidth', 2.0);
set(gca, 'FontSize', 13, 'LineWidth', 1.0, 'Box', 'on');
grid on;
xlim([0 max(nu)])
ylim([0 1]);
xlabel('Image frequency \nu [cycles/mm]');
ylabel('MTF');

legend({'Diffraction-limited','Nyquist freq.', 'With surface error + ripple', 'With surface error only' }, ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', 12);
title('System MTF (Diffraction \times STF)', 'FontWeight', 'normal');

%% 19. Nyquist 주파수에서의 MTF 판정


[~, idx_nyq] = min(abs(nu - nu_nyq));
nu_at_nyq  = nu(idx_nyq);
MTF_at_nyq = MTF_total(idx_nyq);

fprintf('\n[Nyquist MTF 판정]  (ripple 포함)\n');
fprintf('  nu_Nyquist         = %6.2f cycles/mm\n', nu_nyq);
fprintf('  nu(가장 가까운 값) = %6.2f cycles/mm\n', nu_at_nyq);
fprintf('  MTF                = %6.4f\n', MTF_at_nyq);

if MTF_at_nyq >= 0.2
    fprintf('  -> PASS (MTF >= 0.2)\n\n');
else
    fprintf('  -> FAIL (MTF < 0.2)\n\n');
end
