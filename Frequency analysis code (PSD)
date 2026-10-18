%% Four PSD files and correlation analysis
% Power-step maps and variable-by-variable correlation / phase matrices of the two parallel channels
% (19 power steps, 900 s analysis window each).
%
% CHANGES vs v4.4   (every figure and every table of v4.4 is kept, the columns of Band_table.csv / Step_table.csv stay in place;
%                    what changes are the LEVELS of the spectra, see item A)
%   A. NORMALISATION.  v4.4 computed the spectra in absolute units (flow (mL/s)^2/Hz, temperature K^2/Hz) and divided only the
%      amplitude A by the mean flow.  Now (normMode = 'mean') EVERY variable is divided by its own mean over the 900 s window
%      BEFORE any spectrum is computed:   x* = normScale * x / <x>   (normScale = 100: percent of the window mean).
%      <x> = mean|q| for the flow (qRef), the mean in kelvin for the temperature (the records are in degC: tempOffset = 273.15).
%      The PSD of x* is a relative PSD in %^2/Hz (colour bars: dB re 1 %^2/Hz), and sqrt(integral of the PSD over a band) is the
%      standard deviation / mean of the band-limited signal in percent, i.e. the amplitude A of v4.4, which is therefore unchanged
%      (A_slow, A_mid, A_fast, A_high, Atot).  The levels of different steps, channels and variables are now comparable.
%      Scale-free, hence unchanged: f_p of each channel, share, lfShare, prom1 / prom2, coherence, phase, rho, gam, A, the
%      correlation matrices, the regime rule R1.
%      Changed in value: pk_dB, lf_dB, lev_*_dB (levels), the PSD maps Fig_01 ... Fig_04, FigS2, FigS3 (common and differential
%      mode of x*).  Built on the SUM of the two channel spectra (S11+S22: fpc, fpk_b, prom_b; and the ratios sd_b, SumDiff_dB) the
%      channels now weigh by their RELATIVE fluctuation (v4.4: by the absolute one): identical for equal channels, slightly
%      different when the channels have different mean flows.  normMode = 'none' reproduces v4.4 exactly.
%   B. New table tables\Cv_table.csv: standard deviation / mean (percent) of EVERY variable (flow, T_in ... T_w5, P_in, P_out, dP of
%      both channels) in every window (10 Hz copy, linear trend removed) + a console table.  Reference = window mean: |q| for the
%      flow, kelvin for the temperatures (tempOffset), ABSOLUTE pressure for P_in and P_out (presOffset: the records are gauge
%      pressure), dP itself.
%   C. FreqCorr_Results.mat also holds normMode, normScale, tempOffset, presOffset, psdUnitQ, psdUnitT, tRef and Cv; the 3-D scripts
%      (Sec_4_1_2_OriginData_v4, Sec_4_1_2_PSD3D_v3) take the unit of the colour bar from it and smooth the landscape by 0.005
%      decades.  The 2-D maps of this script stay unsmoothed (smoothBins = 1).
%   D. Output folder Fig_4_1_2_v4_5.
%
% CHANGES of v4.4 (kept; vs v4.3)   (every figure of v4.3 is kept; Band_table.csv and Step_table.csv keep all their columns in place, new
%                    columns are added at the END; Regime_rule_check.csv is rebuilt for rule R1)
%   A. Regime rule R1 replaces the v4.2 rule in the console table, in tables\Regime_rule_check.csv and in the last column
%      (ruleId) of Step_table.csv.  The v4.2 rule looked only at the global spectral maximum and called the highest-power step
%      (Mode 3) Stable.  R1 uses the band metrics of v4.3 and is applied in this order:
%        (1) rho_slow < ruleRhoSlow (0)                    -> Mode 3   (the channels move in anti-phase in the slow band)
%        (2) A_slow >= ruleASlow (5 %), larger channel     -> Mode 1   (slow pulsation)
%        (3) f_p >= ruleFpFast (0.4 Hz) in BOTH channels   -> Mode 2   (oscillation near 0.8 Hz)
%        (4) otherwise                                     -> Stable
%      The console prints, for every step, the quantities that decide it and, for every test, how far apart the nearest steps of
%      the two classes are (separation of the classes); Regime_rule_check.csv holds the same numbers.
%   B. Fourth band "high" = 2-100 Hz (fHigh), computed from the RAW 200 Hz flow record with the same Welch segmentation as the
%      other bands (150 s Hann, 50 % overlap, bin width 0.00244 Hz): A1_high, A2_high, rho_high, gam_high, ph_high.  Together with
%      the three bands of v4.3 it closes the variance budget of the 0.01-100 Hz record.  New: the share of the 0.01-100 Hz
%      variance above 2 Hz and above 10 Hz (hfShare2, hfShare10, percent) and psdDiff = largest difference (dB) between the PSD of
%      the 200 Hz record and the PSD of its 10 Hz copy inside 0.01-2 Hz (the 10 Hz copy is what all other band metrics use).
%      File: last columns of tables\Band_table.csv + a second console table.  The 0.01-2 Hz values A_slow, A_mid, A_fast and
%      Atot are band-limited amplitudes, not the total fluctuation of the flow.
%   C. FigS4_band_metrics: fourth panel for the high band (bandFigHigh); the y-axis of the amplitude row no longer jumps to 1000
%      because of the start-up step (0.1 ... 200); shorter y-label.
%   D. FigS2: the dotted vertical line is now the f_p threshold of the Mode 2 test of R1 (0.4 Hz; was f_c = 0.1 Hz of the v4.2 rule).
%
% CHANGES of v4.3 (kept)
%   A. Band metrics.  The flow spectrum 0.01-2 Hz is split into three bands (bandEdges: slow / mid / fast).  For every
%      step and band: the RMS of the flow fluctuation of each channel relative to the mean |q| (A1, A2); the correlation
%      coefficient between the two channels in that band (rho: +1 in phase, -1 anti-phase) with the modulus and phase
%      of the band-integrated cross-spectrum; and the maximum of the common spectrum S11+S22 inside the band (frequency,
%      edge flag, level, coherence / phase / S_sum/S_diff at that frequency, and the maximum against the band median).
%      These do not depend on where the global maximum of the spectrum lies, so a step whose global maximum sits at the
%      lower band edge (e.g. the highest-power step) is still described.  File: tables\Band_table.csv + console table.
%   B. New figure FigS3_channel_coupling: coherence, phase and common-mode / differential-mode spectra between the two
%      channels for one representative step per regime (makeCohSpecFig).
%   C. New figure FigS4_band_metrics: A and rho of the three bands against the 19 power steps (makeBandFig).
%   D. Flow10Hz_windows.mat (exportWindows): the 10 Hz flow and dP of every 900 s window, small file, for checking or
%      redesigning the regime rule outside this script.
%   E. (v4.3 only; replaced by rule R1 in v4.4) The regime rule was the one of v4.2.  It looked at the global maximum only,
%      so it called a step Stable when the maximum lay at the lower band edge in BOTH channels (this happens at the
%      highest power).
%
% CHANGES of v4.2 (kept)
%   1. The PSD maps are drawn WITHOUT frequency smoothing (smoothBins = 1) on a finer frequency grid (nFq = 2400).
%      In v4.1 the Gaussian smoothing along frequency (sigma about 0.015 Hz at 0.8 Hz) was wider than a spectral
%      line (width about 0.01 Hz), so a narrow line was flattened on the map (a line 10 dB above the slow-frequency
%      level was shown about 1 dB BELOW it) and the slow frequencies looked dominant even where the maximum of
%      the spectrum, the value f_p is taken from, is a line near 0.8 Hz.  f_p, the tables and the matrices never
%      used the smoothed map and are unchanged.
%   2. New figure FigS2_regime_spectra (makeSpecFig): the raw flow PSD (200 Hz record) and the temperature PSD of
%      one representative step per regime, both channels, with f_p marked on the curve.  It shows directly where
%      the maximum of the spectrum lies.  specSteps = [] takes the middle step of every regime.
%   3. Step_table.csv has new columns at the END (the old columns are unchanged), per channel 1 / 2:
%      pk_dB (level of the maximum), lf_dB (highest level below fLow_Hz), pkMinusLf_dB (= pk_dB - lf_dB, 0 when the
%      maximum itself lies below fLow_Hz), prom_dB (maximum against the median of its neighbourhood 0.05-0.30 Hz
%      away), lfShare (share of the 0.01-2 Hz variance below fLow_Hz), coh12_fp (coherence between the channels at
%      that channel's f_p); and ruleId (regime given by the rule of item 4: 0 unclassified, 1 Mode 1, 2 Stable,
%      3 Mode 2, 4 Mode 3).
%   4. The decision rule of the regime flow chart (resolved peak -> coherence -> phase -> frequency) is applied
%      to every step and compared with regimeId: console table and tables\Regime_rule_check.csv.
%
% CHANGES of v4.1 (kept)
%   - Fig_01 / Fig_02 (flow PSD maps) come from the RAW 200 Hz flow record, frequency axis up to the Nyquist
%     frequency (fplotQ): the flow stays a low-frequency signal although it was sampled at 200 Hz.
%   - Every figure window follows its size (scaleWindow); PNG / PDF are written at the paper size first.
%   - Output folder is Fig_4_1_2_v4_5, so the earlier figures are not overwritten.
%
% MAIN FIGURES (one figure = one axes; assemble the panels in Word)
%   Fig_01  Ch1 flow PSD map (200 Hz)        Fig_03  Ch1 temperature PSD map
%   Fig_02  Ch2 flow PSD map (200 Hz)        Fig_04  Ch2 temperature PSD map
%   Fig_05  correlation matrix  (rows = variables, columns = the 19 power steps)
%   Fig_06  phase-lag matrix    (same layout; phase of the variable behind the reference flow)
%   Fig_07  mean dP-q characteristic (optional, makeFigDPq)          Fig_S1  coherence matrix (optional)
%   FigS2   spectra of one step per regime, flow (row 1) and temperature (row 2) (optional, makeSpecFig)
%   FigS3   channel coupling of one step per regime: coherence, phase, common / differential mode (makeCohSpecFig)
%   FigS4   band metrics against the 19 power steps: A (row 1) and rho (row 2) of the bands, the 2-100 Hz band last (makeBandFig)
% TABLES (outDir\tables, plain CSV; row labels in column 1, step labels in row 1)
%   Corr_Pearson, Corr_Spearman, Corr_Pearson_p, Corr_Spearman_p, Corr_Pearson_sig, Corr_Spearman_sig,
%   Corr_Neff, Coherence_fp, Phase_deg, PhaseSE_deg, Lag_s, Fref_Hz, Step_table, Regime_rule_check, Band_table
%   (Step_table, Regime_rule_check, Band_table: one row per power step, column names in row 1)
% FILE  outDir\Flow10Hz_windows.mat : q10{step,channel} and dP10{step,channel} (single, 10 Hz, 900 s windows)
%
% ROWS (20 = 2 channels x 10 variables, in this order per channel):
%   T_in, T_out, T_w1 ... T_w5, P_in, P_out, dP = P_in - P_out        (+ optional row 21: q2 vs q1)
% REFERENCE FLOW (refMode): 'own' = each channel's variables against its own flow (Ch1 vs q1, Ch2 vs q2)
%
% DEFINITIONS (all inside the same 900 s windows as Sec. 4.1.2)
%   normalisation (v4.5, normMode = 'mean'): x* = normScale*x/<x> in every window, <x> = mean|q| (flow) or the mean in kelvin (temperature).
%                 All PSDs below (flow of the 10 Hz copy and of the 200 Hz record, temperature) and everything derived from them are
%                 properties of x* (normScale = 100: S in %^2/Hz, V in %^2).  A_c,b = (100/normScale) sqrt(V_c,b) is the number of v4.4.
%                 normMode = 'none': S in (mL/s)^2/Hz or K^2/Hz and A_c,b = 100 sqrt(V_c,b)/mean|q_c| (v4.4).
%   correlation : zero-lag Pearson r (or Spearman rho) between the reference flow and the variable, both
%                 band-passed to corrBand (zero phase).  p-value from the effective sample size
%                 N_eff = N/(1 + 2*sum_k (1-k/N) rho_xx(k) rho_yy(k)) (Bartlett), k = 1..acfLag_s*fs;
%                 significant cells: Benjamini-Hochberg false-discovery rate 5 % over ALL cells (useFDR).
%                 Zero-lag r mixes coupling strength and phase: for a narrow-band relation r ~ gamma*cos(phi).
%   f_p         : frequency of the highest PSD value inside fband (0.01-2 Hz) of the UNSMOOTHED Welch spectrum of the
%                 flow, per channel (f_pc: of the sum of both channels).  It is the position of the spectral MAXIMUM, not
%                 of the largest share of the variance: the share within +-df_pk of f_p is share1 / share2 and the
%                 share below fLow_Hz is lfShare1 / lfShare2 in Step_table.csv.
%   phase       : angle of the cross-spectrum variable/flow at the dominant frequency f_p of the reference
%                 flow (Welch, 150 s Hann segments, 50 % overlap, nfft 4096); positive = variable LAGS flow.
%                 Cells with gamma^2 < gamma^2_crit (alpha = 0.05, K segments) are masked (grey).
%                 In the 'unresolved' steps f_p is the lower band edge (about 0.012 Hz): the phase there is
%                 that of the slow drift, not of an oscillation.
%   inter-ch.   : row 21 = q2 against q1 at the common dominant frequency f_pc; positive = Ch2 lags Ch1.
%   band metrics (Band_table.csv; same Welch spectra as above, df = nfft-bin width; bands are half-open [lo,hi) and
%                 the last one includes its upper edge, so the three bands add up to the 0.01-2 Hz total):
%     V_c,b     = sum of the one-sided flow PSD S_cc over the bins of band b times df  (variance of channel c in band b)
%     A_c,b     = 100 sqrt(V_c,b) / mean|q_c|      relative RMS of the flow fluctuation in the band, percent (Atot: whole 0.01-2 Hz);
%                 with the normalised PSD (v4.5) V_c,b is already relative and A_c,b = (100/normScale) sqrt(V_c,b): the same number
%     rho_b     = sum_b Re(S12) df / sqrt(V_1,b V_2,b),  Re(S12) = (S_sum - S_diff)/4, S_sum / S_diff = PSD of q1+q2 / q1-q2:
%                 correlation coefficient of the two band-limited flows (+1 in phase, -1 anti-phase, 0 unrelated)
%     gam_b     = |sum_b S12 df| / sqrt(V_1,b V_2,b)        ph_b = angle(sum_b S12 df)  (+ : Ch1 leads Ch2)
%     fpk_b     = frequency of the maximum of S11+S22 inside band b;  edge_b = 1 when that maximum is the first or the last
%                 bin of the band (then it is the slope of a larger feature outside the band, not a peak);  lev_b = its
%                 level, 10 log10((S11+S22)/2);  coh_b / phpk_b / sd_b = coherence, phase (deg) and 10 log10(S_sum/S_diff)
%                 between the channels at fpk_b;  prom_b = 10 log10 of the maximum over the median of S11+S22 in the band
%   high band (v4.4; last columns of Band_table.csv): the same V, A, rho, gam, ph for the bins f > fband(2) ... fHigh(2) (2-100 Hz)
%                 of the Welch spectra of the RAW 200 Hz record (same segment length in seconds, same bin width as above, so the
%                 four bands close the variance budget of 0.01-100 Hz; A is normalised with the same mean|q| as the other bands).
%     hfShare2  = 100 V_c(f > fband(2)) / V_c(0.01-100 Hz)     hfShare10 = the same for f > fHigh10 (10 Hz)       (percent)
%     psdDiff_c = max over the bins of fband of |10 log10( S_cc(10 Hz copy) / S_cc(200 Hz record) )|   (dB)
%   regime rule R1 (v4.4): see item A of the changes; f_p is the one defined above (maximum of the PSD in fband, per channel).
% NOTE: thermocouples cannot follow oscillations of about 0.8 Hz, so masked temperature cells in the
%       0.8 Hz regime are a physical result, not a bug.
%
% Needs: Signal Processing Toolbox (pwelch, mscohere, cpsd, butter, filtfilt, decimate, resample,
%        findpeaks), MATLAB R2020a+ (exportgraphics, smoothdata).  No Statistics Toolbox needed.
clc; clear; close all;

%% ===== 0. USER SETTINGS ======================================================
pressFile = 'I:\投稿2\数据\dat数据\20260309双通道\data20260309.mat';
flowFile  = 'I:\投稿3\sci\图片\4.1\Flow_Temperature_Data.mat';
outDir    = fullfile(fileparts(flowFile),'Fig_4_1_2_v4_5');   % figures
tabDir    = fullfile(outDir,'tables');                          % CSV tables
if ~exist(outDir,'dir'), mkdir(outDir); end
if ~exist(tabDir,'dir'), mkdir(tabDir); end

% ---- what is computed / drawn ------------------------------------------------
refMode       = 'own';     % reference flow of the matrices: 'own' | 'ch1' | 'ch2'
corrType      = 'pearson'; % statistic drawn in Fig_05: 'pearson' | 'spearman' (both are written to CSV)
psdTempVar    = 'Tin';     % temperature of the PSD maps: 'Tin' | 'Tout' | 'Tw1' ... 'Tw5'
addInterChRow = true;      % extra last row: q2 against q1 (inter-channel flow coupling)
makeFigDPq    = true;      % Fig_07: mean dP-q characteristic
makeCohFig    = false;     % Fig_S1: coherence matrix (same layout as Fig_06)
makeSpecFig   = true;      % FigS2: flow and temperature spectra of one representative step per regime
specSteps     = [];        % steps shown in FigS2 / FigS3, e.g. [4 10 16 19]; [] = middle step of every regime
makeCohSpecFig = true;     % FigS3: coherence, phase and common / differential-mode spectra between the channels
makeBandFig   = true;      % FigS4: band metrics (A, rho) against the 19 power steps
bandFigHigh   = true;      % FigS4: add the 2-100 Hz band as the last panel (false = the three bands of v4.3 only)
exportWindows = true;     % Flow10Hz_windows.mat: 10 Hz flow and dP of every analysis window (about 3 MB)
scaleWindow   = true;     % open figure windows follow their size when enlarged (false = fixed layout as in v4)
cellText      = true;      % numbers inside the matrix cells
cmapPSD       = 'magma';   % 'magma' | 'viridis' | 'ylgnbu'
cmapCoh       = 'viridis'; % colour map of the optional coherence matrix

% ---- normalisation of the variables (v4.5; see item A of the changes) ----
% 'mean': every variable is divided by its own mean over the analysis window before the spectra are computed (x* = normScale*x/mean(x)):
%         PSD, band amplitudes and map levels are then relative to the mean.   'none': absolute units, as v4.4.
normMode   = 'mean';
normScale  = 100;          % x* in percent of the window mean (100; PSD unit %^2/Hz) or as a fraction (1; PSD unit 1/Hz)
tempOffset = 273.15;       % added to the stored temperature before its mean is taken (records in degC -> kelvin; 0 = mean in the stored unit)
presOffset = 101.33;       % kPa added to P_in, P_out before their mean is taken (records are gauge pressure -> absolute; 0 if absolute); Cv_table only

% ---- sampling ---------------------------------------------------------------
fs_Q = 200;               % flow sampling rate (Hz)
fs_T = 10;                % temperature sampling rate (Hz)
fs_P = [];                % pressure sampling rate (Hz); [] = infer from record length
fs   = 10;                % common analysis rate (Hz)

% [start end] of the analysed window of each power step, seconds (900 s each)
win_s = [ ...
     710  1610;  2510  3410;  4310  5210;  6230  7130;  8090  8990; ...
   10190 11090; 11990 12890; 13790 14690; 15590 16490; 17390 18290; ...
   19610 20510; 21410 22310; 23270 24170; 25070 25970; 26870 27770; ...
   28670 29570; 30470 31370; 32270 33170; 34070 34970];
nW = size(win_s,1);

% heating power per channel (W): the 19 measured values (same as Q_e in the thermodynamic analysis)
P_W = [27.4 68.3 117.6 185.9 236.8 310.8 382.7 442.7 564.9 706.7 849.3 977.8 1062.5 ...
       1164.5 1264.2 1435.5 1542.9 1632.9 1750.3]';
assert(numel(P_W) == nW,'P_W must have one entry per window.');

% regime of each window: 1 Mode 1, 2 Stable, 3 Mode 2, 4 Mode 3  (P1 = start-up transient, flagged below)
regimeId     = [1 1 1 1 1 1  2 2 2 2 2 2 2  3 3 3 3 3  4];
regName      = {'Mode 1','Stable','Mode 2','Mode 3'};
startupSteps = 1;          % steps drawn grey in the step bar (start-up transient)
assert(numel(regimeId) == nW,'regimeId must have one entry per window.');

% ---- spectral estimation (identical to Sec. 4.1.2 text) ---------------------
seg_s = 150;  ovl = 0.5;  nfft = 4096;
fband    = [0.01 2];       % band searched for the dominant frequency (Hz)
fplot    = [0.01 2];       % frequency range of the temperature PSD maps (Hz)
fplotQ   = [0.01 fs_Q/2];  % frequency range of the FLOW PSD maps (Hz): raw 200 Hz record, up to Nyquist
clFlow   = [];             % colour limits of the flow maps (dB); [] = automatic (5th-99.5th percentile of the 0.01-2 Hz band)
alphaSig = 0.05;           % significance level of the coherence
df_pk    = 0.02;           % half-width of the "peak share" window (Hz)
fLow_Hz  = 0.05;           % upper edge of the "slow band" (fband(1) ... fLow_Hz) of the dominance indicators (Hz)

% ---- regime rule R1 (see item A of the changes): thresholds ----
ruleRhoSlow = 0;           % (1) rho of the slow band below this value (channels in anti-phase)          -> Mode 3
ruleASlow   = 5;           % (2) A_slow (percent of mean |q|) of the LARGER channel at or above this value -> Mode 1
ruleFpFast  = 0.4;         % (3) f_p (Hz) at or above this value in BOTH channels                         -> Mode 2; otherwise Stable

% ---- three frequency bands of the band metrics (Hz): slow | mid | fast; first edge = lower end of fband ----
% slow: pulsation of the low-power regime (about 0.015-0.035 Hz);  mid: slow oscillations and their harmonics (0.07-0.4 Hz);
% fast: the 0.4-2 Hz range around the oscillation near 0.8 Hz.  The edges are settings; the number of names must match.
bandEdges = [0.01 0.07 0.4 2];
bandName  = {'slow','mid','fast'};
% ---- fourth band "high" (v4.4), computed from the RAW 200 Hz flow: above the last band edge up to the Nyquist frequency ----
fHigh   = [bandEdges(end) fs_Q/2];   % [2 100] Hz
fHigh10 = 10;                        % the share of the variance above this frequency is reported as well (Hz)

% ---- correlation statistics -------------------------------------------------
corrBand  = fband;         % zero-phase band-pass before the correlation ([] = only linear detrend)
acfLag_s  = 90;            % lags (s) used for the effective sample size
alphaCorr = 0.05;          % significance level of the correlation tests
useFDR    = true;          % Benjamini-Hochberg control over all cells (false = uncorrected p < alphaCorr)

% ---- Type-I pulse detection (steps 1..nTypeI; gives the pulse frequencies quoted in Sec. 4.1.2) ----
nTypeI   = 5;
promFrac = 0.25;           % min. prominence = promFrac*(P95-P5) of the flow signal
minGap_s = 20;             % min. spacing between pulses (s)

% ---- map display (display only; the measured steps are marked by triangles) ----
% smoothBins = 1: no smoothing along frequency (v4.1 used 9, which flattened narrow spectral lines).
% nFq: display cells along the (log) frequency axis; 2400 cells over 4 decades give a cell about as wide as one FFT bin
% (0.00244 Hz) at 0.8 Hz, so a spectral line (width about 4 bins) is drawn with its true height.
nFq = 2400;  nPf = 300;  smoothBins = 1;

% ---- figure geometry (cm) ---------------------------------------------------
mapW = 9.5;   mapH = 7.2;   % PSD maps (two side by side = 19 cm)
matW = 19.0;  matH = 12.6; % matrices (full double-column width)
cellFont = 6;              % font size of the numbers inside the cells (pt)

regCol = [0.90 0.62 0.00; 0.34 0.71 0.91; 0.00 0.62 0.45; 0.80 0.47 0.65];   % Mode 1 / Stable / Mode 2 / Mode 3
cCh    = [0.00 0.45 0.70; 0.84 0.37 0.00];                                    % Ch1 blue, Ch2 vermilion

set(groot,'defaultAxesFontName','Arial','defaultTextFontName','Arial', ...
    'defaultAxesFontSize',8,'defaultAxesLineWidth',0.6,'defaultLineLineWidth',1.0, ...
    'defaultAxesTickDir','in','defaultAxesTickLength',[0.012 0.012],'defaultAxesBox','on');

% ---- variables (rows) ---------------------------------------------------------
vName = {'Tin','Tout','Tw1','Tw2','Tw3','Tw4','Tw5','Pin','Pout','dP'};
vTex  = {'T_{in}','T_{out}','T_{w1}','T_{w2}','T_{w3}','T_{w4}','T_{w5}','P_{in}','P_{out}','\Delta P'};
colT  = [6 7 1 2 3 4 5];   % columns of temCh1/temCh2: 6 inlet fluid, 7 outlet fluid, 1-5 wall sensors
nV    = numel(vName);
iTv   = find(strcmpi(vName,psdTempVar),1);
if isempty(iTv) || iTv > 7, error('psdTempVar must be ''Tin'', ''Tout'' or ''Tw1'' ... ''Tw5''.'); end

% ---- normalisation: switch and unit of the PSD (v4.5) ----
assert(any(strcmpi(normMode,{'mean','none'})),'normMode must be ''mean'' or ''none''.');
assert(any(normScale == [1 100]),'normScale must be 100 (percent) or 1 (fraction).');
doNorm = strcmpi(normMode,'mean');
if doNorm
    if normScale == 100, psdUnitQ = 'dB re 1 %^2/Hz'; else, psdUnitQ = 'dB re 1 Hz^{-1}'; end
    psdUnitT = psdUnitQ;                                   % flow and temperature are both relative now
else
    psdUnitQ = 'dB re (mL/s)^2/Hz';  psdUnitT = 'dB re K^2/Hz';
end

%% ===== 1. LOAD DATA ==========================================================
S    = load(flowFile,'Flo1','Flo2','temCh1','temCh2');
Q    = {double(S.Flo1(:)), double(S.Flo2(:))};                 % flow, mL/s (200 Hz)
Q    = cellfun(@(x) fillmissing(x,'linear'),Q,'UniformOutput',false);

D   = load(pressFile,'data20260309');
PFt = D.data20260309{1,2};
if istable(PFt), PF = PFt{:,:}; else, PF = PFt; end
PF   = double(PF);
Pin  = [PF(:,1), PF(:,5)];   Pout = [PF(:,2), PF(:,6)];

% ---- sampling-rate check (all records are assumed to start at t = 0) ----
dur_s = numel(Q{1})/fs_Q;
if isempty(fs_P)
    fs_P = size(PF,1)/dur_s;
    if abs(fs_P-round(fs_P))/fs_P < 0.01, fs_P = round(fs_P); end
end
fprintf('Record length from flow: %.1f s\n',dur_s);
fprintf('Temperature: %d samples -> %.3f Hz (assumed %g Hz)\n',size(S.temCh1,1),size(S.temCh1,1)/dur_s,fs_T);
fprintf('Pressure   : %d samples -> %.3f Hz (used %g Hz)\n',size(PF,1),size(PF,1)/dur_s,fs_P);
if abs(size(S.temCh1,1)/dur_s - fs_T)/fs_T > 0.02
    warning('Temperature record length does not match fs_T = %g Hz. Check fs_T / time alignment.',fs_T);
end
if abs(size(PF,1)/dur_s - fs_P)/fs_P > 0.02
    warning('Pressure record length does not match fs_P = %g Hz. Set fs_P manually.',fs_P);
end

% ---- common 10 Hz grid:  Qd{c} = flow,  Yd{c} = [T_in T_out T_w1..T_w5 P_in P_out dP] ----
Qd = cellfun(@(x) to10Hz(x,fs_Q,fs),Q,'UniformOutput',false);
Yd = cell(1,2);
for c = 1:2
    Tm = double(S.(sprintf('temCh%d',c)));
    if size(Tm,2) < 7
        error('temCh%d must have 7 columns: 5 wall sensors, inlet fluid, outlet fluid.',c);
    end
    TT = to10Hz(Tm(:,colT),fs_T,fs);
    PP = to10Hz([Pin(:,c), Pout(:,c), Pin(:,c)-Pout(:,c)],fs_P,fs);
    nn = min(size(TT,1),size(PP,1));
    Yd{c} = [TT(1:nn,:), PP(1:nn,:)];
end
need = max(win_s(:))*fs;
if min([numel(Qd{1}) numel(Qd{2}) size(Yd{1},1) size(Yd{2},1)]) < need
    error('At least one record is shorter than the last analysis window.');
end
% band-passed copies for the correlation (filtered once over the whole record, windows are cut afterwards)
Qb = cellfun(@(x) bandpass10(x,corrBand,fs),Qd,'UniformOutput',false);
Yb = cellfun(@(x) bandpass10(x,corrBand,fs),Yd,'UniformOutput',false);
clear S D PFt PF Pin Pout Tm TT PP

%% ===== 2. ROW DEFINITIONS AND SIGN CALIBRATION ================================
nRow   = 2*nV + double(addInterChRow);
rowCh  = [ones(1,nV) 2*ones(1,nV)];
rowVar = [1:nV 1:nV];
switch lower(refMode)
    case 'own', rowRef = rowCh;
    case 'ch1', rowRef = ones(1,2*nV);
    case 'ch2', rowRef = 2*ones(1,2*nV);
    otherwise,  error('refMode must be ''own'', ''ch1'' or ''ch2''.');
end
rowTex = [vTex, vTex];       rowCsv = [strcat('Ch1_',vName), strcat('Ch2_',vName)];
if addInterChRow
    rowTex = [rowTex, {'q_{2} vs q_{1}'}];   rowCsv = [rowCsv, {'q2_vs_q1'}];
end
grp     = [1 nV; nV+1 2*nV];   grpName = {'Ch1','Ch2'};
colCsv  = arrayfun(@(k) sprintf('P%d (%.1f W)',k,P_W(k)),1:nW,'UniformOutput',false);

% phase sign: phSign*angle(cpsd(x,y)) > 0  <=>  x leads y   (calibrated with a synthetic 90 deg lag)
tt = (0:1/fs:300)';  f0 = 0.5;
xa = sin(2*pi*f0*tt);  ya = sin(2*pi*f0*tt - pi/2);               % ya lags xa by 90 deg
[Pt,ft] = cpsd(xa,ya,hann(round(seg_s*fs/3)),[],nfft,fs);
[~,kt]  = min(abs(ft-f0));
phSign  = sign(angle(Pt(kt)));

%% ===== 3. WINDOW-BY-WINDOW ANALYSIS ==========================================
nseg = round(seg_s*fs);  nov = round(ovl*nseg);  w = hann(nseg);
f  = (0:nfft/2)'*fs/nfft;  nf = numel(f);
ib = f>=fband(1) & f<=fband(2);  fb = f(ib);
Kavg  = floor((round(diff(win_s(1,:))*fs) - nov)/(nseg - nov));
gcrit = 1 - alphaSig^(1/(Kavg-1));                               % coherence threshold
Macf  = round(acfLag_s*fs);

% ---- band partition of fband for the band metrics: bins of band b are  lo <= f < hi  (last band: f <= hi) ----
nB = numel(bandEdges) - 1;
assert(numel(bandName) == nB,'bandName must have one entry per band (numel(bandEdges)-1).');
assert(all(diff(bandEdges) > 0) && abs(bandEdges(1)-fband(1)) < 1e-9 && abs(bandEdges(end)-fband(2)) < 1e-9, ...
       'bandEdges must increase and run from fband(1) to fband(2).');
dfBin = f(2) - f(1);                                    % FFT bin width of the 10 Hz copy (Hz)
bIdx  = cell(1,nB);
for b = 1:nB
    if b < nB, bIdx{b} = find(f >= bandEdges(b)*(1-1e-9) & f <  bandEdges(b+1)*(1-1e-9));
    else,      bIdx{b} = find(f >= bandEdges(b)*(1-1e-9) & f <= bandEdges(b+1)*(1+1e-9));
    end
    assert(numel(bIdx{b}) >= 3,'Band %d (%g-%g Hz) contains fewer than 3 FFT bins.',b,bandEdges(b),bandEdges(b+1));
end
assert(isequal(vertcat(bIdx{:}),find(ib)),'The bands must cover fband exactly once.');
iSl = find(strcmp(bandName,'slow'),1);                  % band used by the first test of rule R1
assert(~isempty(iSl),'Rule R1 needs a band named ''slow'' in bandName.');

PSDq = {zeros(nf,nW), zeros(nf,nW)};      % flow PSD of the 10 Hz copy (dominant frequency, coherence, tables)
PSDt = {zeros(nf,nW), zeros(nf,nW)};      % temperature PSD maps (variable psdTempVar)
% flow PSD of the RAW 200 Hz record (Fig_01 / Fig_02): same segment length in seconds, same bin width as above
assert(fplotQ(2) <= fs_Q/2,'fplotQ(2) must not exceed the Nyquist frequency fs_Q/2.');
nseg200 = round(seg_s*fs_Q);  nov200 = round(ovl*nseg200);  w200 = hann(nseg200);
nfft200 = 2*round(nfft*fs_Q/fs/2);
f200    = (0:nfft200/2)'*fs_Q/nfft200;  nf200 = numel(f200);
PSDq200 = {zeros(nf200,nW), zeros(nf200,nW)};
% ---- high band (v4.4): bins of the 200 Hz spectra above the last band edge, and agreement with the 10 Hz grid ----
df200 = f200(2) - f200(1);                                   % bin width of the 200 Hz spectra (Hz)
ib200 = f200 >= fband(1) & f200 <= fband(2);                 % the bins of fband on the 200 Hz grid
assert(abs(df200-dfBin) < 1e-9*dfBin && nnz(ib200) == nnz(ib) && max(abs(f200(ib200)-f(ib))) < 1e-9, ...
       'The 200 Hz and the 10 Hz spectra must have the same frequency bins inside fband (check nfft, seg_s, fs_Q, fs).');
assert(fHigh(1) >= fband(2) && fHigh(2) <= fs_Q/2 && fHigh10 > fHigh(1) && fHigh10 < fHigh(2), ...
       'fHigh must run from fband(2) to at most fs_Q/2, and fHigh10 must lie inside it.');
iHi   = f200 >  fHigh(1)*(1+1e-9) & f200 <= fHigh(2)*(1+1e-9);   % high band:  fHigh(1) < f <= fHigh(2)
iHi10 = f200 >  fHigh10*(1+1e-9)  & f200 <= fHigh(2)*(1+1e-9);   % part of it above fHigh10
HA   = nan(nW,2);  Hrho = nan(nW,1);  Hgam = nan(nW,1);  Hph = nan(nW,1);   % amplitude of each channel / coupling of the channels
HSh2 = nan(nW,2);  HSh10 = nan(nW,2);                        % share (%) of the 0.01-100 Hz variance above fHigh(1) / above fHigh10
PsdDiff = nan(nW,2);                                         % largest |10 log10(S_10Hz copy / S_200Hz record)| inside fband (dB)
Fp = nan(nW,2);  Share = nan(nW,2);  qMean = nan(nW,2);  dpMean = nan(nW,2);
Iq = nan(nW,2);  Aq = nan(nW,2);  Rev = nan(nW,2);  PulseF = nan(nW,2);  PulseT = nan(nW,2);
PkdB = nan(nW,2);  LfdB = nan(nW,2);  PromdB = nan(nW,2);  LfShare = nan(nW,2);   % how dominant is the peak
Coh12fp = nan(nW,2);                                                               % gamma^2_12 at each channel's own f_p
fpc = nan(nW,1);  coh12 = nan(nW,1);  phi12 = nan(nW,1);  sdDB = nan(nW,1);
% ---- band metrics (see DEFINITIONS): one value per step and band ----
qRef  = nan(nW,2);                                                                 % mean |q| of the 10 Hz window (normalisation of A)
BA    = {nan(nW,nB), nan(nW,nB)};                                                  % A: relative RMS of the flow of channel 1 / 2 in the band
Atot  = nan(nW,2);                                                                 % A of the whole band fband
Brho  = nan(nW,nB);  Bgam = nan(nW,nB);  Bph = nan(nW,nB);                         % coupling of the two band-limited flows
Bfpk  = nan(nW,nB);  Bedge = nan(nW,nB); Blev = nan(nW,nB);                        % maximum of S11+S22 inside the band
Bcoh  = nan(nW,nB);  Bphpk = nan(nW,nB); Bsd = nan(nW,nB);  Bprom = nan(nW,nB);
tRef  = nan(nW,2);                      % mean of the PSD temperature of the window, kelvin (normalisation of the temperature)
aQ    = nan(nW,2);                      % A = aQ*sqrt(V): 100/normScale (normalised flow) or 100/mean|q| (absolute flow)
Cv    = nan(nW,2*(nV+1));               % std/mean (percent) of q and of the nV variables, channel 1 then channel 2 (Cv_table.csv)
PSDsum = zeros(nf,nW);  PSDdif = zeros(nf,nW);                                     % PSD of q1+q2 and q1-q2 (FigS3)
C12all = zeros(nf,nW);  PH12all = zeros(nf,nW);                                    % coherence / phase between the channels (FigS3)
q10 = cell(nW,2);  dP10 = cell(nW,2);                                              % 10 Hz windows (Flow10Hz_windows.mat)
Rp = nan(nRow,nW);  Rs = Rp;  pvP = Rp;  pvS = Rp;  Ne = Rp;
Coh = Rp;  LagDeg = Rp;  Fref = Rp;

if doNorm
    fprintf('Normalisation: normMode = ''mean'' -> every variable is divided by its own mean in each window (x* = %g*x/mean(x); flow: mean|q|; temperature: mean in K = stored + %g). PSD unit: %s\n', ...
            normScale,tempOffset,psdUnitQ);
else
    fprintf('Normalisation: normMode = ''none'' -> absolute units as in v4.4 (PSD unit: %s for the flow, %s for the temperature)\n',psdUnitQ,psdUnitT);
end
for k = 1:nW
    id = round(win_s(k,1)*fs)+1 : round(win_s(k,2)*fs);
    ir = round(win_s(k,1)*fs_Q)+1 : round(win_s(k,2)*fs_Q);
    xq = cell(1,2);  xr = cell(1,2);
    for c = 1:2
        % ---- normalisation (v4.5): x* = normScale*x/<x> over this window; flow: <x> = mean|q|, temperature: mean in kelvin ----
        qRef(k,c) = mean(abs(Qd{c}(id)));                              % mean |q| of the 10 Hz window
        tRef(k,c) = mean(Yd{c}(id,iTv)) + tempOffset;                  % mean of the PSD temperature of the window (K)
        assert(qRef(k,c) > 0 && tRef(k,c) > 0,'Step %d, channel %d: the mean of the flow or of the temperature is not positive (tempOffset?).',k,c);
        if doNorm
            nQ = normScale/qRef(k,c);   nT = normScale/tRef(k,c);   aQ(k,c) = 100/normScale;
        else
            nQ = 1;                     nT = 1;                     aQ(k,c) = 100/qRef(k,c);
        end
        xq{c} = nQ*detrend(Qd{c}(id),1);
        Sqq   = pwelch(xq{c},w,nov,nfft,fs);
        PSDq{c}(:,k) = Sqq;
        PSDt{c}(:,k) = pwelch(nT*detrend(Yd{c}(id,iTv),1),w,nov,nfft,fs);
        Sb = Sqq(ib);
        [~,j] = max(Sb);  fpk = fb(j);  Fp(k,c) = fpk;                 % f_p = maximum of the (unsmoothed) PSD in the band
        inn = abs(f-fpk) <= df_pk & ib;
        Share(k,c)  = trapz(f(inn),Sqq(inn))/trapz(fb,Sb);
        % how dominant is that maximum?  (same PSD Sb as for f_p)
        il = fb <= fLow_Hz;                                             % slow band  fband(1) ... fLow_Hz
        ob = abs(fb-fpk) >= 0.05 & abs(fb-fpk) <= 0.30;                 % neighbourhood of the maximum
        PkdB(k,c)    = 10*log10(Sb(j));                                 % level of the maximum
        LfdB(k,c)    = 10*log10(max(Sb(il)));                           % highest level in the slow band
        LfShare(k,c) = trapz(fb(il),Sb(il))/trapz(fb,Sb);               % share of the band variance in the slow band
        if any(ob), PromdB(k,c) = 10*log10(Sb(j)/median(Sb(ob))); end   % maximum against its neighbourhood
        qMean(k,c)  = mean(Qd{c}(id));
        dpMean(k,c) = mean(Yd{c}(id,nV));
        q10{k,c}    = single(Qd{c}(id));   dP10{k,c} = single(Yd{c}(id,nV));
        % intensity of the flow oscillation (raw flow, as in Sec. 4.1.2)
        qraw = Q{c}(ir);
        xr{c} = nQ*detrend(qraw,1);                                               % raw 200 Hz window, linear trend removed
        PSDq200{c}(:,k) = pwelch(xr{c},w200,nov200,nfft200,fs_Q);              % PSD of the raw 200 Hz record
        qref = mean(abs(qraw));
        Iq(k,c) = std(qraw)/qref;
        Aq(k,c) = (pct(qraw,95)-pct(qraw,5))/qref;
        Rev(k,c) = mean(qraw < 0);
        % standard deviation / mean of q and of every variable of this channel, percent (10 Hz copy, linear trend removed; Cv_table.csv)
        cvJ = (c-1)*(nV+1);
        Cv(k,cvJ+1) = 100*std(detrend(Qd{c}(id),1))/qRef(k,c);
        for cvV = 1:nV
            cvY = Yd{c}(id,cvV);
            if cvV <= 7, cvOff = tempOffset; elseif cvV <= 9, cvOff = presOffset; else, cvOff = 0; end
            Cv(k,cvJ+1+cvV) = 100*std(detrend(cvY,1))/abs(mean(cvY) + cvOff);
        end
        % Type-I pulse statistics
        if k <= nTypeI
            xpl = Qd{c}(id);
            prm = promFrac*(pct(xpl,95)-pct(xpl,5));
            [~,loc] = findpeaks(xpl,fs,'MinPeakProminence',prm,'MinPeakDistance',minGap_s);
            PulseF(k,c) = numel(loc)/(win_s(k,2)-win_s(k,1));
            if numel(loc) >= 3, PulseT(k,c) = mean(diff(loc)); end
        end
    end

    % ---- high band (fHigh) from the RAW 200 Hz record: amplitude of each channel and coupling of the channels ----
    % same Welch segmentation as the 10 Hz copy (150 s Hann, 50 %, bin width 0.00244 Hz); the cross-spectrum is one-sided like the PSD
    P12h = cpsd(xr{1},xr{2},w200,nov200,nfft200,fs_Q);
    S1h  = PSDq200{1}(:,k);  S2h = PSDq200{2}(:,k);
    Vh   = [sum(S1h(iHi))   sum(S2h(iHi))]*df200;             % variance of channel 1 / 2 above fHigh(1)
    Vh10 = [sum(S1h(iHi10)) sum(S2h(iHi10))]*df200;           % ... of which above fHigh10
    Vlo  = [sum(S1h(ib200)) sum(S2h(ib200))]*df200;           % variance in fband (from the 200 Hz spectra)
    X12h = sum(P12h(iHi))*df200;                              % band integral of the cross-spectrum
    HA(k,:)    = aQ(k,:).*sqrt(Vh);                           % percent of the mean |q| (same normalisation as the other bands)
    Hrho(k)    = real(X12h)/sqrt(Vh(1)*Vh(2));
    Hgam(k)    = abs(X12h)/sqrt(Vh(1)*Vh(2));
    Hph(k)     = phSign*angle(X12h)*180/pi;                   % + : Ch1 leads Ch2
    HSh2(k,:)  = 100*Vh./(Vlo + Vh);
    HSh10(k,:) = 100*Vh10./(Vlo + Vh);
    for c = 1:2                                               % does the 10 Hz copy have the spectrum of the 200 Hz record?
        PsdDiff(k,c) = max(abs(10*log10(max(PSDq{c}(ib,k),realmin)./max(PSDq200{c}(ib200,k),realmin))));
    end

    % ---- inter-channel flow coupling at the common dominant frequency f_pc ----
    S11 = PSDq{1}(:,k);  S22 = PSDq{2}(:,k);
    [~,j] = max(S11(ib)+S22(ib));  fpc(k) = fb(j);
    [~,jc12] = min(abs(f-fpc(k)));
    C12 = mscohere(xq{1},xq{2},w,nov,nfft,fs);
    P12 = cpsd(xq{1},xq{2},w,nov,nfft,fs);
    coh12(k) = C12(jc12);
    for c = 1:2, [~,jj] = min(abs(f-Fp(k,c)));  Coh12fp(k,c) = C12(jj); end       % coherence at each channel's own f_p
    phi12(k) = phSign*angle(P12(jc12))*180/pi;                    % + : Ch1 leads Ch2
    Psum  = pwelch(xq{1}+xq{2},w,nov,nfft,fs);
    Pdiff = pwelch(xq{1}-xq{2},w,nov,nfft,fs);
    sdDB(k) = 10*log10(Psum(jc12)/Pdiff(jc12));                   % > 0 in-phase, < 0 anti-phase
    PSDsum(:,k) = Psum;  PSDdif(:,k) = Pdiff;  C12all(:,k) = C12;  PH12all(:,k) = phSign*angle(P12)*180/pi;

    % ---- band metrics: amplitude of each channel, coupling of the channels, maximum of the common spectrum ----
    Sc = S11 + S22;                                               % common spectrum
    for b = 1:nB
        ii  = bIdx{b};
        V1  = sum(S11(ii))*dfBin;   V2 = sum(S22(ii))*dfBin;      % variance of channel 1 / 2 in the band
        Rb  = sum(Psum(ii) - Pdiff(ii))/4*dfBin;                  % band integral of Re(S12) = (S_sum - S_diff)/4
        X12 = sum(P12(ii))*dfBin;                                 % band integral of the cross-spectrum
        BA{1}(k,b) = aQ(k,1)*sqrt(V1);   BA{2}(k,b) = aQ(k,2)*sqrt(V2);   % percent
        Brho(k,b)  = Rb/sqrt(V1*V2);
        Bgam(k,b)  = abs(X12)/sqrt(V1*V2);
        Bph(k,b)   = phSign*angle(X12)*180/pi;                    % + : Ch1 leads Ch2
        [~,j] = max(Sc(ii));  jj = ii(j);
        Bfpk(k,b)  = f(jj);   Bedge(k,b) = double(j == 1 || j == numel(ii));
        Blev(k,b)  = 10*log10(Sc(jj)/2);                          % mean of the two channel PSDs at the maximum
        Bcoh(k,b)  = C12(jj);  Bphpk(k,b) = phSign*angle(P12(jj))*180/pi;
        Bsd(k,b)   = 10*log10(Psum(jj)/Pdiff(jj));
        Bprom(k,b) = 10*log10(Sc(jj)/median(Sc(ii)));             % maximum against the median of the band
    end
    Atot(k,:) = aQ(k,:).*sqrt([sum(S11(ib)) sum(S22(ib))]*dfBin);

    % ---- matrix rows: variable (or q2) against the reference flow ----
    for r = 1:nRow
        if r <= 2*nV
            c = rowCh(r);  v = rowVar(r);  rc = rowRef(r);
            ysp  = detrend(Yd{c}(id,v),1);   xsp = xq{rc};         % spectral analysis (linear detrend)
            ycr  = Yb{c}(id,v);              xcr = Qb{rc}(id);     % correlation (band-passed)
            fr   = Fp(k,rc);
        else
            ysp  = xq{2};  xsp = xq{1};  ycr = Qb{2}(id);  xcr = Qb{1}(id);  fr = fpc(k);
        end
        [~,jc] = min(abs(f-fr));
        Cxy = mscohere(ysp,xsp,w,nov,nfft,fs);
        Pyx = cpsd(ysp,xsp,w,nov,nfft,fs);
        Coh(r,k)    = Cxy(jc);
        LagDeg(r,k) = -phSign*angle(Pyx(jc))*180/pi;              % + : variable lags the reference flow
        Fref(r,k)   = f(jc);
        [Rp(r,k),pvP(r,k),Ne(r,k)] = corrNeff(xcr,ycr,Macf);
        [Rs(r,k),pvS(r,k)]         = corrNeff(rankAvg(xcr),rankAvg(ycr),Macf);
    end
end
if any(tRef(:) > 473)                                  % temperature records that are already in kelvin
    warning('The temperature records look like kelvin (mean + tempOffset > 473 K): set tempOffset = 0, otherwise the normalised temperature PSD and Cv are wrong.');
end
unres = Fp <= fb(1);                                   % peak at the lower band edge = unresolved
PkLf  = PkdB - LfdB;                                   % dB, >= 0 : the maximum of the band lies this far above the slow band
LagS  = LagDeg./(360*Fref);                            % time lag, s (+ : variable lags)

% ---- regime rule R1, applied to every step (the tests are applied in this order; the first one that fires decides) ----
%   (1) rho_slow < ruleRhoSlow                  -> Mode 3     (anti-phase channels in the slow band)
%   (2) A_slow (larger channel) >= ruleASlow    -> Mode 1     (slow pulsation)
%   (3) f_p >= ruleFpFast in BOTH channels      -> Mode 2     (oscillation near 0.8 Hz)
%   (4) otherwise                               -> Stable
ASlow  = max(BA{1}(:,iSl),BA{2}(:,iSl));               % A_slow of the larger channel (percent)
RhoSl  = Brho(:,iSl);                                  % correlation of the two channels in the slow band
ruleId   = zeros(nW,1);                                % 0 unclassified, 1 Mode 1, 2 Stable, 3 Mode 2, 4 Mode 3
ruleTest = zeros(nW,1);                                % number (1-4) of the test that decided
for k = 1:nW
    if     RhoSl(k) < ruleRhoSlow,                     ruleId(k) = 4;  ruleTest(k) = 1;
    elseif ASlow(k) >= ruleASlow,                      ruleId(k) = 1;  ruleTest(k) = 2;
    elseif all(Fp(k,:) >= ruleFpFast),                 ruleId(k) = 3;  ruleTest(k) = 3;
    else,                                              ruleId(k) = 2;  ruleTest(k) = 4;
    end
end
PhSE  = min(180,sqrt(max(1-Coh,0)./(2*max(Coh,eps)*(Kavg-1)))*180/pi);   % approx. 1-sigma phase error, deg

if useFDR
    sigP = bhFDR(pvP,alphaCorr);   sigS = bhFDR(pvS,alphaCorr);
else
    sigP = pvP < alphaCorr;        sigS = pvS < alphaCorr;
end

%% ===== 4. TABLES =============================================================
writeMatCsv(fullfile(tabDir,'Corr_Pearson.csv'),    Rp,        rowCsv,colCsv);
writeMatCsv(fullfile(tabDir,'Corr_Spearman.csv'),   Rs,        rowCsv,colCsv);
writeMatCsv(fullfile(tabDir,'Corr_Pearson_p.csv'),  pvP,       rowCsv,colCsv);
writeMatCsv(fullfile(tabDir,'Corr_Spearman_p.csv'), pvS,       rowCsv,colCsv);
writeMatCsv(fullfile(tabDir,'Corr_Pearson_sig.csv'),  double(sigP),rowCsv,colCsv);
writeMatCsv(fullfile(tabDir,'Corr_Spearman_sig.csv'), double(sigS),rowCsv,colCsv);
writeMatCsv(fullfile(tabDir,'Corr_Neff.csv'),       Ne,        rowCsv,colCsv);
writeMatCsv(fullfile(tabDir,'Coherence_fp.csv'),    Coh,       rowCsv,colCsv);
writeMatCsv(fullfile(tabDir,'Phase_deg.csv'),       LagDeg,    rowCsv,colCsv);
writeMatCsv(fullfile(tabDir,'PhaseSE_deg.csv'),     PhSE,      rowCsv,colCsv);
writeMatCsv(fullfile(tabDir,'Lag_s.csv'),           LagS,      rowCsv,colCsv);
writeMatCsv(fullfile(tabDir,'Fref_Hz.csv'),         Fref,      rowCsv,colCsv);
stepNames = {'P_W','regime','fp1_Hz','fp2_Hz','fpc_Hz','share1','share2','coh12','phi12_deg','SumDiff_dB', ...
             'Iq1','Iq2','Aq1','Aq2','rev1','rev2','pulseF1_Hz','pulseF2_Hz','pulseT1_s','pulseT2_s', ...
             'qMean1','qMean2','dPmean1','dPmean2','fp1_unresolved','fp2_unresolved', ...
             'pk1_dB','pk2_dB','lf1_dB','lf2_dB','pkMinusLf1_dB','pkMinusLf2_dB','prom1_dB','prom2_dB', ...
             'lfShare1','lfShare2','coh12_fp1','coh12_fp2','ruleId'};
stepM = [P_W, regimeId(:), Fp(:,1), Fp(:,2), fpc, Share(:,1), Share(:,2), coh12, phi12, sdDB, ...
         Iq(:,1), Iq(:,2), Aq(:,1), Aq(:,2), Rev(:,1), Rev(:,2), PulseF(:,1), PulseF(:,2), ...
         PulseT(:,1), PulseT(:,2), qMean(:,1), qMean(:,2), dpMean(:,1), dpMean(:,2), ...
         double(unres(:,1)), double(unres(:,2)), ...
         PkdB(:,1), PkdB(:,2), LfdB(:,1), LfdB(:,2), PkLf(:,1), PkLf(:,2), PromdB(:,1), PromdB(:,2), ...
         LfShare(:,1), LfShare(:,2), Coh12fp(:,1), Coh12fp(:,2), ruleId];
writeTableCsv(fullfile(tabDir,'Step_table.csv'),stepNames,stepM);

% ---- Band_table.csv : band metrics, one row per step (A in percent; fpk / edge / lev / coh / phpk / sd / prom: maximum of S11+S22 in the band) ----
bandNames = {'step','P_W','regime','A1_tot_pct','A2_tot_pct'};
bandM     = [(1:nW)', P_W, regimeId(:), Atot(:,1), Atot(:,2)];
for b = 1:nB
    nm = bandName{b};
    bandNames = [bandNames, {['A1_' nm '_pct'],['A2_' nm '_pct'],['rho_' nm],['gam_' nm],['ph_' nm '_deg'], ...
                             ['fpk_' nm '_Hz'],['edge_' nm],['lev_' nm '_dB'],['coh_' nm],['phpk_' nm '_deg'], ...
                             ['sd_' nm '_dB'],['prom_' nm '_dB']}]; %#ok<AGROW>
    bandM = [bandM, BA{1}(:,b), BA{2}(:,b), Brho(:,b), Bgam(:,b), Bph(:,b), Bfpk(:,b), Bedge(:,b), Blev(:,b), ...
             Bcoh(:,b), Bphpk(:,b), Bsd(:,b), Bprom(:,b)]; %#ok<AGROW>
end
% v4.4: high band (2-100 Hz, raw 200 Hz record), variance shares above 2 / 10 Hz and 200 Hz-vs-10 Hz check, at the END
bandNames = [bandNames, {'A1_high_pct','A2_high_pct','rho_high','gam_high','ph_high_deg', ...
                         'hfShare2_1_pct','hfShare2_2_pct','hfShare10_1_pct','hfShare10_2_pct','psdDiff1_dB','psdDiff2_dB'}];
bandM     = [bandM, HA(:,1), HA(:,2), Hrho, Hgam, Hph, HSh2(:,1), HSh2(:,2), HSh10(:,1), HSh10(:,2), ...
             PsdDiff(:,1), PsdDiff(:,2)];
writeTableCsv(fullfile(tabDir,'Band_table.csv'),bandNames,bandM);

% ---- Cv_table.csv (v4.5): standard deviation / mean (percent) of every variable in every window ----
cvNames = {'step','P_W','regime'};
for c = 1:2
    cvNames = [cvNames, strcat(sprintf('Ch%d_',c),[{'q'}, vName],'_cv_pct')]; %#ok<AGROW>
end
writeTableCsv(fullfile(tabDir,'Cv_table.csv'),cvNames,[(1:nW)', P_W, regimeId(:), Cv]);

% ---- 10 Hz flow and dP of every analysis window (for checking / redesigning the regime rule outside this script) ----
if exportWindows
    fs10 = fs;
    save(fullfile(outDir,'Flow10Hz_windows.mat'),'q10','dP10','win_s','P_W','regimeId','fs10','-v7');
end
save(fullfile(outDir,'FreqCorr_Results.mat'),'Rp','Rs','pvP','pvS','sigP','sigS','Ne','Coh','LagDeg','PhSE', ...
     'LagS','Fref','rowCsv','colCsv','P_W','regimeId','Fp','fpc','coh12','phi12','sdDB','PSDq','PSDt','f', ...
     'PSDq200','f200','win_s','gcrit','Kavg','refMode','corrBand','vName','unres','PkdB','LfdB','PromdB', ...
     'LfShare','Coh12fp','ruleId','ruleTest','ASlow','RhoSl','ruleRhoSlow','ruleASlow','ruleFpFast', ...
     'bandEdges','bandName','BA','Atot','Brho','Bgam','Bph','Bfpk','Bedge','Blev', ...
     'Bcoh','Bphpk','Bsd','Bprom','qRef','PSDsum','PSDdif','C12all','PH12all', ...
     'fHigh','fHigh10','HA','Hrho','Hgam','Hph','HSh2','HSh10','PsdDiff', ...
     'normMode','normScale','tempOffset','presOffset','psdUnitQ','psdUnitT','tRef','Cv');
fprintf('Coherence threshold (alpha=%.2f, K=%d): gamma^2_crit = %.3f\n',alphaSig,Kavg,gcrit);
fprintf('Correlation: %d of %d cells significant (Pearson), %d (Spearman); FDR control = %d\n', ...
        nnz(sigP),nnz(~isnan(Rp)),nnz(sigS),useFDR);
big = abs(Rp-Rs) > 0.2 & ~isnan(Rp);
fprintf('Cells with |r_Pearson - rho_Spearman| > 0.2: %d (check these for spikes / non-linearity)\n',nnz(big));
fprintf('Phase cells above the coherence threshold: %d of %d\n',nnz(Coh >= gcrit),nnz(~isnan(Coh)));

% ---- regime rule R1 against regimeId (typed above) ----
ruleAbbr = {'uncl.','Mode 1','Stable','Mode 2','Mode 3'};              % index = ruleId + 1
isStart  = ismember((1:nW)',startupSteps(:));
okRule   = (ruleId(:) == regimeId(:)) | isStart;
fprintf('\n==== Regime rule R1 applied to every step (v4.4) ====\n');
fprintf('tests in this order, the first one that fires decides:\n');
fprintf('  (1) rho_slow < %g                          -> Mode 3  (the channels move in anti-phase in the slow band)\n',ruleRhoSlow);
fprintf('  (2) A_slow >= %g %% (larger of the channels) -> Mode 1  (slow pulsation)\n',ruleASlow);
fprintf('  (3) f_p >= %g Hz in BOTH channels          -> Mode 2  (oscillation near 0.8 Hz)\n',ruleFpFast);
fprintf('  (4) otherwise                              -> Stable\n');
fprintf('rho_slow = correlation of the two channels in the %s band (%g-%g Hz);  A1, A2 = RMS of the flow of channel 1 / 2 in that band,\n', ...
        bandName{iSl},bandEdges(iSl),bandEdges(iSl+1));
fprintf('%% of mean |q|;  fp1, fp2 = frequency (Hz) of the PSD maximum in %g-%g Hz (n.r. = maximum at the lower band edge = no resolved peak);\n',fband(1),fband(2));
fprintf('test = number of the test that decided\n');
fprintf(' step  rho_slow    A1_slow   A2_slow      fp1     fp2   test   rule    typed   check\n');
for k = 1:nW
    if unres(k,1), s1 = '  n.r.'; else, s1 = sprintf('%6.4f',Fp(k,1)); end
    if unres(k,2), s2 = '  n.r.'; else, s2 = sprintf('%6.4f',Fp(k,2)); end
    if isStart(k), mk = 'start-up'; elseif okRule(k), mk = 'ok'; else, mk = '<< DIFFERENT'; end
    fprintf(' P%-3d   %6.2f     %7.2f   %7.2f   %s  %s     %d    %-6s  %-6s  %s\n', ...
            k,RhoSl(k),BA{1}(k,iSl),BA{2}(k,iSl),s1,s2,ruleTest(k),ruleAbbr{ruleId(k)+1},regName{regimeId(k)},mk);
end
if all(okRule)
    fprintf('R1 reproduces the typed regimes in all %d steps (start-up steps are not compared).\n',nW-nnz(isStart));
else
    fprintf('%d step(s) differ: the rule or regimeId has to be revised.\n',nnz(~okRule));
end

% separation of the classes: for every test, the nearest steps of the two classes on either side of its threshold
isC   = ~isStart;                                                       % classified steps (start-up steps excluded)
mM3 = isC & regimeId(:) == 4;  mM1 = isC & regimeId(:) == 1;  mSt = isC & regimeId(:) == 2;  mM2 = isC & regimeId(:) == 3;
fpMin = min(Fp,[],2);                                                   % smaller f_p of the two channels
fprintf('Separation of the classes (typed regimes, start-up steps excluded; a wide gap = the thresholds are not critical):\n');
ovTxt = {'', '   << THE CLASSES OVERLAP: no threshold separates them'};
if any(mM3) && any(isC & ~mM3)
    r3 = max(RhoSl(mM3));  r0 = min(RhoSl(isC & ~mM3));
    fprintf('  (1) rho_slow : Mode 3 steps  max %8.3f      | all other steps       min %8.3f      (threshold %g)%s\n', ...
            r3,r0,ruleRhoSlow,ovTxt{1+(r3 >= r0)});
end
if any(mM1) && any(mSt | mM2)
    a1 = min(ASlow(mM1));  a2 = max(ASlow(mSt | mM2));
    fprintf('  (2) A_slow   : Mode 1 steps  min %8.2f %%    | Stable, Mode 2 steps  max %8.2f %%    (threshold %g %%; ratio %.2f; Mode 3 is decided by test 1)%s\n', ...
            a1,a2,ruleASlow,a1/a2,ovTxt{1+(a1 <= a2)});
end
if any(mM2) && any(isC & ~mM2)
    f2 = min(fpMin(mM2));  f0 = max(fpMin(isC & ~mM2));
    fprintf('  (3) f_p      : Mode 2 steps  min %8.4f Hz   | all other steps      max %8.4f Hz   (threshold %g Hz; ratio %.1f; f_p = smaller of the two channels)%s\n', ...
            f2,f0,ruleFpFast,f2/f0,ovTxt{1+(f2 <= f0)});
end
kb = find(diff(ruleId(:)) ~= 0)';                                       % boundaries given by the rule: middle between two steps
kb = kb(~isStart(kb)' & ~isStart(kb+1)');
if ~isempty(kb)
    fprintf('Boundaries given by R1 (middle between the two steps, heating power per channel):\n');
    for kk = kb
        fprintf('  P%d|P%d  %s -> %s  at %.1f W\n',kk,kk+1,ruleAbbr{ruleId(kk)+1},ruleAbbr{ruleId(kk+1)+1},(P_W(kk)+P_W(kk+1))/2);
    end
end
chkNames = {'step','P_W','typedRegime','ruleId','match','ruleTest','rho_slow','A1_slow_pct','A2_slow_pct', ...
            'fp1_Hz','fp2_Hz','fp1_unresolved','fp2_unresolved', ...
            'fpc_Hz','coh12','SumDiff_dB','coh12_fp1','coh12_fp2','pkMinusLf1_dB','pkMinusLf2_dB','prom1_dB','prom2_dB', ...
            'lfShare1','lfShare2'};
chkM = [(1:nW)', P_W, regimeId(:), ruleId, double(okRule), ruleTest, RhoSl, BA{1}(:,iSl), BA{2}(:,iSl), ...
        Fp(:,1), Fp(:,2), double(unres(:,1)), double(unres(:,2)), ...
        fpc, coh12, sdDB, Coh12fp(:,1), Coh12fp(:,2), PkLf(:,1), PkLf(:,2), PromdB(:,1), PromdB(:,2), ...
        LfShare(:,1), LfShare(:,2)];
writeTableCsv(fullfile(tabDir,'Regime_rule_check.csv'),chkNames,chkM);

% ---- band metrics of every step (console copy of tables\Band_table.csv) ----
fprintf('\n==== Band metrics of the flow (tables\\Band_table.csv) ====\n');
fprintf('bands:');
for b = 1:nB, fprintf('   %s %g-%g Hz',bandName{b},bandEdges(b),bandEdges(b+1)); end
fprintf('\n');
fprintf('A1, A2 = RMS of the flow fluctuation of channel 1 / 2 in the band as %% of mean |q|   (Atot = whole %g-%g Hz)\n',fband(1),fband(2));
fprintf('rho    = correlation coefficient between the two channels in the band: +1 in phase, -1 anti-phase, 0 unrelated\n');
fprintf('fpk    = frequency (Hz) of the maximum of S11+S22 inside the band;  * = that maximum is the first or last bin of the band\n');
hd1 = sprintf(' step regime   Atot1 Atot2');   hd2 = sprintf('                    %%      %%');
for b = 1:nB
    hd1 = [hd1 sprintf(' |%-26s',[' ' bandName{b} ' band'])];                        %#ok<AGROW>
    hd2 = [hd2 sprintf(' |    A1     A2    rho    fpk ')];                             %#ok<AGROW>
end
fprintf('%s\n%s\n',hd1,hd2);
for k = 1:nW
    ln = sprintf(' P%-3d %-7s %5.1f %5.1f',k,regName{regimeId(k)},Atot(k,1),Atot(k,2));
    for b = 1:nB
        if Bedge(k,b) > 0, st = '*'; else, st = ' '; end
        ln = [ln sprintf(' |%6.2f %6.2f %6.2f %6.3f%s',BA{1}(k,b),BA{2}(k,b),Brho(k,b),Bfpk(k,b),st)]; %#ok<AGROW>
    end
    fprintf('%s\n',ln);
end

% ---- high band from the raw 200 Hz record, variance above 2 / 10 Hz, 200 Hz-vs-10 Hz agreement (v4.4; last columns of Band_table.csv) ----
fprintf('\n==== High band %g-%g Hz from the raw %g Hz record, and agreement of the 10 Hz copy (last columns of tables\\Band_table.csv) ====\n', ...
        fHigh(1),fHigh(2),fs_Q);
fprintf('A1h, A2h = RMS of the flow fluctuation of channel 1 / 2 in the high band as %% of mean |q|;  rho_h = correlation of the channels in that band\n');
fprintf('share    = part of the %g-%g Hz variance of the channel that lies above %g Hz / above %g Hz (%%)\n',fband(1),fHigh(2),fHigh(1),fHigh10);
fprintf('dPSD     = largest difference (dB) between the PSD of the 10 Hz copy and the PSD of the 200 Hz record inside %g-%g Hz\n',fband(1),fband(2));
fprintf(' step regime    A1h   A2h  rho_h |  share >%g Hz  |  share >%g Hz | dPSD (dB)\n',fHigh(1),fHigh10);
fprintf('                  %%     %%        |   Ch1    Ch2   |   Ch1    Ch2  |  Ch1   Ch2\n');
for k = 1:nW
    fprintf(' P%-3d %-7s %5.2f %5.2f %6.2f | %5.1f  %5.1f  | %5.1f  %5.1f  | %5.2f %5.2f\n',k,regName{regimeId(k)}, ...
            HA(k,1),HA(k,2),Hrho(k),HSh2(k,1),HSh2(k,2),HSh10(k,1),HSh10(k,2),PsdDiff(k,1),PsdDiff(k,2));
end
fprintf('largest difference between the 10 Hz copy and the 200 Hz record over all steps: %.2f dB (Ch1), %.2f dB (Ch2)\n', ...
        max(PsdDiff(:,1)),max(PsdDiff(:,2)));

% ---- standard deviation / mean of the variables (v4.5; tables\Cv_table.csv) ----
iCv = [1, 1+find(strcmp(vName,'Tin'),1), 1+find(strcmp(vName,'Tout'),1), 1+find(strcmp(vName,'Pin'),1), 1+find(strcmp(vName,'dP'),1)];
fprintf('\n==== Standard deviation / mean of the variables, percent (tables\\Cv_table.csv: all %d variables of both channels) ====\n',2*(nV+1));
fprintf('10 Hz copy of each window, linear trend removed; mean of the flow = mean|q|, of T in kelvin (stored + %g), of P_in / P_out as absolute pressure (stored + %g kPa), of dP its own mean\n',tempOffset,presOffset);
fprintf('                       Ch1                                        Ch2\n');
fprintf(' step regime        q     Tin    Tout     Pin      dP          q     Tin    Tout     Pin      dP\n');
for k = 1:nW
    fprintf(' P%-3d %-7s',k,regName{regimeId(k)});
    for c = 1:2
        fprintf(' %8.3f%7.3f%7.3f%8.3f%8.3f ',Cv(k,(c-1)*(nV+1)+iCv));
    end
    fprintf('\n');
end

%% ===== 5. FIGURES 1-4 : PSD power-frequency maps =============================
regOrder = unique(regimeId(:)','stable');                         % regimes in order of appearance
iChg     = find(diff(regimeId(:)) ~= 0);
Pbound   = (P_W(iChg) + P_W(iChg+1))/2;                           % regime boundaries (mid-point of two steps)
Pedge    = [P_W(1); Pbound(:); P_W(end)];
regC     = zeros(1,numel(regOrder));
for i = 1:numel(regOrder), regC(i) = mean(Pedge(i:i+1)); end
if any(strcmpi(cmapPSD,{'magma','viridis'})), ovlCol = [1 1 1]; else, ovlCol = [0 0 0]; end

Zq = cell(1,2);  Zt = cell(1,2);
for c = 1:2
    [fqQ,Pf,Zq{c}] = buildMap(P_W,f200,PSDq200{c},fplotQ,nFq,nPf,smoothBins,true);    % flow: raw 200 Hz record, up to Nyquist
    [fq ,~ ,Zt{c}] = buildMap(P_W,f,   PSDt{c},   fplot, nFq,nPf,smoothBins,false);   % temperature: 10 Hz record
end
if isempty(clFlow)                                 % colour scale from the 0.01-2 Hz band only, so that the
    inb = fqQ <= fband(2)*(1+1e-9);                % featureless high-frequency part does not set the range
    zz  = [reshape(Zq{1}(inb,:),[],1); reshape(Zq{2}(inb,:),[],1)];
    clq = [pct(zz,5) pct(zz,99.5)];
else
    clq = clFlow;
end
zz = [Zt{1}(:); Zt{2}(:)];  clt = [pct(zz,5) pct(zz,99.5)];

om = struct('W',mapW,'H',mapH,'cmap',cmapPSD,'Pbound',Pbound,'regC',regC,'regNames',{regName(regOrder)}, ...
            'ovl',ovlCol);
for c = 1:2
    o = om;  o.cl = clq;  o.cbLab = ['PSD (' psdUnitQ ')'];
    o.tag = sprintf('Ch%d, q',c);
    o.fpY = Fp(:,c);  o.fpOK = ~unres(:,c);  o.fpStyle = 'peak';
    drawPSDMap(outDir,sprintf('Fig_%02d_Ch%d_flow_PSD_map',c,c),P_W,fqQ,Pf,Zq{c},o);
end
for c = 1:2
    o = om;  o.cl = clt;  o.cbLab = ['PSD (' psdUnitT ')'];
    o.tag = sprintf('Ch%d, %s',c,vTex{iTv});
    o.fpY = Fp(:,c);  o.fpOK = ~unres(:,c);  o.fpStyle = 'ref';
    drawPSDMap(outDir,sprintf('Fig_%02d_Ch%d_%s_PSD_map',2+c,c,vName{iTv}),P_W,fq,Pf,Zt{c},o);
end

%% ===== 6. FIGURES 5-6 (and S1) : matrices ====================================
colLab = arrayfun(@(p) sprintf('%d',round(p)),P_W,'UniformOutput',false);
mo = struct('W',matW,'H',matH,'rowTex',{rowTex},'grp',grp,'grpName',{grpName},'colLab',{colLab}, ...
            'regimeId',regimeId(:)','startup',startupSteps,'regName',{regName},'regCol',regCol, ...
            'cellText',cellText,'fontCell',cellFont,'xlab','Heating power (W)');

% ---- Fig_05 correlation matrix ----
if strcmpi(corrType,'spearman')
    Rdraw = Rs;  sigDraw = sigS;  tt5 = 'Spearman \rho';  nm5 = 'spearman';
else
    Rdraw = Rp;  sigDraw = sigP;  tt5 = 'Pearson r';       nm5 = 'pearson';
end
o = mo;  o.cmap = 'rdbu';  o.cl = [-1 1];  o.cbTicks = -1:0.5:1;  o.cbTickLab = {};
o.cbLab = [tt5 ' with the reference flow (zero lag)'];
o.fmt = @(v) regexprep(sprintf('%.2f',v),'^(-?)0\.','$1.');
o.strips = struct('label',{},'vals',{},'na',{});
drawMatrixFig(outDir,['Fig_05_correlation_matrix_' nm5],Rdraw,sigDraw & ~isnan(Rdraw), ...
              ~sigDraw & ~isnan(Rdraw),o);

% ---- Fig_06 phase-lag matrix (+ f_p strips) ----
if strcmpi(refMode,'own')
    stp = struct('label',{'f_{p}, Ch1','f_{p}, Ch2'},'vals',{Fp(:,1)',Fp(:,2)'},'na',{unres(:,1)',unres(:,2)'});
else
    rc = rowRef(1);
    stp = struct('label',{sprintf('f_{p}, Ch%d',rc)},'vals',{Fp(:,rc)'},'na',{unres(:,rc)'});
end
okPh = Coh >= gcrit;
o = mo;  o.cmap = 'twilight';  o.cl = [-180 180];  o.cbTicks = -180:90:180;
o.cbTickLab = {'-180','-90','0','90','180'};
o.cbLab = ['Phase lag behind the reference flow at f_p (' char(176) ')'];
o.fmt = @(v) sprintf('%d',round(v));
o.strips = stp;
drawMatrixFig(outDir,'Fig_06_phase_lag_matrix',LagDeg,okPh,false(size(okPh)),o);

% ---- Fig_S1 coherence matrix (optional) ----
if makeCohFig
    o = mo;  o.cmap = cmapCoh;  o.cl = [0 1];  o.cbTicks = 0:0.25:1;  o.cbTickLab = {};
    o.cbLab = ['Coherence \gamma^{2} at f_p (grey: below ' sprintf('%.2f',gcrit) ')'];
    o.fmt = @(v) regexprep(sprintf('%.2f',v),'^(-?)0\.','$1.');
    o.strips = stp;
    drawMatrixFig(outDir,'FigS1_coherence_matrix',Coh,okPh,~okPh & ~isnan(Coh),o);
end

%% ===== 7. FIGURE 7 : mean dP-q characteristic (optional) =====================
if makeFigDPq
    [fig,ax] = newAxesFig(9,7,[1.30 1.10 0.35 0.30]);
    mkReg = {'o','s','^','d'};
    hC = gobjects(1,2);
    for c = 1:2
        plot(ax,qMean(:,c),dpMean(:,c),'-','Color',cCh(c,:),'LineWidth',0.8);
        for r = 1:4
            ii = regimeId(:) == r;
            plot(ax,qMean(ii,c),dpMean(ii,c),mkReg{r},'Color',cCh(c,:), ...
                 'MarkerFaceColor',cCh(c,:),'MarkerSize',4.5,'LineStyle','none');
        end
        hC(c) = plot(ax,NaN,NaN,'-','Color',cCh(c,:));
        text(ax,qMean(1,c),  dpMean(1,c),  ' P1', 'FontSize',6.5,'Color',cCh(c,:));
        text(ax,qMean(nW,c), dpMean(nW,c), ' P19','FontSize',6.5,'Color',cCh(c,:));
    end
    hR = gobjects(1,4);
    for r = 1:4
        hR(r) = plot(ax,NaN,NaN,mkReg{r},'Color',[0.35 0.35 0.35],'MarkerFaceColor',[0.35 0.35 0.35],'MarkerSize',4.5);
    end
    xlabel(ax,'Mean flow rate, q (mL/s)');  ylabel(ax,'Mean \Delta P (as stored, kPa)');
    legend(ax,[hC hR],[{'Ch1','Ch2'}, regName],'Box','off','FontSize',6.5,'NumColumns',2,'Location','best');
    saveFig(fig,outDir,'Fig_07_dP_q_characteristic',false);
end

%% ===== 7a. one representative step per regime (FigS2, FigS3): the middle step of the regime, start-up steps excluded =====
if isempty(specSteps)
    for i = 1:numel(regOrder)
        kk = find(regimeId(:)' == regOrder(i) & ~ismember(1:nW,startupSteps));
        if ~isempty(kk), specSteps(end+1) = kk(ceil(numel(kk)/2)); end %#ok<SAGROW>
    end
end

%% ===== 7b. FIGURE S2 : spectra of one representative step per regime ========
% where is the maximum of the spectrum?  Row 1: flow PSD of the raw 200 Hz record (0.01 ... 100 Hz); row 2: temperature
% PSD (10 Hz record).  Circles: f_p of the flow (maximum of the PSD in the grey band fband); dashed lines in row 2:
% the same f_p on the temperature spectrum.  Dotted line: ruleFpFast (0.4 Hz), the lower limit of f_p for Mode 2 in rule R1.
if makeSpecFig
    nS  = numel(specSteps);
    sp  = struct('W',19,'H',10.6,'tag','FreqCorr_v4_5','k',specSteps,'P',P_W(specSteps)', ...
                 'rid',regimeId(specSteps),'regName',{regName(regimeId(specSteps))},'regCol',regCol,'cCh',cCh, ...
                 'fband',fband,'fc',ruleFpFast,'tLab',vTex{iTv});
    sp.ZF = cell(2,nS);  sp.ZT = cell(2,nS);
    sp.fp = zeros(2,nS);  sp.ok = false(2,nS);  sp.pkF = nan(2,nS);
    for c = 1:2
        [sp.fqF,ZF] = cellSpectrum(f200,PSDq200{c}(:,specSteps),[fplotQ(1) fplotQ(2)],3000,true);
        [sp.fqT,ZT] = cellSpectrum(f,   PSDt{c}(:,specSteps),   [fplot(1) fs/2],       900,false);
        for j = 1:nS
            sp.ZF{c,j} = ZF(:,j);   sp.ZT{c,j} = ZT(:,j);
            sp.fp(c,j) = Fp(specSteps(j),c);
            sp.ok(c,j) = ~unres(specSteps(j),c);
            sp.pkF(c,j) = 10*log10(PSDq200{c}(round(Fp(specSteps(j),c)/df200)+1,specSteps(j)));   % level of the maximum
        end
    end
    drawSpecFig(outDir,'FigS2_regime_spectra',sp);
end

%% ===== 7c. FIGURE S3 : coupling between the two channels, one representative step per regime ===========
% row 1 coherence, row 2 phase, row 3 spectra of the common mode (q1+q2)/2 and the differential mode (q1-q2)/2
if makeCohSpecFig
    sp3 = struct('W',19,'H',14.8,'tag','FreqCorr_v4_5','k',specSteps,'P',P_W(specSteps)', ...
                 'rid',regimeId(specSteps),'regName',{regName(regimeId(specSteps))},'regCol',regCol, ...
                 'f',fb,'fband',fband,'gcrit',gcrit,'bandEdges',bandEdges,'bandName',{bandName});
    sp3.C    = C12all(ib,specSteps);                                     % coherence between the channels
    sp3.PH   = PH12all(ib,specSteps);                                    % phase (deg), + : Ch1 leads Ch2
    sp3.Ssum = 10*log10(max(PSDsum(ib,specSteps)/4,1e-30));              % PSD of the common mode (q1+q2)/2, dB
    sp3.Sdif = 10*log10(max(PSDdif(ib,specSteps)/4,1e-30));              % PSD of the differential mode (q1-q2)/2, dB
    if doNorm                                                            % legend of the spectra of row 3
        sp3.legC = 'common (q^{*}_{1}+q^{*}_{2})/2';  sp3.legD = 'differential (q^{*}_{1}-q^{*}_{2})/2';
    else
        sp3.legC = 'common (q_1+q_2)/2';              sp3.legD = 'differential (q_1-q_2)/2';
    end
    drawCohFig(outDir,'FigS3_channel_coupling',sp3);
end

%% ===== 7d. FIGURE S4 : band metrics against the 19 power steps ===================
if makeBandFig
    bf = struct('W',19,'H',9.6,'tag','FreqCorr_v4_5','P',P_W(:)','rid',regimeId(:)','startup',startupSteps, ...
                'A1',BA{1},'A2',BA{2},'rho',Brho,'bandEdges',bandEdges,'bandName',{bandName}, ...
                'regName',{regName},'regCol',regCol,'cCh',cCh);
    if bandFigHigh                                                   % fourth panel: the 2-100 Hz band (200 Hz record)
        bf.A1 = [bf.A1, HA(:,1)];   bf.A2 = [bf.A2, HA(:,2)];   bf.rho = [bf.rho, Hrho];
        bf.bandEdges = [bandEdges, fHigh(2)];   bf.bandName = [bandName, {'high'}];
    end
    drawBandFig(outDir,'FigS4_band_metrics',bf);
end

%% ===== 8. OPEN FIGURE WINDOWS FOLLOW THEIR SIZE ==============================
% (all files are already written at the paper size; this only changes the windows that stay open on screen)
if scaleWindow
    hF = findall(0,'Type','figure','Tag','FreqCorr_v4_5');       % only the figures made by this script
    for i = 1:numel(hF), makeScalable(hF(i)); end
end

fprintf('All figures written to: %s\nAll tables written to : %s\n',outDir,tabDir);

%% ===== local functions =======================================================
function y = to10Hz(x,fsx,fsout)
% anti-aliased resampling of every column of x to the common analysis rate
    x  = fillmissing(double(x),'linear');
    nc = size(x,2);  r = fsx/fsout;
    for j = nc:-1:1
        xj = x(:,j);
        if abs(r-1) < 1e-6
            yj = xj;
        elseif abs(r-20) < 1e-6
            yj = decimate(decimate(xj,4),5);            % 200 -> 10 Hz, zero phase
        else
            [p,q] = rat(fsout/fsx,1e-6);
            yj = resample(xj,p,q);
        end
        if j == nc, y = zeros(numel(yj),nc); end
        y(:,j) = yj;
    end
end

function y = bandpass10(x,band,fs)
% zero-phase band-pass of every column: 2nd-order high-pass + 4th-order low-pass Butterworth (filtfilt)
    y = x;
    if isempty(band), return; end
    if band(1) > 0
        [b,a] = butter(2,band(1)/(fs/2),'high');  y = filtfilt(b,a,y);
    end
    if band(2) < fs/2
        [b,a] = butter(4,band(2)/(fs/2),'low');   y = filtfilt(b,a,y);
    end
end

function v = pct(x,p)
% percentile without the Statistics Toolbox
    x = sort(x(~isnan(x)));  n = numel(x);
    v = x(max(1,min(n,round(p/100*n))));
end

function a = acfFFT(x,M)
% normalised (biased) autocorrelation, lags 0..M
    N  = numel(x);  nf = 2^nextpow2(2*N-1);
    X  = fft(x,nf);
    c  = real(ifft(abs(X).^2));
    a  = c(1:M+1)/c(1);
end

function [r,p,neff] = corrNeff(x,y,M)
% zero-lag correlation coefficient with a two-sided p-value based on the effective sample size
%   N_eff = N/(1 + 2*sum_{k=1..M} (1-k/N) rho_xx(k) rho_yy(k))          (Bartlett)
%   t = r*sqrt((N_eff-2)/(1-r^2)),  p = betainc((N_eff-2)/(N_eff-2+t^2),(N_eff-2)/2,1/2)
    x = x(:) - mean(x);   y = y(:) - mean(y);   N = numel(x);
    sxx = x'*x;  syy = y'*y;
    if sxx <= 0 || syy <= 0
        r = NaN;  p = NaN;  neff = NaN;  return;
    end
    r  = (x'*y)/sqrt(sxx*syy);
    ax = acfFFT(x,M);   ay = acfFFT(y,M);
    k  = (1:M)';
    s  = sum((1 - k/N).*ax(2:end).*ay(2:end));
    neff = N/(1 + 2*s);
    neff = min(max(neff,3),N);
    nu = neff - 2;
    r2 = min(r^2,1-1e-12);
    t2 = r2*nu/(1-r2);
    p  = betainc(nu/(nu+t2),nu/2,0.5);
end

function r = rankAvg(x)
% ranks with ties sharing the mean rank (Statistics Toolbox not needed)
    x = x(:);  n = numel(x);
    [xs,ix] = sort(x);
    newG = [true; xs(2:end) ~= xs(1:end-1)];
    r = zeros(n,1);
    if all(newG)                                   % no ties (the usual case for filtered data)
        r(ix) = (1:n)';
    else
        gid = cumsum(newG);
        gmn = accumarray(gid,(1:n)')./accumarray(gid,1);
        r(ix) = gmn(gid);
    end
end

function sig = bhFDR(p,q)
% Benjamini-Hochberg false-discovery-rate control; NaN p-values are never significant
    sig = false(size(p));
    ok  = find(~isnan(p));
    [ps,ix] = sort(p(ok));
    m = numel(ps);
    kmax = find(ps(:) <= q*(1:m)'/m,1,'last');
    if ~isempty(kmax), sig(ok(ix(1:kmax))) = true; end
end

function cm = getCmap(name,n)
% perceptually uniform / ColorBrewer maps, interpolated from anchor colours
    if nargin < 2, n = 256; end
    switch lower(name)
        case 'magma'
            A = [
        0.0015 0.0005 0.0139;
        0.0396 0.0311 0.1335;
        0.1131 0.0655 0.2768;
        0.2117 0.0620 0.4186;
        0.3167 0.0717 0.4854;
        0.4147 0.1104 0.5047;
        0.5128 0.1482 0.5076;
        0.6136 0.1818 0.4985;
        0.7164 0.2150 0.4753;
        0.8169 0.2559 0.4365;
        0.9043 0.3196 0.3881;
        0.9609 0.4183 0.3596;
        0.9867 0.5356 0.3822;
        0.9961 0.6537 0.4462;
        0.9969 0.7696 0.5349;
        0.9924 0.8843 0.6401;
        0.9871 0.9914 0.7495;
                ];
        case 'viridis'
            A = [
        0.2670 0.0049 0.3294;
        0.2788 0.0621 0.3866;
        0.2832 0.1208 0.4406;
        0.2788 0.1755 0.4834;
        0.2680 0.2235 0.5120;
        0.2504 0.2743 0.5331;
        0.2297 0.3224 0.5457;
        0.2105 0.3637 0.5522;
        0.1906 0.4071 0.5561;
        0.1727 0.4488 0.5579;
        0.1577 0.4859 0.5580;
        0.1419 0.5265 0.5560;
        0.1276 0.5669 0.5506;
        0.1197 0.6038 0.5414;
        0.1263 0.6441 0.5253;
        0.1579 0.6838 0.5017;
        0.2080 0.7187 0.4729;
        0.2815 0.7552 0.4326;
        0.3692 0.7889 0.3829;
        0.4587 0.8164 0.3297;
        0.5655 0.8424 0.2629;
        0.6785 0.8637 0.1895;
        0.7833 0.8793 0.1254;
        0.8963 0.8936 0.0963;
        0.9932 0.9062 0.1439;
                ];
        case 'twilight'      % cyclic (twilight_shifted): light = 0 deg, dark = +-180 deg
            A = [
        0.1874 0.0771 0.2162;
        0.2374 0.0674 0.3011;
        0.3055 0.0940 0.4337;
        0.3506 0.1660 0.5615;
        0.3672 0.2688 0.6496;
        0.3724 0.3677 0.7006;
        0.3841 0.4614 0.7309;
        0.4180 0.5490 0.7488;
        0.4867 0.6342 0.7626;
        0.5830 0.7096 0.7793;
        0.7030 0.7782 0.8092;
        0.8190 0.8335 0.8550;
        0.8857 0.8500 0.8857;
        0.8656 0.8144 0.8077;
        0.8294 0.7355 0.6733;
        0.8003 0.6410 0.5373;
        0.7759 0.5355 0.4241;
        0.7445 0.4347 0.3558;
        0.6981 0.3383 0.3221;
        0.6358 0.2494 0.3121;
        0.5566 0.1727 0.3142;
        0.4538 0.1162 0.3097;
        0.3455 0.0857 0.2803;
        0.2477 0.0687 0.2373;
        0.1849 0.0794 0.2131;
                ];
        case 'rdbu'          % ColorBrewer RdBu, blue = negative, red = positive
            A = [
        0.0196 0.1882 0.3804;
        0.1294 0.4000 0.6745;
        0.2627 0.5765 0.7647;
        0.5725 0.7725 0.8706;
        0.8196 0.8980 0.9412;
        0.9686 0.9686 0.9686;
        0.9922 0.8588 0.7804;
        0.9569 0.6471 0.5098;
        0.8392 0.3765 0.3020;
        0.6980 0.0941 0.1686;
        0.4039 0.0000 0.1216;
                ];
        case 'ylgnbu'
            A = [255 255 217; 199 233 180; 65 182 196; 34 94 168; 8 29 88]/255;
        otherwise
            error('Unknown colour map: %s',name);
    end
    cm = interp1(linspace(0,1,size(A,1)),A,linspace(0,1,n));
    cm = min(max(cm,0),1);
end

function [fig,ax] = newAxesFig(W,H,mar)
% one figure = one axes with exact margins; mar = [left bottom right top] in cm
    fig = figure('Units','centimeters','Position',[2 2 W H],'Color','w', ...
                 'PaperUnits','centimeters','PaperSize',[W H],'Tag','FreqCorr_v4_5');
    ax  = axes(fig,'Units','centimeters','Position',[mar(1) mar(2) W-mar(1)-mar(3) H-mar(2)-mar(4)]);
    hold(ax,'on');  box(ax,'on');
end

function [fq,Zl] = cellSpectrum(f,S,fplot,nFq,avgCells)
% spectra (columns of S, linear power) on a log-frequency grid of nFq display cells, in dB
%   avgCells = false : grid from the first to the last FFT bin inside fplot, linear interpolation of the dB spectrum
%   avgCells = true  : grid from fplot(1) to fplot(2) exactly; a display cell that contains >= 2 FFT bins shows the
%                      MEAN POWER of those bins (needed when the range spans decades, e.g. 0.01 ... 100 Hz:
%                      point sampling would alias the noise floor and could miss narrow lines)
    if ~avgCells
        ip = f>=fplot(1) & f<=fplot(2);
        f1 = f(find(ip,1));  f2 = f(find(ip,1,'last'));    % first / last FFT bin inside the range
        fq = logspace(log10(f1),log10(f2),nFq)';
        fq(1) = f1;  fq(end) = f2;                         % avoid NaN from round-off at the ends
        Z  = 10*log10(max(S(ip,:),1e-30));
        Zl = interp1(f(ip),Z,fq,'linear');
    else
        fq = logspace(log10(fplot(1)),log10(fplot(2)),nFq)';
        fq(1) = fplot(1);  fq(end) = fplot(2);
        Zl = interp1(f,10*log10(max(S,1e-30)),fq,'linear');
        rr = (fplot(2)/fplot(1))^(1/(nFq-1));              % ratio of neighbouring grid frequencies
        df = f(2) - f(1);
        j1 = max(1,ceil(fq/sqrt(rr)/df - 1e-9) + 1);       % first FFT bin of every display cell
        j2 = min(numel(f),ceil(fq*sqrt(rr)/df - 1e-9));    % last FFT bin of every display cell
        nb = j2 - j1 + 1;
        cs = [zeros(1,size(S,2)); cumsum(S,1)];            % running sum -> mean power of bins j1 ... j2
        for i = find(nb(:)' >= 2)
            Zl(i,:) = 10*log10(max((cs(j2(i)+1,:) - cs(j1(i),:))/nb(i),1e-30));
        end
    end
end

function [fq,Pf,Zs] = buildMap(P,f,S,fplot,nFq,nPf,smoothBins,avgCells)
% display map: log-frequency grid (cellSpectrum), optional Gaussian smoothing along frequency (smoothBins > 1;
% 1 = none, so that narrow spectral lines keep their height), shape-preserving (pchip) interpolation between
% the measured power steps
    [fq,Zl] = cellSpectrum(f,S,fplot,nFq,avgCells);
    if smoothBins > 1, Zl = smoothdata(Zl,1,'gaussian',smoothBins); end
    Pf = linspace(P(1),P(end),nPf);
    Zs = interp1(P,Zl.',Pf,'pchip').';
end

function axC = drawCbar(fig,cmp,cl,ticks,tickLab,lab,posCm)
% colour bar drawn as a separate axes at an exact position (cm); no dependence on colorbar() layout rules
    axC = axes(fig,'Units','centimeters','Position',posCm);
    nc  = size(cmp,1);
    image(axC,[0 1],[cl(1) cl(2)],reshape(cmp,[nc 1 3]));
    set(axC,'YDir','normal','XTick',[],'YAxisLocation','right','TickDir','out','Box','on', ...
            'Layer','top','FontSize',6.5,'LineWidth',0.5,'TickLength',[0.12 0.12]/max(posCm(3:4)));
    xlim(axC,[0 1]);  ylim(axC,cl);
    if ~isempty(ticks), set(axC,'YTick',ticks); end
    if ~isempty(tickLab), set(axC,'YTickLabel',tickLab); end
    ylabel(axC,lab,'FontSize',7);
end

function drawPSDMap(outDir,name,P_W,fq,Pf,Z,o)
% power-frequency map: PSD (dB) over heating power and (log) frequency
    mar = [1.40 1.15 2.00 0.85];                       % left bottom right top, cm
    [fig,ax] = newAxesFig(o.W,o.H,mar);
    nW = numel(P_W);
    pcolor(ax,Pf,fq,Z);  shading(ax,'interp');
    colormap(ax,getCmap(o.cmap,256));
    set(ax,'CLim',o.cl,'YScale','log','Layer','top');
    view(ax,2);
    for b = 1:numel(o.Pbound)                          % regime boundaries
        plot3(ax,[o.Pbound(b) o.Pbound(b)],[fq(1) fq(end)],[1 1],'--','Color',o.ovl,'LineWidth',0.5);
    end
    ok = o.fpOK(:);
    if strcmp(o.fpStyle,'peak')                        % dominant frequency of the flow (resolved steps)
        plot3(ax,P_W(ok),o.fpY(ok),ones(nnz(ok),1),'o','MarkerSize',3.8, ...
              'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',0.5);
    else                                               % temperature map: flow f_p as a reference
        plot3(ax,P_W(ok),o.fpY(ok),ones(nnz(ok),1),'o','MarkerSize',3.8, ...
              'MarkerFaceColor','none','MarkerEdgeColor',o.ovl,'LineWidth',0.7);
    end
    rg = fq(end)/fq(1);                                % frequency range (ratio): positions below are fractions of the log axis
    plot3(ax,P_W,fq(1)*rg^0.019*ones(nW,1),ones(nW,1),'^','MarkerSize',3, ...
          'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',0.3);   % measured steps
    xl = [P_W(1) P_W(end)];
    xlim(ax,xl);  ylim(ax,[fq(1) fq(end)]);
    if fq(end) > 5                                     % wide range (flow map up to Nyquist): decade ticks + minor ticks
        yt = [0.001 0.01 0.1 1 10 100 1000];
        set(ax,'YMinorTick','on');
    else
        yt = [0.02 0.05 0.1 0.2 0.5 1 2];
    end
    yt  = yt(yt >= fq(1)*(1-1e-6) & yt <= fq(end)*(1+1e-6));
    ytl = arrayfun(@(v) sprintf('%g',v),yt,'UniformOutput',false);
    set(ax,'XTick',0:250:2000,'YTick',yt,'YTickLabel',ytl);
    xlabel(ax,'Heating power (W)');  ylabel(ax,'Frequency (Hz)');
    text(ax,xl(1)+0.025*diff(xl),fq(end)/rg^0.022,1,o.tag,'Color',o.ovl,'FontSize',8,'FontWeight','bold', ...
         'VerticalAlignment','top','HorizontalAlignment','left');
    for i = 1:numel(o.regC)                            % regime names above the panel
        xn = (o.regC(i)-xl(1))/diff(xl);
        text(ax,xn,1.012,o.regNames{i},'Units','normalized','HorizontalAlignment','center', ...
             'VerticalAlignment','bottom','FontSize',6.5,'FontWeight','bold','Color',[0.25 0.25 0.25], ...
             'Clipping','off');
    end
    axW = o.W-mar(1)-mar(3);  axH = o.H-mar(2)-mar(4);
    drawCbar(fig,getCmap(o.cmap,256),o.cl,[],{},o.cbLab,[mar(1)+axW+0.30 mar(2) 0.28 axH]);
    saveFig(fig,outDir,name,true);                     % true = raster PDF for maps
end

function drawSpecFig(outDir,name,sp)
% spectra of one representative step per regime (columns).  Row 1: flow PSD of the raw 200 Hz record; row 2:
% temperature PSD (10 Hz record); both channels.  Circles: f_p of the flow (maximum of the PSD inside the grey band);
% dashed lines in row 2: the same f_p on the temperature spectrum.  Dotted line: the f_p threshold of the Mode 2 test of R1.
    nCol = numel(sp.k);
    mar = [1.60 1.15 0.40 1.45];  gx = 0.45;  gy = 0.62;                  % left bottom right top, gaps (cm)
    pw  = (sp.W - mar(1) - mar(3) - (nCol-1)*gx)/nCol;
    ph  = (sp.H - mar(2) - mar(4) - gy)/2;
    fig = figure('Units','centimeters','Position',[2 2 sp.W sp.H],'Color','w','PaperUnits','centimeters', ...
                 'PaperSize',[sp.W sp.H],'Tag',sp.tag);
    yl = zeros(2,2);  ys = zeros(2,1);
    for r = 1:2                                                           % common y-limits of every row
        if r == 1, ZZ = sp.ZF; else, ZZ = sp.ZT; end
        v   = cell2mat(reshape(ZZ,[],1));
        top = max(v);  lo = pct(v,0.5);                                   % 0.5th percentile: ignore single deep notches
        yl(r,:) = [floor((min(lo,top-40)-2)/10)*10, ceil((top+1)/10)*10];
        ys(r)   = 20*(1 + (diff(yl(r,:)) > 90));                          % tick step, dB
    end
    xtAll = [0.01 0.1 1 10 100];
    for r = 1:2
        for j = 1:nCol
            ax = axes(fig,'Units','centimeters', ...
                      'Position',[mar(1)+(j-1)*(pw+gx), mar(2)+(2-r)*(ph+gy), pw, ph]);
            hold(ax,'on');  box(ax,'on');
            if r == 1, fq = sp.fqF;  ZZ = sp.ZF; else, fq = sp.fqT;  ZZ = sp.ZT; end
            patch(ax,[sp.fband(1) sp.fband(2) sp.fband(2) sp.fband(1)],yl(r,[1 1 2 2]),[0.93 0.93 0.93], ...
                  'EdgeColor','none');                                    % band searched for f_p
            h = gobjects(1,2);
            for c = 1:2
                h(c) = plot(ax,fq,ZZ{c,j},'-','Color',sp.cCh(c,:),'LineWidth',0.6);
            end
            hfc = plot(ax,sp.fc*[1 1],yl(r,:),':','Color',[0.30 0.30 0.30],'LineWidth',0.7);
            for c = 1:2
                if ~sp.ok(c,j), continue; end
                if r == 1
                    plot(ax,sp.fp(c,j),sp.pkF(c,j),'o','MarkerSize',4,'MarkerFaceColor',sp.cCh(c,:), ...
                         'MarkerEdgeColor','w','LineWidth',0.4);
                else
                    plot(ax,sp.fp(c,j)*[1 1],yl(r,:),'--','Color',sp.cCh(c,:),'LineWidth',0.7);
                end
            end
            set(ax,'XScale','log','Layer','top','FontSize',6.5,'LineWidth',0.5);
            xlim(ax,[fq(1) fq(end)]);  ylim(ax,yl(r,:));
            xt = xtAll(xtAll >= fq(1)*(1-1e-6) & xtAll <= fq(end)*(1+1e-6));
            if j > 1,     xt = xt(xt > fq(1)*(1+1e-6));  end                 % end ticks only at the outer panel edges
            if j < nCol,  xt = xt(xt < fq(end)*(1-1e-6)); end
            set(ax,'XTick',xt,'XTickLabel',arrayfun(@(v) sprintf('%g',v),xt,'UniformOutput',false));
            set(ax,'YTick',ceil(yl(r,1)/ys(r))*ys(r) : ys(r) : yl(r,2));
            if j > 1, set(ax,'YTickLabel',{}); end
            if r == 2, xlabel(ax,'Frequency (Hz)','FontSize',7); end
            if j == 1
                if r == 1, ylabel(ax,'Flow PSD (dB)','FontSize',7);
                else,      ylabel(ax,[sp.tLab ' PSD (dB)'],'FontSize',7);
                end
            end
            if r == 1
                text(ax,0.5,1+0.88/ph,sp.regName{j},'Units','normalized','HorizontalAlignment','center', ...
                     'VerticalAlignment','bottom','FontSize',7.5,'FontWeight','bold', ...
                     'Color',0.75*sp.regCol(sp.rid(j),:),'Clipping','off');
                text(ax,0.5,1+0.50/ph,sprintf('P%d, %.0f W',sp.k(j),sp.P(j)),'Units','normalized', ...
                     'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',6.5, ...
                     'Color',[0.25 0.25 0.25],'Clipping','off');
                haL = {'right','left'};
                for c = 1:2                                               % f_p of both channels, above the panel
                    if sp.ok(c,j), s = sprintf('Ch%d %.3g Hz',c,sp.fp(c,j));
                    else,          s = sprintf('Ch%d n.r.',c);
                    end
                    text(ax,0.5+(c-1.5)*0.04,1+0.10/ph,s,'Units','normalized', ...
                         'HorizontalAlignment',haL{c},'VerticalAlignment','bottom','FontSize',6, ...
                         'Color',sp.cCh(c,:),'Clipping','off');
                end
            end
            if j == 1
                if r == 1
                    legend(ax,[h hfc],{'Ch1','Ch2',sprintf('%g Hz',sp.fc)},'Box','off','FontSize',6,'Location','southwest');
                else
                    hdash = plot(ax,NaN,NaN,'--','Color',[0.30 0.30 0.30],'LineWidth',0.7);
                    legend(ax,[h hdash],{'Ch1','Ch2','flow f_p'},'Box','off','FontSize',6,'Location','southwest');
                end
            end
        end
    end
    saveFig(fig,outDir,name,false);
end

function drawCohFig(outDir,name,sp)
% coupling between the two channels for one representative step per regime (columns); frequency axis = fband
%   row 1 : coherence gamma^2 between q1 and q2 (dashed: critical value, alpha = 0.05)
%   row 2 : phase of the cross-spectrum (+ : Ch1 leads Ch2); coloured dots where gamma^2 >= critical, grey dots below it
%   row 3 : PSD of the common mode (q1+q2)/2 (black) and of the differential mode (q1-q2)/2 (grey): common mode above
%           differential mode = the channels move in phase, below = anti-phase
% dotted vertical lines: edges of the bands of the band metrics
    nCol = numel(sp.k);  nRow = 3;
    mar = [1.60 1.15 0.40 1.20];  gx = 0.45;  gy = 0.55;                  % left bottom right top, gaps (cm)
    pw  = (sp.W - mar(1) - mar(3) - (nCol-1)*gx)/nCol;
    ph  = (sp.H - mar(2) - mar(4) - (nRow-1)*gy)/nRow;
    fig = figure('Units','centimeters','Position',[2 2 sp.W sp.H],'Color','w','PaperUnits','centimeters', ...
                 'PaperSize',[sp.W sp.H],'Tag',sp.tag);
    v   = [sp.Ssum(:); sp.Sdif(:)];
    top = max(v);  lo = pct(v,0.5);                                       % 0.5th percentile: ignore single deep notches
    yl3 = [floor((min(lo,top-40)-2)/10)*10, ceil((top+1)/10)*10];
    ys3 = 20*(1 + (diff(yl3) > 90));
    yl  = {[0 1.12], [-195 195], yl3};
    yt  = {0:0.5:1, -180:90:180, ceil(yl3(1)/ys3)*ys3:ys3:yl3(2)};
    for r = 1:nRow
        for j = 1:nCol
            ax = axes(fig,'Units','centimeters', ...
                      'Position',[mar(1)+(j-1)*(pw+gx), mar(2)+(nRow-r)*(ph+gy), pw, ph]);
            hold(ax,'on');  box(ax,'on');
            col = 0.75*sp.regCol(sp.rid(j),:);
            ylr = yl{r};
            if r == 1
                plot(ax,sp.f,sp.C(:,j),'-','Color',col,'LineWidth',0.6);
                plot(ax,sp.fband,sp.gcrit*[1 1],'--','Color',[0.30 0.30 0.30],'LineWidth',0.6);
            elseif r == 2
                ok = sp.C(:,j) >= sp.gcrit;
                plot(ax,sp.fband,[0 0],'-','Color',[0.78 0.78 0.78],'LineWidth',0.4);
                plot(ax,sp.f(~ok),sp.PH(~ok,j),'.','Color',[0.76 0.76 0.76],'MarkerSize',3);
                plot(ax,sp.f(ok), sp.PH(ok,j), '.','Color',col,'MarkerSize',5);
            else
                hd = plot(ax,sp.f,sp.Sdif(:,j),'-','Color',[0.62 0.62 0.62],'LineWidth',0.7);
                hs = plot(ax,sp.f,sp.Ssum(:,j),'-','Color',[0 0 0],'LineWidth',0.6);
            end
            for e = sp.bandEdges(2:end-1)                                 % band edges
                plot(ax,[e e],ylr,':','Color',[0.45 0.45 0.45],'LineWidth',0.5);
            end
            set(ax,'XScale','log','Layer','top','FontSize',6.5,'LineWidth',0.5);
            xlim(ax,sp.fband);  ylim(ax,ylr);
            set(ax,'XTick',[0.01 0.1 1],'YTick',yt{r});
            if r < nRow
                set(ax,'XTickLabel',{});
            else
                set(ax,'XTickLabel',{'0.01','0.1','1'});
                xlabel(ax,'Frequency (Hz)','FontSize',7);
            end
            if j > 1, set(ax,'YTickLabel',{}); end
            if j == 1
                if r == 1,     ylabel(ax,'Coherence \gamma^{2}_{12}','FontSize',7);
                elseif r == 2, ylabel(ax,'Phase \phi_{12} (deg)','FontSize',7);
                else,          ylabel(ax,'PSD (dB)','FontSize',7);
                end
            end
            if r == 1
                text(ax,0.5,1+0.50/ph,sp.regName{j},'Units','normalized','HorizontalAlignment','center', ...
                     'VerticalAlignment','bottom','FontSize',7.5,'FontWeight','bold', ...
                     'Color',0.75*sp.regCol(sp.rid(j),:),'Clipping','off');
                text(ax,0.5,1+0.10/ph,sprintf('P%d, %.0f W',sp.k(j),sp.P(j)),'Units','normalized', ...
                     'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',6.5, ...
                     'Color',[0.25 0.25 0.25],'Clipping','off');
                for b = 1:numel(sp.bandName)                              % band names inside the top of row 1
                    text(ax,sqrt(sp.bandEdges(b)*sp.bandEdges(b+1)),1.065,sp.bandName{b}, ...
                         'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',5.5, ...
                         'Color',[0.35 0.35 0.35]);
                end
            end
            if r == nRow && j == 1
                legend(ax,[hs hd],{sp.legC,sp.legD},'Box','off','FontSize',5.5, ...
                       'Location','southwest');
            end
        end
    end
    saveFig(fig,outDir,name,false);
end

function drawBandFig(outDir,name,bf)
% band metrics against the power steps (one category per step).  Row 1: RMS of the flow fluctuation in the band as % of
% mean |q| (both channels, log axis); row 2: correlation coefficient between the channels in the band (+1 in phase,
% -1 anti-phase).  Background = regime of the step; the start-up step is grey and drawn with open markers.
    nB   = numel(bf.bandName);  nP = numel(bf.P);  nRow = 2;
    fsz  = 6.5 - 0.5*(nB > 3);                                              % tick-label size: smaller panels, smaller labels
    mar  = [1.55 1.15 0.40 1.30];  gx = 0.55;  gy = 0.50;                 % left bottom right top, gaps (cm)
    pw   = (bf.W - mar(1) - mar(3) - (nB-1)*gx)/nB;
    ph   = (bf.H - mar(2) - mar(4) - gy)/nRow;
    fig  = figure('Units','centimeters','Position',[2 2 bf.W bf.H],'Color','w','PaperUnits','centimeters', ...
                  'PaperSize',[bf.W bf.H],'Tag',bf.tag);
    A    = [bf.A1(:); bf.A2(:)];   A = A(A > 0 & ~isnan(A));
    yTop = 1.05*max(A);  dec = 10^floor(log10(yTop));                      % top of the axis: next 1-2-5 value above the data
    nice = dec*[1 2 5 10];   yTop = nice(find(yTop <= nice,1));
    yl1  = [10^floor(log10(min(A))), yTop];
    if yl1(2) < yl1(1)*9.99, yl1(2) = yl1(1)*10; end                       % at least one decade
    ytk  = 10.^(ceil(log10(yl1(1)*(1-1e-9))):floor(log10(yl1(2)*(1+1e-9))));   % ticks at the powers of ten inside the axis
    ur   = unique(bf.rid,'stable');
    isSt = ismember(1:nP,bf.startup);
    mk   = {'o','s'};
    xt   = 1:3:nP;
    xtl  = arrayfun(@(k) sprintf('%d',round(bf.P(k))),xt,'UniformOutput',false);
    tint = @(rid) 1 - 0.22*(1 - bf.regCol(rid,:));                         % light background colour of a regime
    for r = 1:nRow
        if r == 1, yl = yl1; else, yl = [-1.08 1.08]; end
        for b = 1:nB
            ax = axes(fig,'Units','centimeters', ...
                      'Position',[mar(1)+(b-1)*(pw+gx), mar(2)+(nRow-r)*(ph+gy), pw, ph]);
            hold(ax,'on');  box(ax,'on');
            for i = 1:numel(ur)                                            % regime background
                kk = find(bf.rid == ur(i));
                patch(ax,[min(kk)-0.5 max(kk)+0.5 max(kk)+0.5 min(kk)-0.5],yl([1 1 2 2]),tint(ur(i)),'EdgeColor','none');
            end
            for k = find(isSt)                                             % start-up step
                patch(ax,[k-0.5 k+0.5 k+0.5 k-0.5],yl([1 1 2 2]),[0.87 0.87 0.87],'EdgeColor','none');
            end
            if r == 1
                for c = 1:2
                    if c == 1, y = bf.A1(:,b)'; else, y = bf.A2(:,b)'; end
                    plot(ax,1:nP,y,'-','Color',bf.cCh(c,:),'LineWidth',0.6);
                    plot(ax,find(~isSt),y(~isSt),mk{c},'Color',bf.cCh(c,:),'MarkerFaceColor',bf.cCh(c,:), ...
                         'MarkerSize',3.2,'LineStyle','none');
                    plot(ax,find(isSt),y(isSt),mk{c},'Color',bf.cCh(c,:),'MarkerFaceColor','w', ...
                         'MarkerSize',3.2,'LineStyle','none');
                end
                set(ax,'YScale','log');
            else
                plot(ax,[0.5 nP+0.5],[0 0],'-','Color',[0.45 0.45 0.45],'LineWidth',0.4);
                y = bf.rho(:,b)';
                plot(ax,1:nP,y,'-','Color',[0.2 0.2 0.2],'LineWidth',0.6);
                plot(ax,find(~isSt),y(~isSt),'o','Color',[0.2 0.2 0.2],'MarkerFaceColor',[0.2 0.2 0.2], ...
                     'MarkerSize',3.2,'LineStyle','none');
                plot(ax,find(isSt),y(isSt),'o','Color',[0.2 0.2 0.2],'MarkerFaceColor','w', ...
                     'MarkerSize',3.2,'LineStyle','none');
            end
            set(ax,'Layer','top','FontSize',fsz,'LineWidth',0.5);
            xlim(ax,[0.5 nP+0.5]);  ylim(ax,yl);
            if r == 1
                set(ax,'YTick',ytk,'YTickLabel',arrayfun(@(v) sprintf('%g',v),ytk,'UniformOutput',false));
                set(ax,'XTick',xt,'XTickLabel',{});
                text(ax,0.5,1+0.10/ph,sprintf('%s band, %g-%g Hz',bf.bandName{b},bf.bandEdges(b),bf.bandEdges(b+1)), ...
                     'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom', ...
                     'FontSize',7,'FontWeight','bold','Clipping','off');
            else
                set(ax,'YTick',-1:0.5:1,'XTick',xt,'XTickLabel',xtl);
                xlabel(ax,'Heating power (W)','FontSize',7);
            end
            if b > 1
                set(ax,'YTickLabel',{});
            elseif r == 1
                ylabel(ax,'A (% of mean |q|)','FontSize',7);
            else
                ylabel(ax,'\rho between channels','FontSize',7);
            end
        end
    end
    % key (regimes, start-up step, channels) in the top margin
    axL = axes(fig,'Units','centimeters','Position',[0 0 bf.W bf.H]);
    hold(axL,'on');  axis(axL,'off');  xlim(axL,[0 bf.W]);  ylim(axL,[0 bf.H]);
    yy = bf.H - 0.38;  xs = 0.55;  slot = 2.55;
    for i = 1:numel(ur)
        patch(axL,[xs xs+0.30 xs+0.30 xs],[yy-0.15 yy-0.15 yy+0.15 yy+0.15],tint(ur(i)), ...
              'EdgeColor',[0.5 0.5 0.5],'LineWidth',0.4);
        text(axL,xs+0.40,yy,bf.regName{ur(i)},'FontSize',6.5,'VerticalAlignment','middle');
        xs = xs + slot;
    end
    patch(axL,[xs xs+0.30 xs+0.30 xs],[yy-0.15 yy-0.15 yy+0.15 yy+0.15],[0.87 0.87 0.87], ...
          'EdgeColor',[0.5 0.5 0.5],'LineWidth',0.4);
    text(axL,xs+0.40,yy,'start-up step','FontSize',6.5,'VerticalAlignment','middle');
    xs = xs + slot + 0.4;
    for c = 1:2
        plot(axL,[xs xs+0.45],[yy yy],'-','Color',bf.cCh(c,:),'LineWidth',0.6);
        plot(axL,xs+0.225,yy,mk{c},'Color',bf.cCh(c,:),'MarkerFaceColor',bf.cCh(c,:),'MarkerSize',3.2);
        text(axL,xs+0.58,yy,sprintf('Ch%d',c),'FontSize',6.5,'VerticalAlignment','middle');
        xs = xs + 1.6;
    end
    saveFig(fig,outDir,name,false);
end

function drawMatrixFig(outDir,name,Z,okC,greyTxt,o)
% heat-map matrix: rows = variables, columns = power steps.
%   okC     : cells that are coloured (significant / coherent); the others stay light grey
%   greyTxt : grey cells that still show their value in grey text
    [nR,nC] = size(Z);
    nS   = numel(o.strips);
    rowH = 0.30;                                                   % cm, one strip row
    mar  = [2.55 1.30 2.15 0.62+(nS+1)*rowH];                      % left bottom right top, cm
    axW  = o.W - mar(1) - mar(3);   axH = o.H - mar(2) - mar(4);
    cm2x = nC/axW;                                                 % data units per cm along x
    [fig,ax] = newAxesFig(o.W,o.H,mar);
    cmp = getCmap(o.cmap,256);

    % true-colour image: coloured cells from the colour map, all other cells light grey
    cGrey = [0.93 0.93 0.93];
    nc  = size(cmp,1);
    C   = repmat(reshape(cGrey,[1 1 3]),[nR nC 1]);
    for i = 1:nR
        for j = 1:nC
            if okC(i,j) && ~isnan(Z(i,j))
                ix = round(1 + (min(max(Z(i,j),o.cl(1)),o.cl(2))-o.cl(1))/(o.cl(2)-o.cl(1))*(nc-1));
                C(i,j,:) = reshape(cmp(ix,:),[1 1 3]);
            end
        end
    end
    image(ax,C);
    set(ax,'YDir','reverse','Layer','top','TickLength',[0 0],'FontSize',6.5);
    xlim(ax,[0.5 nC+0.5]);  ylim(ax,[0.5 nR+0.5]);
    for j = 0.5:1:nC+0.5, plot(ax,[j j],[0.5 nR+0.5],'-','Color','w','LineWidth',0.35); end
    for i = 0.5:1:nR+0.5, plot(ax,[0.5 nC+0.5],[i i],'-','Color','w','LineWidth',0.35); end
    xs = find(diff(o.regimeId) ~= 0) + 0.5;                        % regime boundaries (vertical lines)
    if ~isempty(o.startup), xs = unique([xs, max(o.startup)+0.5]); end
    for xb = xs
        plot(ax,[xb xb],[0.5 nR+0.5],'-','Color',[0.15 0.15 0.15],'LineWidth',0.6);
    end
    for g = 1:size(o.grp,1)-1                                      % thicker line between the channel blocks
        yb = o.grp(g,2)+0.5;
        plot(ax,[0.5 nC+0.5],[yb yb],'-','Color',[0.15 0.15 0.15],'LineWidth',0.8);
    end
    if nR > o.grp(end,2)
        yb = o.grp(end,2)+0.5;
        plot(ax,[0.5 nC+0.5],[yb yb],'-','Color',[0.15 0.15 0.15],'LineWidth',0.8);
    end

    % numbers inside the cells
    if o.cellText
        for i = 1:nR
            for j = 1:nC
                if isnan(Z(i,j)), continue; end
                if okC(i,j)
                    ix  = round(1 + (min(max(Z(i,j),o.cl(1)),o.cl(2))-o.cl(1))/(o.cl(2)-o.cl(1))*(nc-1));
                    lum = cmp(ix,:)*[0.299; 0.587; 0.114];
                    if lum < 0.5, tc = [1 1 1]; else, tc = [0.1 0.1 0.1]; end
                elseif greyTxt(i,j)
                    tc = [0.60 0.60 0.60];
                else
                    continue;
                end
                text(ax,j,i,o.fmt(Z(i,j)),'HorizontalAlignment','center','VerticalAlignment','middle', ...
                     'FontSize',o.fontCell,'Color',tc);
            end
        end
    end

    set(ax,'XTick',1:nC,'XTickLabel',o.colLab,'YTick',1:nR,'YTickLabel',o.rowTex);
    xlabel(ax,o.xlab,'FontSize',8);
    for g = 1:size(o.grp,1)                                        % channel brackets and names
        plot(ax,0.5-1.28*cm2x*[1 1],[o.grp(g,1)-0.38 o.grp(g,2)+0.38],'-','Color',[0.15 0.15 0.15], ...
             'LineWidth',0.8,'Clipping','off');
        text(ax,0.5-1.62*cm2x,mean(o.grp(g,:)),o.grpName{g},'Rotation',90,'HorizontalAlignment','center', ...
             'VerticalAlignment','middle','FontSize',7.5,'FontWeight','bold','Clipping','off');
    end

    % ---- strips above the matrix: step / regime bar (top) and, optionally, f_p rows ----
    axS = axes(fig,'Units','centimeters','Position',[mar(1) mar(2)+axH+0.05 axW (nS+1)*rowH]);
    hold(axS,'on');  axis(axS,'off');
    xlim(axS,[0.5 nC+0.5]);  ylim(axS,[0 nS+1]);
    for m = 1:nS                                                   % f_p rows, bottom to top
        s = o.strips(nS-m+1);
        for k = 1:nC
            if s.na(k), fc = [0.90 0.90 0.90]; else, fc = [0.97 0.97 0.97]; end
            patch(axS,[k-0.5 k+0.5 k+0.5 k-0.5],[m-1 m-1 m m],fc,'EdgeColor','w','LineWidth',0.5);
            if s.na(k)
                txt = 'n.r.';  tcol = [0.55 0.55 0.55];
            else
                if s.vals(k) < 0.1, txt = sprintf('%.3f',s.vals(k)); txt = txt(2:end);
                else, txt = sprintf('%.2f',s.vals(k)); end
                tcol = [0.1 0.1 0.1];
            end
            text(axS,k,m-0.5,txt,'HorizontalAlignment','center','VerticalAlignment','middle', ...
                 'FontSize',max(4.5,o.fontCell-0.8),'Color',tcol);
        end
        text(axS,0.5-0.12*cm2x,m-0.5,s.label,'HorizontalAlignment','right','VerticalAlignment','middle', ...
             'FontSize',6,'Clipping','off');
    end
    for k = 1:nC                                                   % regime bar
        if any(k == o.startup)
            fc = [0.78 0.78 0.78];
        else
            fc = o.regCol(o.regimeId(k),:);
        end
        patch(axS,[k-0.5 k+0.5 k+0.5 k-0.5],[nS nS nS+1 nS+1],fc,'EdgeColor','w','LineWidth',0.5);
        text(axS,k,nS+0.5,sprintf('%d',k),'HorizontalAlignment','center','VerticalAlignment','middle', ...
             'FontSize',5.5,'FontWeight','bold','Color','w');
    end
    ur = unique(o.regimeId,'stable');
    for i = 1:numel(ur)
        kk = find(o.regimeId == ur(i) & ~ismember(1:nC,o.startup));
        if isempty(kk), continue; end
        text(axS,mean(kk),nS+1.18,o.regName{ur(i)},'HorizontalAlignment','center','VerticalAlignment','bottom', ...
             'FontSize',6.5,'FontWeight','bold','Color',0.75*o.regCol(ur(i),:),'Clipping','off');
    end
    text(axS,0.5-0.12*cm2x,nS+0.5,'Step','HorizontalAlignment','right','VerticalAlignment','middle', ...
         'FontSize',6,'Clipping','off');

    drawCbar(fig,cmp,o.cl,o.cbTicks,o.cbTickLab,o.cbLab,[mar(1)+axW+0.35 mar(2) 0.28 axH]);
    saveFig(fig,outDir,name,false);
end

function writeMatCsv(file,M,rowNames,colNames)
% numeric matrix with row / column labels as plain CSV (opens directly in Excel / Origin); NaN -> empty cell
    fid = fopen(file,'wt');
    if fid < 0, error('Cannot write %s',file); end
    fprintf(fid,'variable,%s\n',strjoin(colNames,','));
    for i = 1:size(M,1)
        s = rowNames{i};
        for j = 1:size(M,2)
            if isnan(M(i,j)), s = [s ',']; else, s = [s sprintf(',%.6g',M(i,j))]; end %#ok<AGROW>
        end
        fprintf(fid,'%s\n',s);
    end
    fclose(fid);
end

function writeTableCsv(file,names,M)
% columns with names as plain CSV; NaN -> empty cell
    fid = fopen(file,'wt');
    if fid < 0, error('Cannot write %s',file); end
    fprintf(fid,'%s\n',strjoin(names,','));
    for i = 1:size(M,1)
        s = '';
        for j = 1:size(M,2)
            if isnan(M(i,j)), t = ''; else, t = sprintf('%.6g',M(i,j)); end
            if j == 1, s = t; else, s = [s ',' t]; end %#ok<AGROW>
        end
        fprintf(fid,'%s\n',s);
    end
    fclose(fid);
end

function saveFig(fig,outDir,name,isMap)
% PNG is always written; PDF has fallbacks so one failing format never stops the script.
% isMap = true -> the PDF embeds a 600 dpi bitmap (a vector PDF of an interpolated surface is huge).
    if nargin < 4, isMap = false; end
    try
        exportgraphics(fig,fullfile(outDir,[name '.png']),'Resolution',600);
    catch ME
        warning('PNG export failed: %s',ME.message);
    end
    try
        if isMap
            exportgraphics(fig,fullfile(outDir,[name '.pdf']),'ContentType','image','Resolution',600);
        else
            exportgraphics(fig,fullfile(outDir,[name '.pdf']),'ContentType','vector');
        end
    catch ME1
        warning('exportgraphics PDF failed (%s). Trying print ...',ME1.message);
        try
            fig.PaperUnits = 'centimeters';
            fig.PaperSize  = fig.Position(3:4);
            fig.PaperPosition = [0 0 fig.Position(3:4)];
            print(fig,fullfile(outDir,[name '_print.pdf']),'-dpdf','-r600');
        catch ME2
            warning('PDF skipped (%s). Use the PNG instead.',ME2.message);
        end
    end
end

function makeScalable(fig)
% Called at the very end of the script, AFTER all PNG / PDF files are written, so the exported figures keep
% their exact paper layout.  From here on the open figure window follows its size:
%   - axes use relative (normalized) positions, so the layout fills the window;
%   - all fonts are tied to the axes height (FontUnits normalized);
%   - line widths, marker sizes and the legend font are rescaled in rescaleFig (SizeChangedFcn).
% Everything is wrapped in try / catch: if a release does not support a step, the figure just stays as it is.
    try
        hA = findall(fig,'Type','axes');
        set(hA,'Units','normalized');
        set(hA,'FontUnits','normalized');
        set(findall(fig,'Type','text'),'FontUnits','normalized');
        pos = get(fig,'Position');
        b.h0 = pos(4);
        b.hM = findall(fig,'Type','line');                                    % markers
        b.hL = [hA(:); b.hM(:); findall(fig,'Type','patch')];                 % line widths
        b.hG = findall(fig,'Type','legend');                                  % legend font
        b.lw = zeros(numel(b.hL),1);  for j = 1:numel(b.hL), b.lw(j) = get(b.hL(j),'LineWidth'); end
        b.ms = zeros(numel(b.hM),1);  for j = 1:numel(b.hM), b.ms(j) = get(b.hM(j),'MarkerSize'); end
        b.fs = zeros(numel(b.hG),1);  for j = 1:numel(b.hG), b.fs(j) = get(b.hG(j),'FontSize'); end
        setappdata(fig,'scaleBase',b);
        set(fig,'SizeChangedFcn',@(src,~) rescaleFig(src));
    catch ME
        warning('Window scaling not available (%s). The figure is unchanged.',ME.message);
    end
end

function rescaleFig(fig)
% SizeChangedFcn: line widths, marker sizes and legend font grow / shrink with the window height
    try
        b = getappdata(fig,'scaleBase');
        if isempty(b), return; end
        pos = get(fig,'Position');
        s   = max(0.3,pos(4)/b.h0);
        if ~isempty(b.hL), set(b.hL,{'LineWidth'},num2cell(b.lw*s)); end
        if ~isempty(b.hM), set(b.hM,{'MarkerSize'},num2cell(b.ms*s)); end
        if ~isempty(b.hG), set(b.hG,{'FontSize'},num2cell(b.fs*s)); end
    catch
    end
end
