# hf-gps-tec receiver methodology

This document describes the high-frequency (HF) pseudorandom-noise-
(PRN-) coded beacon network that `hf-gps-tec` is built to receive,
the receiver digital-signal-processing (DSP) chain, the output
schema, and the open specification gaps that prevent the recorder
from locking to real over-the-air signals at this scaffolding stage.

The reference design follows Hysell, Baumgarten, Milla, Valdez &
Kuyeng (2018, *J. Geophys. Res.: Space Physics* 123:6851–6864) §2,
with the network-topology and observables update from Aricoche &
Hysell (2024, *J. Geophys. Res.: Machine Learning & Computation*
1:e2024JH000270).

`README.md` carries the entry-door overview and install commands.
`docs/OVERVIEW.md` carries the project summary, transmit and
receive architecture, data products, and HamSCI / ionospheric-
science rationale.  This file is the deep technical reference and
the spec-gap tracking document.

## Contents

1. [The Hysell network](#1-the-hysell-network)
2. [Signal waveform](#2-signal-waveform)
3. [Receiver DSP chain](#3-receiver-dsp-chain)
4. [Observables and output schema](#4-observables-and-output-schema)
5. [Scientific value as opportunistic ionospheric input](#5-scientific-value-as-opportunistic-ionospheric-input)
6. [Open gaps (what JRO still needs to supply)](#6-open-gaps-what-jro-still-needs-to-supply)
7. [References](#7-references)

---

## 1. The beacon network

Continuous-wave HF transmitters operated under Dr. David Hysell's
leadership at Cornell University.  The first generation of the
network operated in Peru and was documented in Hysell et al.
(2018) and Aricoche & Hysell (2024) — three transmit sites (Ancon,
Sicaya, Ica) feeding six JRO-administered receive sites (Jicamarca,
Huancayo, Mala, La Merced, Barranca, Oroya).  That deployment is
no longer the active transmit infrastructure; the receive-side DSP
chain described in this document is reused unchanged for the new
deployment.

**Current transmit sites** (North American deployment):

| Site                                                | Status      | Approximate location |
|-----------------------------------------------------|-------------|------|
| Fairbanks, Alaska                                   | operational | ≈ 64.84°N, −147.72°E |
| Ithaca, New York (Cornell University)               | planned     | ≈ 42.45°N, −76.47°E  |

Coordinates above are city-centre approximations pending
confirmation from the network operator (see `data/stations.toml`).

Each transmit site radiates **0.5 W continuous power per frequency**
into inverted-V antennas (per Hysell 2018 §2; carry-over from the
Peru deployment unless and until the operator publishes a
revision).  Both transmit frequencies are emitted simultaneously
from each site.

**Receive sites.**  Not centrally administered in the current
deployment.  Any sigmond station with HF reception infrastructure
can run `hf-gps-tec` and contribute observations; the local
station identifies itself via the `[station]` block in the
recorder config.  The Peru deployment used dual-antenna receive
sites (northeast-southwest plus northwest-southeast) to permit
polarization and arrival-angle (interferometric) measurement; this
scaffolding processes a single antenna per site.

Timing at every transmit site is disciplined by the Global
Positioning System (GPS), anchoring the transmit code-epoch grid
to Coordinated Universal Time (UTC) at sub-microsecond accuracy.
For the pseudorange observable to be quantitatively meaningful at
the receiver, the receive-side sample clock should be similarly
GPS-disciplined — typically via a GPS-disciplined oscillator
(GPSDO) feeding the front-end's reference input, as the sibling
`gpsdo-monitor` sigmond client expects.  A free-running receive
clock would still permit Doppler measurement but the absolute
pseudorange would drift.

## 2. Signal waveform

Both frequencies — **2.72 MHz and 3.64 MHz** — carry a
**unique-per-transmitter binary phase-shift-keyed (BPSK) signal**
modulated by a pseudorandom noise (PRN) code.

Per Hysell 2018 §2:

| Parameter                  | Value          | Derived |
|----------------------------|----------------|---------|
| Chip duration              | 10 µs          | Null-to-null bandwidth ≈ 100 kHz |
| Compression ratio          | 10,000         | Code length = 10,000 chips |
| Code repetition period     | 100 ms         | = 10,000 × 10 µs |
| Modulation                 | Binary phase   | BPSK; per-chip phase ∈ {0, π} |
| Code gain per Doppler bin  | 1 × 10⁶ (60 dB)| 10⁴ (chips) × 10² (coherent reps) |

The code per transmitter is **unique**; the receiver discriminates
co-channel transmitters by correlating against each transmitter's
distinct replica in parallel.  See §6 for what still needs to be
specified about the code structure.

## 3. Receiver DSP chain

The reference receiver in Hysell 2018 §2 samples directly at
10 mega-samples per second (MS/s) at an intermediate frequency
(IF) between 2.72 and 3.64 MHz, then digitally down-converts and
decimates.  `hf-gps-tec` instead lets `radiod` (ka9q-radio) own
the in-phase / quadrature (I/Q) down-conversion and decimation:
one ka9q-radio channel per Tx frequency, ≈100 kilo-samples per
second (kS/s) wideband complex I/Q, ±50 kHz around the carrier.
The remaining DSP runs in Python.

The chain matches Hysell §2 stage-for-stage:

```
ka9q-radio channel @ 2.72 MHz (or 3.64 MHz), 100 kS/s complex I/Q
  │
  ▼
Frame to 100 ms blocks (10,000 samples) aligned to Coordinated
Universal Time (UTC) epoch grid
  │   See §6 gap #2 — UTC alignment protocol presumed but unconfirmed.
  ▼
PRN correlator bank (one replica per known Tx on this frequency)
  │   fast Fourier transform (FFT)-based circular cross-correlation:
  │     r_n[k] = IFFT( FFT(rx_frame) · conj(FFT(replica_n)) )
  │   → complex range profile, 10,000 bins × 1500 m/bin
  ▼
Coherent integrator (100 successive range profiles → 10 s)
  │   Stack into a 100 × 10,000 complex matrix
  │   FFT along axis 0 (slow-time) → range-Doppler matrix
  │   Doppler resolution = 1 / 10 s = 0.1 Hz
  │   Doppler ambiguity = 1 / 100 ms = ±5 Hz
  │   Post-coherent gain ≈ 60 dB total (40 dB code + 20 dB Doppler)
  ▼
Incoherent integrator (6 × 10 s power averaging → 1 min)
  │   |range-Doppler|² averaged over 6 coherent windows
  ▼
First-hop detector
  │   Scan range bins outward from short range; find first bin
  │   exceeding (noise_floor + snr_threshold_db) where snr is the
  │   signal-to-noise ratio (SNR).
  │   → pseudorange (km) = bin_index × 1.5
  │   In that range bin, compute Doppler first moment
  │   → Doppler shift (Hz)
  │   Peak power in that bin → amplitude (dB above noise floor)
  ▼
Per-minute record emitted to JSON Lines (JSONL) + HamSCI sink
```

### 3.1 Channel configuration via ka9q-python

Each frequency is opened as a separate `MultiStream` subscription.
The recorder explicitly overrides the `iq`-preset audio filter via
`ensure_channel(low_edge_hz=-50000, high_edge_hz=+50000)` so the full
PRN bandwidth survives — equivalent to the wideband-filter wiring
codar-sounder uses for the CODAR chirp band.

### 3.2 Doppler ambiguity vs ionospheric reality

At 3.64 MHz, ±5 Hz of Doppler ambiguity corresponds to ±200 m/s
line-of-sight velocity, comfortably above any expected ionospheric
Doppler in the equatorial F region (which Hysell 2024 reports as
peaking near 30 m/s during the prereversal enhancement).  The 0.1 Hz
Doppler bin is ≈4 mm/s, finer than needed but cheap.

### 3.3 Code-free detection mode

The DSP chain above describes the **locked mode** that requires
the per-transmitter PRN code.  Until that specification is
supplied by the network operator (§6 Gap 1), the daemon runs in
**code-free mode** instead, implemented in
`core/detect_codeless.py` and `core/codeless_pipeline.py`.

The discriminating property the detector exploits is that every
PRN-coded beacon of this family repeats its waveform exactly
every 100 ms by construction.  The normalised lagged
autocorrelation at lag τ = one code period is therefore:

```
r(τ) = ⟨ s(t) · s*(t + τ) ⟩ / ⟨ |s(t)|² ⟩
```

For Gaussian noise alone this averages to ≈ 0 (with variance
∝ 1/N).  For a periodic-at-τ signal with power P_s in noise
power P_n,

```
|r(τ)| ≈ P_s / (P_s + P_n)
arg r(τ) = −2π · f_d · τ          → Doppler shift f_d
```

A reference autocorrelation at a non-code-period lag (default
137 ms, chosen to be coprime with the code period) gives the
noise floor.  The detection statistic is the dB ratio of the
two magnitudes; default threshold 6 dB.

```
60-s buffer of contiguous I/Q at 100 kS/s   (≈ 6 million samples)
  │
  ▼
lagged_autocorr(buf, lag = 10,000 samples)  ← r(τ = code period)
lagged_autocorr(buf, lag = 13,700 samples)  ← r(τ = reference)
  │
  ▼
autocorr_db = 20·log10(|r_code| / |r_ref|)
doppler_hz  = −arg(r_code) / (2π · 0.1)
snr_estimate_db = 10·log10(|r_code| / (1 − |r_code|))
  │
  ▼
One JSONL record per integration window emitted to
  /var/lib/hf-gps-tec/<radiod>/codeless/YYYY/MM/DD.jsonl
plus an additive row in sigmond's hf_gps_tec_codeless.spots
```

**What this mode gives:**

- Confirmation that a PRN-coded beacon is present on the channel
  — the autocorrelation magnitude is far above the floor when
  any such beacon is reaching the receiver.
- Doppler shift, unambiguous within ±5 Hz (= ±200 m/s
  line-of-sight at 3.64 MHz).
- Rough received-SNR estimate.
- Band-power time series — directly useful as a propagation
  monitor.

**What this mode cannot give:**

- Per-transmitter identification.  All PRN codes share the
  100 ms period, so co-band transmitters sum into a single
  detection statistic.
- Pseudorange (group delay).  That requires a code-phase
  reference, which the autocorrelation does not provide.

**Sensitivity.**  For input SNR = −20 dB (signal is 1 % of total
received power), |r_code| ≈ 0.01 and the noise floor at
N = 6 × 10⁶ samples is ≈ 4 × 10⁻⁴ — about 28 dB of headroom over
threshold.  Marginal-propagation paths (input SNR closer to
−30 dB) yield ~8 dB headroom; longer integration windows
(parameter `codeless_integration_seconds`) buy additional
detection margin at the cost of cadence.

**Mode selection.**  The daemon picks between locked and
code-free modes from `[mode] mode` in the config:

- `auto` — code-free when `correlate.PRN_IS_STUB` is `True`,
  locked otherwise.  This is the default and the recommended
  setting; the daemon automatically upgrades to locked mode the
  moment the real PRN code is wired in.
- `codeless` — always code-free.
- `locked` — always locked (will produce records derived from
  a fake code if the PRN is still stubbed; not useful).

`inventory --json` reports both the configured and resolved mode
per instance under `mode_configured` / `mode_resolved`.

## 4. Observables and output schema

One JSONL record per (transmitter, receiver, frequency) per minute,
written to:

```
/var/lib/hf-gps-tec/<radiod_id>/YYYY/MM/DD.jsonl
```

and, when sigmond's local sink is writable, mirrored as one row in
the `hf_gps_tec.spots` table of `/var/lib/sigmond/sink.db`.

### Fields

| Field            | Meaning |
|------------------|---|
| `time`           | UTC timestamp at the end of the 1-min incoherent window. |
| `tx_id`          | Transmitter site identifier (`ANCON`, `SICAYA`, `ICA`, …). |
| `rx_id`          | This receiver's station identifier. |
| `radiod_id`      | The radiod this recorder bound to. |
| `frequency_hz`   | Centre frequency (2.72e6 or 3.64e6). |
| `pseudorange_km` | Group delay × c / 2, first-hop X-mode (Hysell 2018 §2). |
| `doppler_hz`     | First moment of Doppler spectrum in the first-hop range bin. |
| `amplitude_db`   | Peak power in the first-hop bin, dB above incoherent noise floor. |
| `snr_db`         | First-hop peak SNR. |
| `n_hops`         | 1 (first-hop only at v0.1.0; multi-hop is future work). |
| `lock_quality`   | 0–1 heuristic (peak prominence + slow-time phase consistency). |
| `noise_floor_db` | Estimated incoherent noise floor in the range-Doppler matrix. |
| `processing_version` | `hf-gps-tec` version string. |
| `contract_version`   | `0.7`. |

The schema is upstream-compatible with the `.out.mod` text format
that Hysell's inversion code (`focus.c`) consumes — emitting that
flavour from the same per-minute record is a small additional sink
(deferred to a follow-up).

## 5. Scientific value as opportunistic ionospheric input

The Hysell network meets the criteria for a useful opportunistic HF
propagation source for HamSCI:

1. **Known transmit geometry.**  Each transmit site's location is
   public to sub-degree accuracy and stable.
2. **Known frequencies and codes.**  Two stable carriers per site,
   each with a fixed PRN code.
3. **GPS-disciplined timing.**  Every transmission carries an
   absolute UTC reference at sub-microsecond accuracy.
4. **Continuous operation.**  0.5 W of continuous wave per
   frequency, 24/7.
5. **Rich per-frame telemetry.**  Pseudorange, Doppler, amplitude,
   and (with dual-antenna sites) polarization and arrival angle —
   more per-spot information than WSPR or FT8 carries.
6. **Both endpoints known when a Tx is decoded** — each detection
   nails down a complete great-circle propagation path with known
   geometry on both ends.
7. **Two-frequency diversity per Tx.**  Pseudorange and Doppler are
   different moments of the ionospheric electron-density profile;
   two frequencies give two independent constraints.

The recorder's outputs are sized to feed Hysell's existing regional
inversion (Aricoche & Hysell 2024) directly — group delay, Doppler,
amplitude at 1-min cadence is exactly the input format that
`focus.c` ingests.

## 6. Open gaps (what the network operator still needs to supply)

The receiver chain is fully specified by Hysell 2018 §2.  The
**transmitted waveform** is not — these are the gaps blocking
real-signal **locked-mode** lock.  Code-free detection (§3.3)
does not need any of the items below: it confirms beacon
presence and recovers Doppler without code knowledge, and is the
default operating mode while these gaps remain open.

### Gap 1 — PRN code specification per transmitter

The Hysell papers say only "unique pseudorandom binary phase code
with compression ratio 10,000."  Almost certainly a maximal-length
sequence (m-sequence) from a 14-stage linear-feedback shift register
truncated to 10,000 chips, or a 14-stage Gold code pair, but the
**generator polynomial and seed per transmitter** are not in the
papers.

What we need, per (transmitter, frequency):
- Generator polynomial (e.g. `x^14 + x^13 + x^8 + x^4 + 1`).
- Seed / initial register state.
- Whether the code is the same across the two frequencies of one
  transmitter or different.
- Chip-clock phase polarity convention (i.e. does a `0` bit produce
  phase `0` or phase `π`).

### Gap 2 — UTC code-epoch alignment

Presumed: each 100 ms code period starts on a 100-ms-aligned UTC
tick (with the network being GPS-disciplined this is the natural
choice).  Should be confirmed.

What we need:
- The UTC boundary the code repeats on (e.g. `t_chip0 mod 100 ms = 0`,
  or some other offset).
- Direction of phase advance through the code (forward in time vs
  reversed).

### Gap 3 — Amplitude calibration reference (lower priority)

Hysell 2024 uses amplitude (in dB) but does not specify whether it is
calibrated to an absolute reference, relative to in-band noise, or
referenced to a synthetic peak.  Since `hf-gps-tec` will
report dB above its own noise floor, this only matters when comparing
across receivers; the inversion uses each receiver's series
internally.

### What's *not* a gap

- The receiver DSP chain (Hysell 2018 §2 is fully specified).
- The output schema (specified by what `focus.c` reads).
- The network topology (current sites listed in §1; new sites can
  be added by editing `/etc/hf-gps-tec/stations.toml`).
- The radiod channel configuration (matches the codar-sounder
  wideband-IQ pattern).

Once gaps 1 and 2 are closed, only `core/correlate.py`'s
`generate_prn_code(tx_id, freq_id)` needs to change — every other
stage of the pipeline is waveform-agnostic.

## 7. References

- Hysell, D. L., Baumgarten, Y., Milla, M. A., Valdez, A., & Kuyeng, K.
  (2018).  "Ionospheric Specification and Space Weather Forecasting
  with an HF Beacon Network in the Peruvian Sector."  *J. Geophys.
  Res.: Space Physics* 123, 6851–6864.
  [doi:10.1029/2018JA025648](https://doi.org/10.1029/2018JA025648).
- Aricoche, J. A. & Hysell, D. L. (2024).  "Ionospheric Radio Beacon
  Signal Analysis and Parameter Estimation Using Automatic
  Differentiation."  *J. Geophys. Res.: Machine Learning &
  Computation* 1, e2024JH000270.
  [doi:10.1029/2024JH000270](https://doi.org/10.1029/2024JH000270).
- Hysell, D. L., Milla, M. A., & Vierinen, J. (2016).  "A multistatic
  HF beacon network for ionospheric specification in the Peruvian
  sector."  *Radio Science* 51, 392–401.
  doi:10.1002/2016RS005951.  Earlier network reference.
