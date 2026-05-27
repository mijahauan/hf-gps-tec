"""PRN replica generation + FFT-based cross-correlation.

The correlator itself (FFT-based circular cross-correlation) is fully
implemented and tested.  The PRN code generator is a **stub** —
deterministic ±1 sequence keyed on transmitter ID — pending the
per-transmitter PRN specification from the JRO network operators.

When the real spec arrives, replace `generate_prn_code()` and set
``PRN_IS_STUB = False``.  Nothing else in the pipeline needs to change.

References:
  - Hysell et al. (2018, *JGR Space Physics* 123:6851–6864), §2:
    "unique pseudorandom binary phase code with a compression ratio
    of 10,000."
  - docs/RECEIVER.md §6 — current gap list.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


# ---------------------------------------------------------------------------
# PRN code generator — STUB
# ---------------------------------------------------------------------------

#: True until the real per-Tx PRN spec is wired in.  Contract `validate
#: --json` surfaces this as a warning so operators can't accidentally
#: deploy a stub-built recorder thinking it's live.
PRN_IS_STUB: bool = True


def generate_prn_code(tx_id: str, frequency_hz: int, n_chips: int = 10_000) -> np.ndarray:
    """Return an ``n_chips`` long ±1 sequence representing the PRN code.

    Currently returns a deterministic m-sequence-like ±1 sequence keyed
    on a stable hash of ``(tx_id, frequency_hz)``.  This is **shape-
    correct** but **does not match** any real transmitter — the
    correlator pipeline can be exercised end-to-end with synthetic data
    using these stubs.

    JRO-spec arrival:
      - Replace this function body with the per-Tx generator polynomial
        and seed.
      - Set ``PRN_IS_STUB = False`` above.
      - Update tests in tests/test_correlate.py accordingly.
    """
    # Stable, reproducible per-(tx_id, freq) seed.
    seed = (hash((tx_id.upper(), int(frequency_hz))) & 0xFFFF_FFFF)
    rng = np.random.default_rng(seed)
    bits = rng.integers(0, 2, size=n_chips, dtype=np.int8)
    return np.where(bits == 0, np.int8(-1), np.int8(1))


# ---------------------------------------------------------------------------
# Replica bank
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Replica:
    """One precomputed transmitter replica for FFT-based correlation."""
    tx_id: str
    frequency_hz: int
    chips: np.ndarray                # int8, ±1, shape (n_chips,)
    fft_conj: np.ndarray             # complex64, shape (n_samples,)


class ReplicaBank:
    """Bank of precomputed replicas for one receive frequency.

    The replicas are upsampled (chip-rate → sample-rate) and their
    conjugate FFTs are cached.  Correlation against a received frame
    is then a single multiply + IFFT per replica.
    """

    def __init__(self, n_samples: int, samples_per_chip: int):
        if n_samples % samples_per_chip != 0:
            raise ValueError(
                f"n_samples ({n_samples}) must be a multiple of "
                f"samples_per_chip ({samples_per_chip})"
            )
        self.n_samples = int(n_samples)
        self.samples_per_chip = int(samples_per_chip)
        self.n_chips = self.n_samples // self.samples_per_chip
        self._replicas: dict[tuple[str, int], Replica] = {}

    def add(self, tx_id: str, frequency_hz: int) -> Replica:
        key = (tx_id.upper(), int(frequency_hz))
        if key in self._replicas:
            return self._replicas[key]
        chips = generate_prn_code(tx_id, frequency_hz, n_chips=self.n_chips)
        # Upsample by sample-and-hold across `samples_per_chip` samples.
        upsampled = np.repeat(chips.astype(np.float32), self.samples_per_chip)
        # Precompute conj(FFT(replica)) for FFT-based circular correlation.
        replica = Replica(
            tx_id=key[0],
            frequency_hz=key[1],
            chips=chips,
            fft_conj=np.conj(np.fft.fft(upsampled)).astype(np.complex64),
        )
        self._replicas[key] = replica
        return replica

    def add_many(self, tx_ids: Iterable[str], frequency_hz: int) -> list[Replica]:
        return [self.add(t, frequency_hz) for t in tx_ids]

    def __iter__(self):
        return iter(self._replicas.values())

    def __len__(self) -> int:
        return len(self._replicas)


# ---------------------------------------------------------------------------
# Correlator
# ---------------------------------------------------------------------------


def correlate(rx_frame: np.ndarray, replica: Replica) -> np.ndarray:
    """FFT-based circular cross-correlation of one rx frame with one replica.

    rx_frame : complex64, shape (n_samples,)
    replica  : Replica with fft_conj precomputed.

    Returns a complex64 range profile of shape (n_chips,) — one bin per
    chip (i.e. one bin per range-resolution cell).
    """
    if rx_frame.dtype != np.complex64:
        rx_frame = rx_frame.astype(np.complex64, copy=False)
    if rx_frame.shape[0] != replica.fft_conj.shape[0]:
        raise ValueError(
            f"rx_frame length {rx_frame.shape[0]} != replica length "
            f"{replica.fft_conj.shape[0]}"
        )
    spectrum = np.fft.fft(rx_frame)
    corr_fft = spectrum * replica.fft_conj
    corr = np.fft.ifft(corr_fft).astype(np.complex64)
    # Decimate by samples_per_chip to get one bin per range cell.
    samples_per_chip = corr.shape[0] // replica.chips.shape[0]
    if samples_per_chip == 1:
        return corr
    return corr[::samples_per_chip].copy()


def correlate_bank(rx_frame: np.ndarray, bank: ReplicaBank) -> dict[str, np.ndarray]:
    """Correlate one rx frame against every replica in the bank.

    Returns a dict mapping tx_id → complex range profile.
    """
    spectrum = np.fft.fft(rx_frame.astype(np.complex64, copy=False))
    samples_per_chip = bank.samples_per_chip
    out: dict[str, np.ndarray] = {}
    for rep in bank:
        corr = np.fft.ifft(spectrum * rep.fft_conj).astype(np.complex64)
        if samples_per_chip > 1:
            corr = corr[::samples_per_chip].copy()
        out[rep.tx_id] = corr
    return out
