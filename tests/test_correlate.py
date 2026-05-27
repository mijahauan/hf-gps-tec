"""Correlator tests.

Exercises the FFT-based circular cross-correlation against synthetic
signals built from the stub PRN code generator.  These tests verify
that the pipeline shape is correct independent of which PRN sequence
is in use, so they will keep passing once the real spec replaces the
stub.
"""

from __future__ import annotations

import numpy as np

from hf_gps_tec.core import correlate as cc


def test_replica_bank_dimensions() -> None:
    bank = cc.ReplicaBank(n_samples=10_000, samples_per_chip=1)
    rep = bank.add("FAIRBANKS", 2_720_000)
    assert rep.chips.shape == (10_000,)
    assert rep.fft_conj.shape == (10_000,)
    assert set(np.unique(rep.chips)).issubset({-1, 1})


def test_stub_codes_are_distinct_per_tx() -> None:
    """Different Tx IDs must produce different stub codes — otherwise the
    correlator would lump them together."""
    a = cc.generate_prn_code("FAIRBANKS", 2_720_000)
    b = cc.generate_prn_code("CORNELL", 2_720_000)
    assert a.shape == b.shape == (10_000,)
    assert not np.array_equal(a, b)


def test_correlator_peak_at_zero_lag() -> None:
    """A noise-free copy of the replica must produce a sharp peak at lag 0."""
    n_samples = 10_000
    bank = cc.ReplicaBank(n_samples=n_samples, samples_per_chip=1)
    rep = bank.add("FAIRBANKS", 2_720_000)
    # Transmit a clean copy of the replica.
    tx = rep.chips.astype(np.complex64)
    profile = cc.correlate(tx, rep)
    assert profile.shape == (10_000,)
    peak_bin = int(np.argmax(np.abs(profile)))
    assert peak_bin == 0, f"expected peak at lag 0, got bin {peak_bin}"
    # Peak should sit far above the noise-floor sidelobes.
    peak = np.abs(profile[0])
    median_sidelobe = float(np.median(np.abs(profile[1:])))
    assert peak / max(median_sidelobe, 1e-12) > 50.0


def test_correlator_resolves_delay() -> None:
    """A circularly-shifted replica must produce a peak at the shift bin."""
    n_samples = 10_000
    bank = cc.ReplicaBank(n_samples=n_samples, samples_per_chip=1)
    rep = bank.add("FAIRBANKS", 2_720_000)
    shift = 137
    tx = np.roll(rep.chips.astype(np.complex64), shift)
    profile = cc.correlate(tx, rep)
    peak_bin = int(np.argmax(np.abs(profile)))
    assert peak_bin == shift


def test_correlate_bank_keys_by_tx_id() -> None:
    bank = cc.ReplicaBank(n_samples=10_000, samples_per_chip=1)
    bank.add("FAIRBANKS", 2_720_000)
    bank.add("CORNELL", 2_720_000)
    rx = np.zeros(10_000, dtype=np.complex64)
    profiles = cc.correlate_bank(rx, bank)
    assert set(profiles.keys()) == {"FAIRBANKS", "CORNELL"}
    for prof in profiles.values():
        assert prof.shape == (10_000,)


def test_prn_is_stub_flag_set() -> None:
    """While the PRN spec hasn't landed, this flag must be True so the
    contract surface warns operators."""
    assert cc.PRN_IS_STUB is True
