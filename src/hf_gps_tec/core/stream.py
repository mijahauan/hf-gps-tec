"""HfGpsTecSource — wideband I/Q ingest from radiod via ka9q-python.

One channel per Tx frequency, ≈100 kS/s complex I/Q.  Frames are
emitted at the code-period boundary (100 ms / 10,000 samples by
default), timestamped with the RTP-anchored UTC of the frame's first
sample (per the sigmond timing-authority invariant: never use the
host wall clock).

ka9q-python is the only mandatory runtime dependency for live capture.
For tests we use the synthetic `IqFrameSource` protocol below, which
the daemon also accepts so unit tests can feed canned data without
touching the network.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Iterator, Optional, Protocol

import numpy as np


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class IqFrame:
    """One code-period frame of complex I/Q with UTC anchor."""
    frequency_hz: int
    sample_rate_hz: int
    samples: np.ndarray            # complex64, shape (n_samples,)
    timestamp_utc: datetime
    rtp_anchor_ns: Optional[int] = None   # RTP-anchored UTC ns (preferred)
    radiod_id: str = ""


class IqFrameSource(Protocol):
    """Anything that yields IqFrames for one frequency."""

    def frames(self) -> Iterator[IqFrame]: ...
    def close(self) -> None: ...


# ---------------------------------------------------------------------------
# Live source — ka9q-python MultiStream
# ---------------------------------------------------------------------------


@dataclass
class HfGpsTecSource:
    """Live RTP I/Q source for one Tx frequency.

    Lazy-imports ka9q-python so the rest of the package (config,
    contract, tests) can be exercised without it.
    """
    radiod_status: str        # mDNS hostname of the radiod
    frequency_hz: int
    sample_rate_hz: int = 100_000
    filter_guard_hz: int = 1500
    frame_n_samples: int = 10_000
    client_id: str = "hf-gps-tec"
    radiod_id: str = ""
    _stream: object = field(default=None, init=False, repr=False)

    def open(self) -> None:
        try:
            from ka9q import MultiStream  # type: ignore[import-not-found]
        except ImportError as exc:
            raise RuntimeError(
                "ka9q-python is required for live I/Q capture; "
                "install with `uv sync` or run with a stub IqFrameSource."
            ) from exc

        nyquist = self.sample_rate_hz // 2
        low_edge = -(nyquist - self.filter_guard_hz)
        high_edge = +(nyquist - self.filter_guard_hz)

        logger.info(
            "opening ka9q channel: status=%s freq=%d sr=%d filter=±%d Hz",
            self.radiod_status, self.frequency_hz, self.sample_rate_hz,
            nyquist - self.filter_guard_hz,
        )
        self._stream = MultiStream(
            status=self.radiod_status,
            frequency=self.frequency_hz,
            sample_rate=self.sample_rate_hz,
            preset="iq",
            client_id=self.client_id,
            low_edge=low_edge,
            high_edge=high_edge,
        )

    def frames(self) -> Iterator[IqFrame]:
        """Yield one IqFrame per code period."""
        if self._stream is None:
            self.open()
        assert self._stream is not None  # noqa: S101

        buffer = np.empty(0, dtype=np.complex64)
        first_sample_ns: Optional[int] = None

        for chunk in self._stream:  # type: ignore[attr-defined]
            samples = np.asarray(chunk.samples, dtype=np.complex64)
            chunk_ns = getattr(chunk, "rtp_anchor_ns", None)
            if first_sample_ns is None and chunk_ns is not None:
                first_sample_ns = chunk_ns
            buffer = np.concatenate([buffer, samples]) if buffer.size else samples

            while buffer.size >= self.frame_n_samples:
                frame_samples = buffer[: self.frame_n_samples].copy()
                buffer = buffer[self.frame_n_samples :]
                if first_sample_ns is not None:
                    ts = datetime.fromtimestamp(first_sample_ns / 1e9, tz=timezone.utc)
                else:
                    # Fallback: host clock.  Inferior to RTP-anchored — emit a
                    # warning so this doesn't go unnoticed in production.
                    ts = datetime.now(tz=timezone.utc)
                    logger.warning(
                        "no RTP anchor available for frame; using host clock "
                        "(timing-authority invariant violation)"
                    )
                yield IqFrame(
                    frequency_hz=self.frequency_hz,
                    sample_rate_hz=self.sample_rate_hz,
                    samples=frame_samples,
                    timestamp_utc=ts,
                    rtp_anchor_ns=first_sample_ns,
                    radiod_id=self.radiod_id,
                )
                if first_sample_ns is not None:
                    first_sample_ns += int(
                        self.frame_n_samples * 1e9 / self.sample_rate_hz
                    )

    def close(self) -> None:
        if self._stream is not None:
            try:
                self._stream.close()  # type: ignore[attr-defined]
            except Exception:
                logger.exception("ka9q stream close failed")
            self._stream = None
