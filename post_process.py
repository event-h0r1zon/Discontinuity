from __future__ import annotations

import argparse
import glob
import os
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


@dataclass(frozen=True)
class Snapshot:
    t: float
    x: np.ndarray
    rho: np.ndarray
    u: np.ndarray
    p: np.ndarray
    e: np.ndarray


def read_snapshot_time(csv_path: str) -> float:
    with open(csv_path, "r", encoding="utf-8") as f:
        line = f.readline().strip()
    parts = [s.strip() for s in line.split(",")]
    if len(parts) < 2 or parts[0] != "t":
        raise ValueError(f"{csv_path}: first line must be 't,<value>'")
    return float(parts[1])


def load_snapshot(csv_path: str) -> Snapshot:
    t = read_snapshot_time(csv_path)
    # CSV layout:
    # line 1: t,<value>
    # line 2: i,x,rho,u,p,E
    data = np.loadtxt(csv_path, delimiter=",", skiprows=2)
    if data.ndim == 1:
        data = data[None, :]  # handle single-row edge case

    x = data[:, 1]
    rho = data[:, 2]
    u = data[:, 3]
    p = data[:, 4]
    e = data[:, 5]
    return Snapshot(t=t, x=x, rho=rho, u=u, p=p, e=e)


def discover_snapshots(data_dir: str, pattern: str) -> List[Tuple[float, str]]:
    paths = sorted(glob.glob(os.path.join(data_dir, pattern)))
    if not paths:
        raise FileNotFoundError(f"No snapshot files found: {os.path.join(data_dir, pattern)}")

    t_and_paths = [(read_snapshot_time(p), p) for p in paths]
    t_and_paths.sort(key=lambda tp: tp[0])  # sort by time
    return t_and_paths


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def save_static_plot(snap: Snapshot, out_path: str) -> None:
    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(8, 10))
    axes[0].plot(snap.x, snap.rho, lw=2)
    axes[1].plot(snap.x, snap.u, lw=2)
    axes[2].plot(snap.x, snap.p, lw=2)
    axes[3].plot(snap.x, snap.e, lw=2)

    axes[0].set_ylabel("rho")
    axes[1].set_ylabel("u")
    axes[2].set_ylabel("p")
    axes[3].set_ylabel("e")
    axes[3].set_xlabel("x")

    fig.suptitle(f"State at t ≈ {snap.t:.6f}")
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def save_animation(
    snaps: List[Snapshot],
    out_path: str,
    fps: int,
) -> None:
    x = snaps[0].x

    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(8, 10))
    (line_rho,) = axes[0].plot(x, snaps[0].rho, lw=2)
    (line_u,) = axes[1].plot(x, snaps[0].u, lw=2)
    (line_p,) = axes[2].plot(x, snaps[0].p, lw=2)
    (line_e,) = axes[3].plot(x, snaps[0].e, lw=2)

    axes[0].set_ylabel("rho")
    axes[1].set_ylabel("u")
    axes[2].set_ylabel("p")
    axes[3].set_ylabel("e")
    axes[3].set_xlabel("x")

    # Set y-limits from all frames for stable axes
    rho_all = np.concatenate([s.rho for s in snaps])
    u_all = np.concatenate([s.u for s in snaps])
    p_all = np.concatenate([s.p for s in snaps])
    e_all = np.concatenate([s.e for s in snaps])
    axes[0].set_ylim(rho_all.min() * 0.98, rho_all.max() * 1.02)
    axes[1].set_ylim(u_all.min() * 0.98, u_all.max() * 1.02)
    axes[2].set_ylim(p_all.min() * 0.98, p_all.max() * 1.02)
    axes[3].set_ylim(e_all.min() * 0.98, e_all.max() * 1.02)

    title = fig.suptitle(f"t = {snaps[0].t:.6f}")

    def update(i: int):
        s = snaps[i]
        line_rho.set_ydata(s.rho)
        line_u.set_ydata(s.u)
        line_p.set_ydata(s.p)
        line_e.set_ydata(s.e)
        title.set_text(f"t = {s.t:.6f}")
        return (line_rho, line_u, line_p, line_e, title)

    anim = FuncAnimation(fig, update, frames=len(snaps), interval=1000 / fps, blit=False)
    fig.tight_layout()

    # Prefer MP4 via ffmpeg; fallback to GIF via Pillow if ffmpeg isn't available
    ext = os.path.splitext(out_path)[1].lower()
    try:
        if ext == ".mp4":
            from matplotlib.animation import FFMpegWriter

            anim.save(out_path, writer=FFMpegWriter(fps=fps))
        elif ext == ".gif":
            from matplotlib.animation import PillowWriter

            anim.save(out_path, writer=PillowWriter(fps=fps))
        else:
            raise ValueError("Animation output must end with .mp4 or .gif")
    finally:
        plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data", help="Directory containing step_*.csv files")
    ap.add_argument("--pattern", default="step_*.csv", help="Glob pattern for snapshots")
    ap.add_argument("--out-dir", default="data", help="Where to write plots/animations")
    ap.add_argument("--t-target", type=float, default=0.25, help="Target time for static plot")
    ap.add_argument("--stride", type=int, default=1, help="Use every Nth frame for animation")
    ap.add_argument("--fps", type=int, default=30, help="Animation FPS")
    ap.add_argument("--anim", default="animation.mp4", help="Animation filename (.mp4 or .gif)")
    args = ap.parse_args()

    ensure_dir(args.out_dir)

    t_and_paths = discover_snapshots(args.data_dir, args.pattern)
    t_and_paths = t_and_paths[:: max(1, args.stride)]

    # Load all snapshots (simple + fast enough for small 1D runs)
    snaps = [load_snapshot(p) for (_, p) in t_and_paths]

    # Static plot at closest time to target
    times = np.array([s.t for s in snaps])
    idx = int(np.argmin(np.abs(times - args.t_target)))
    snap_closest = snaps[idx]

    static_path = os.path.join(args.out_dir, f"state_t_{snap_closest.t:.6f}.png")
    save_static_plot(snap_closest, static_path)
    print(f"Wrote {static_path}")

    anim_path = os.path.join(args.out_dir, args.anim)
    save_animation(snaps, anim_path, fps=args.fps)
    print(f"Wrote {anim_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())