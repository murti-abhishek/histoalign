"""3x3 homogeneous matrices in the pixel-corner convention, and alignment-file I/O.

Convention: pixel (row i, col j) has its centre at x = j + 0.5, y = i + 0.5. An alignment matrix
maps H&E level-0 px -> Xenium morphology level-0 px; Xenium microns = Xenium px * pixel_size.
"""
from pathlib import Path

import numpy as np


def S(sx, sy=None):
    sy = sx if sy is None else sy
    return np.array([[sx, 0, 0], [0, sy, 0], [0, 0, 1]], float)


def T(tx, ty):
    return np.array([[1, 0, tx], [0, 1, ty], [0, 0, 1]], float)


def to_idx(M):
    """Corner-convention matrix -> OpenCV pixel-index convention."""
    return T(-0.5, -0.5) @ M @ T(0.5, 0.5)


def apply(M, pts):
    pts = np.asarray(pts, float)
    return (pts @ M[:2, :2].T) + M[:2, 2]


def decompose(M):
    """Rotation (deg), scales, shear-ish summary of the 2x2 part."""
    A = M[:2, :2]
    sx, sy = np.linalg.norm(A[:, 0]), np.linalg.norm(A[:, 1])
    rot = np.degrees(np.arctan2(A[1, 0], A[0, 0]))
    return dict(rotation_deg=float(rot), scale_x=float(sx), scale_y=float(sy),
                det=float(np.linalg.det(A)))


def load_matrix(path):
    """Read a 3x3 matrix from a CSV (Xenium Explorer or numpy layout, LF or CRLF)."""
    txt = Path(path).read_text().replace("\r", "").strip().splitlines()
    M = np.array([[float(v) for v in line.split(",")] for line in txt[:3]])
    if M.shape != (3, 3):
        raise ValueError(f"{path}: expected a 3x3 matrix, got shape {M.shape}")
    return M


def explorer_csv_bytes(G):
    """The matrix as Xenium Explorer writes it: CRLF line endings, no newline after the last row.
    Explorer's loader is strict about this layout."""
    G = np.asarray(G, float)
    rws = [",".join(repr(float(v)) for v in G[0]), ",".join(repr(float(v)) for v in G[1]), "0,0,1"]
    return "\r\n".join(rws).encode()


def write_explorer_csv(G, path):
    """Write an H&E level-0 px -> Xenium morphology level-0 px matrix for Xenium Explorer."""
    Path(path).write_bytes(explorer_csv_bytes(G))
