"""Coarse alignment: orientation search and SIFT + RANSAC similarity fits."""
import cv2
import numpy as np


# --------------------------------------------------------------------------------------------
# Orientation: the 8 elements of the dihedral group acting on an image of size (h, w)
# --------------------------------------------------------------------------------------------
def dihedral(img, k, flip):
    """Returns (oriented image, M) with M mapping original corner coords -> oriented corner coords."""
    h, w = img.shape[:2]
    M = np.eye(3)
    out = img
    if flip:
        out = out[:, ::-1]
        M = np.array([[-1, 0, w], [0, 1, 0], [0, 0, 1]], float) @ M
    for _ in range(k % 4):  # np.rot90 (counter-clockwise): (x, y) -> (y, W - x)
        W = out.shape[1]
        out = np.rot90(out)
        M = np.array([[0, 1, 0], [-1, 0, W], [0, 0, 1]], float) @ M
    return np.ascontiguousarray(out), M


def sift_matches(a_u8, b_u8, nfeat=20000, ratio=0.8):
    sift = cv2.SIFT_create(nfeatures=nfeat)
    ka, da = sift.detectAndCompute(a_u8, None)
    kb, db = sift.detectAndCompute(b_u8, None)
    if da is None or db is None or len(ka) < 3 or len(kb) < 3:
        return np.zeros((0, 2)), np.zeros((0, 2))
    m = cv2.BFMatcher(cv2.NORM_L2).knnMatch(da, db, k=2)
    good = [p[0] for p in m if len(p) == 2 and p[0].distance < ratio * p[1].distance]
    pa = np.array([ka[g.queryIdx].pt for g in good]) + 0.5  # -> corner convention
    pb = np.array([kb[g.trainIdx].pt for g in good]) + 0.5
    return pa, pb


def fit_similarity(src, dst, thr):
    if len(src) < 3:
        return None, 0
    A, inl = cv2.estimateAffinePartial2D(src, dst, method=cv2.RANSAC, ransacReprojThreshold=thr,
                                         maxIters=10000, confidence=0.999)
    if A is None:
        return None, 0
    M = np.vstack([A, [0, 0, 1]])
    s = np.sqrt(abs(np.linalg.det(M[:2, :2])))
    if not (0.8 < s < 1.25):
        return M, 0
    return M, int(inl.sum())
