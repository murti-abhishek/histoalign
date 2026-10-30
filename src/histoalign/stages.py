"""The four pipeline stages. Each reads and writes its results under `opts.out` (cached in `out/cache`)."""
import json
import re
import time
from pathlib import Path

import cv2
import numpy as np

from ._util import HistoalignError, log
from .coarse import dihedral, fit_similarity, sift_matches
from .combine import benchmark_points, fit_global_from_regions, load_refined, pct, piecewise_to_he, region_lookup
from .io import Pyramid, he_mpp, image_stem, xenium_info
from .layout import choose_layout, regions_from_blocks, regions_from_components, section_mask
from .nuclei import qc_region_nuclei, xenium_nuclei
from .preprocess import clahe, dapi_prep, hematoxylin, ncc, norm01, optical_density, overlay, tissue_mask, warp
from .refine import he_onto_grid, mi_affine, phase_shift, shift_img, summarize, tile_residuals, tiles_absolute
from .transforms import S, apply, decompose, write_explorer_csv


# --------------------------------------------------------------------------------------------
# Stage 1: global alignment + tissue regions
# --------------------------------------------------------------------------------------------
def stage_global(a):
    cache = a.out / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    px, dapi_path, exp = xenium_info(a.xenium)
    mpp = he_mpp(a.he, a.he_mpp)
    log(f"Xenium pixel size {px} um | DAPI: {dapi_path.name}")
    log(f"H&E pixel size {mpp:.4f} um | {Path(a.he).name}")
    D = Pyramid(dapi_path, px, is_ome=False, align=a.pyramid_align, name="DAPI")
    Hp = Pyramid(a.he, mpp, align=a.pyramid_align, name="H&E")
    log(f"DAPI L0 {D.shapes[0]}  ({D.shapes[0][1] * px / 1000:.1f} x {D.shapes[0][0] * px / 1000:.1f} mm)")
    log(f"H&E  L0 {Hp.shapes[0][:2]}  ({Hp.shapes[0][1] * mpp / 1000:.1f} x {Hp.shapes[0][0] * mpp / 1000:.1f} mm)")

    r = a.coarse_um
    dimg, Md = D.grid(r)
    hrgb, Mh = Hp.grid(r)
    dn = dapi_prep(dimg)
    hh = hematoxylin(hrgb)
    hod = optical_density(hrgb)
    hmask = hod > 0.08
    hn = norm01(hh, hmask, 1, 99.5)
    del hrgb, hh
    log(f"coarse grid {r} um/px: DAPI {dn.shape}, H&E {hn.shape}")

    # ---- orientation: SIFT is rotation-invariant, so only the mirror flip has to be searched;
    #      the similarity fit recovers any rotation angle (not just multiples of 90 degrees) ----
    f = 0.5
    ds_ = cv2.resize(dn, None, fx=f, fy=f, interpolation=cv2.INTER_AREA)
    hs_ = cv2.resize(hn, None, fx=f, fy=f, interpolation=cv2.INTER_AREA)
    du8 = clahe(ds_)
    scores = []
    for flip in (False, True):
        ho, Oh = dihedral(hs_, 0, flip)
        pa, pb = sift_matches(clahe(ho), du8, nfeat=a.sift_features)
        M, n = fit_similarity(pa, pb, thr=a.ransac_um / (r / f))
        rot = decompose(M)["rotation_deg"] if M is not None else None
        scores.append(dict(flip=flip, inliers=n, matches=len(pa), rotation_deg=rot,
                           M_search=None if M is None else (M @ Oh).tolist(),
                           _pa=apply(np.linalg.inv(Oh), pa) if len(pa) else pa, _pb=pb))
        log(f"  flip={flip}: {n} RANSAC inliers / {len(pa)} matches"
            + (f", rotation {rot:.1f} deg" if rot is not None and n else ""))
    scores.sort(key=lambda s: -s["inliers"])
    best = scores[0]
    if best["inliers"] < 12:
        raise HistoalignError("No orientation produced a reliable match. Check that the H&E and Xenium are the same "
                 "section, or pass a manual initial transform.")
    ratio = best["inliers"] / max(1, scores[1]["inliers"])
    log(f"orientation: flip={best['flip']}, rotation {best['rotation_deg']:.2f} deg "
        f"({best['inliers']} inliers, {ratio:.1f}x the alternative)")

    # ---- full affine from the search matches (original, unflipped H&E coords), moved to the coarse grid.
    #      Per-piece refinement at higher resolution takes care of the remaining error. ----
    # exact per-axis scale of each half-res image (cv2.resize rounds sizes, so it is not exactly f)
    fh = np.array([hs_.shape[1] / hn.shape[1], hs_.shape[0] / hn.shape[0]])
    fd = np.array([ds_.shape[1] / dn.shape[1], ds_.shape[0] / dn.shape[0]])
    pa, pb = best.pop("_pa") / fh, best.pop("_pb") / fd  # half-res -> coarse grid (corner convention)
    for s_ in scores:
        s_.pop("_pa", None), s_.pop("_pb", None)
    Msearch = S(1 / fd[0], 1 / fd[1]) @ np.array(best["M_search"]) @ S(fh[0], fh[1])
    resid = np.linalg.norm(apply(Msearch, pa) - pb, axis=1)
    keep = resid < (2 * a.ransac_um / r)
    A, inl = cv2.estimateAffine2D(pa[keep], pb[keep], method=cv2.RANSAC,
                                  ransacReprojThreshold=a.ransac_um / r / 2, maxIters=10000, confidence=0.999)
    if A is None or inl.sum() < 12:
        raise HistoalignError("Affine fit failed at the coarse grid.")
    A = np.vstack([A, [0, 0, 1]])
    inl = inl.ravel().astype(bool)
    rms = float(np.sqrt(np.mean(np.sum((apply(A, pa[keep][inl]) - pb[keep][inl]) ** 2, axis=1))) * r)
    log(f"global affine: {int(inl.sum())} inliers of {int(keep.sum())} matches, "
        f"RMS of inlier keypoints {rms:.1f} um")

    G = np.linalg.inv(Md) @ A @ Mh  # H&E L0 px -> DAPI L0 px (any flip is inside A)

    # ---- tissue layout and regions ----
    Mdi = np.linalg.inv(Md)
    m_dapi = tissue_mask(dn, sigma_px=max(1, 8 / r), close_px=max(1, int(30 / r)),
                         min_area_px=a.min_piece_mm2 * 1e6 / r ** 2)
    layout, why = choose_layout(m_dapi, a.layout, a.auto_section_frac)
    log(f"layout: {layout} ({why})")
    if layout == "tma":
        m = m_dapi
        pieces, lab = regions_from_components(m, r)
    else:
        hod_on_dapi = warp(hod, Mh @ np.linalg.inv(G) @ Mdi, dn.shape)
        m, mask_info = section_mask(dimg, m_dapi, hod, hod_on_dapi, r, a)
        log(f"section mask {mask_info['section_mask_mm2']:.1f} mm2 (DAPI-only mask "
            f"{mask_info['dapi_mask_mm2']:.1f} mm2, Xenium imaged field {mask_info['imaged_field_mm2']:.1f} mm2)")
        pieces, lab = regions_from_blocks(m, r, a.block_um, a.min_block_tissue)
    pad = a.pad_um / px
    for p in pieces:
        x0, y0, x1, y1 = p["bbox_grid"]
        c = apply(Mdi, [[x0, y0], [x1, y1]])
        p["core_dapi0"] = [float(c[0, 0]), float(c[0, 1]), float(c[1, 0]), float(c[1, 1])]
        p["bbox_dapi0"] = [float(c[0, 0] - pad), float(c[0, 1] - pad), float(c[1, 0] + pad), float(c[1, 1] + pad)]
    kind = "tissue pieces" if layout == "tma" else f"{a.block_um / 1000:g} mm blocks"
    log(f"{len(pieces)} regions ({kind}), tissue {m.sum() * r * r / 1e6:.1f} mm2")

    np.save(cache / "pieces_labels_coarse.npy", lab.astype(np.int32))
    np.save(cache / "tissue_mask_coarse.npy", m.astype(bool))
    lab_map = {int(p["label"]): p["id"] for p in pieces}

    # ---- whole-slide overlay with the global transform ----
    hw = warp(hn, Mh @ np.linalg.inv(G) @ Mdi, dn.shape)
    ov = overlay(dn, norm01(hw, m_dapi, 1, 99.5))
    sc = min(1.0, 2000 / max(ov.shape[:2]))
    ov = cv2.resize(ov, None, fx=sc, fy=sc, interpolation=cv2.INTER_AREA)
    if layout == "section":
        cnt, _ = cv2.findContours(cv2.resize(m.astype(np.uint8), (ov.shape[1], ov.shape[0]),
                                             interpolation=cv2.INTER_NEAREST), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        cv2.drawContours(ov, cnt, -1, (0, 200, 255), 1)
    for p in pieces:
        x0, y0, x1, y1 = (np.array(p["bbox_grid"]) * sc).astype(int)
        if layout == "section":
            cv2.rectangle(ov, (x0, y0), (x1, y1), (255, 255, 0), 1)
        cv2.putText(ov, p["id"], (x0 + 2, y0 + 14), cv2.FONT_HERSHEY_SIMPLEX, 0.45, (255, 255, 0), 1, cv2.LINE_AA)
    cv2.imwrite(str(a.out / "01_global_overlay.png"), cv2.cvtColor(ov, cv2.COLOR_RGB2BGR))

    res = dict(xenium=str(Path(a.xenium).resolve()), he=str(Path(a.he).resolve()), xenium_px_um=px,
               he_mpp=mpp, dapi_path=str(dapi_path), coarse_um=r, orientation=best,
               orientation_scores=scores, orientation_margin=ratio, G=G.tolist(),
               G_summary=decompose(G), global_inliers=int(inl.sum()),
               global_keypoint_rms_um=rms, M_dapi0_to_coarse=Md.tolist(), pieces=pieces,
               label_to_piece=lab_map, layout=layout, layout_reason=why,
               block_um=a.block_um if layout == "section" else None,
               tissue_mask_mm2=float(m.sum() * r * r / 1e6),
               dapi_mask_mm2=float(m_dapi.sum() * r * r / 1e6))
    (cache / "global.json").write_text(json.dumps(res, indent=1))
    log("wrote 01_global_overlay.png and cache/global.json")


def stage_refine(a):
    cache = a.out / "cache"
    g = json.loads((cache / "global.json").read_text())
    G = np.array(g["G"])
    px, mpp = g["xenium_px_um"], g["he_mpp"]
    D = Pyramid(g["dapi_path"], px, is_ome=False, align=a.pyramid_align, name="DAPI")
    Hp = Pyramid(a.he, mpp, align=a.pyramid_align, name="H&E")
    pieces = g["pieces"]
    sel = range(len(pieces))
    if a.pieces:
        lo, _, hi = a.pieces.partition(":")
        sel = range(int(lo), int(hi) if hi else int(lo) + 1)
    (a.out / "pieces").mkdir(exist_ok=True)
    for i in sel:
        if i >= len(pieces):
            break
        p = pieces[i]
        fout = cache / f"piece_{p['id']}.json"
        if fout.exists() and not a.force:
            log(f"{p['id']}: cached, skipping")
            continue
        t0 = time.time()
        bb = p["bbox_dapi0"]
        # ---- refine at fine_um ----
        r = a.fine_um
        dimg, Mdg = D.grid(r, bb)
        dn = dapi_prep(dimg)
        dmask = tissue_mask(dn, sigma_px=max(1, 6 / r), close_px=max(1, int(20 / r)), min_area_px=2000 / r ** 2)
        hrgb = he_onto_grid(Hp, G, Mdg, dn.shape, r)
        hn = norm01(hematoxylin(hrgb), dmask, 1, 99.5)
        ncc_g = ncc(cv2.GaussianBlur(dn, (0, 0), 2), cv2.GaussianBlur(hn, (0, 0), 2), dmask)
        # translation init by phase correlation on blurred images, accepted only if it helps
        db, hb = cv2.GaussianBlur(dn, (0, 0), 3), cv2.GaussianBlur(hn, (0, 0), 3)
        sx, sy = phase_shift(db, hb)
        init = (0.0, 0.0)
        if np.hypot(sx, sy) * r < a.max_shift_um and ncc(db, shift_img(hb, sx, sy), dmask) > ncc(db, hb, dmask):
            init = (sx, sy)
        Mc, mi0, mi1, stop = mi_affine(dn, hn, dmask, r, init)
        # Mc maps DAPI-grid corner coords -> current H&E-on-grid corner coords
        Ginv_new = np.linalg.inv(G) @ np.linalg.inv(Mdg) @ Mc @ Mdg
        Gp = np.linalg.inv(Ginv_new)
        ctr = apply(np.linalg.inv(Mdg), [[dn.shape[1] / 2, dn.shape[0] / 2]])[0]
        move_um = float(np.linalg.norm(apply(np.linalg.inv(Gp), [ctr])[0] - apply(np.linalg.inv(G), [ctr])[0]) * mpp)
        dec = decompose(Mc)
        accepted = (mi1 < mi0) and (move_um < a.max_shift_um) and abs(dec["det"] - 1) < 0.08
        hn_r = norm01(hematoxylin(he_onto_grid(Hp, Gp, Mdg, dn.shape, r)), dmask, 1, 99.5)
        ncc_r = ncc(cv2.GaussianBlur(dn, (0, 0), 2), cv2.GaussianBlur(hn_r, (0, 0), 2), dmask)
        if accepted and ncc_r < ncc_g:
            accepted = False
        Guse = Gp if accepted else G
        # overlay: global | refined
        ov = np.concatenate([overlay(dn, hn), np.full((dn.shape[0], 8, 3), 255, np.uint8), overlay(dn, hn_r)], axis=1)
        sc = min(1.0, 1400 / ov.shape[1])
        ov = cv2.resize(ov, None, fx=sc, fy=sc, interpolation=cv2.INTER_AREA)
        cv2.imwrite(str(a.out / "pieces" / f"{p['id']}_overlay_global_vs_refined.png"), cv2.cvtColor(ov, cv2.COLOR_RGB2BGR))

        # ---- residual QC at qc_um (nuclear scale) ----
        rq = a.qc_um
        dq, Mdq = D.grid(rq, bb)
        dqn = dapi_prep(dq)
        qmask = cv2.resize(dmask.astype(np.uint8), (dqn.shape[1], dqn.shape[0]), interpolation=cv2.INTER_NEAREST) > 0
        res = {}
        for name, GG in (("global", G), ("refined", Gp)):
            hq = norm01(hematoxylin(he_onto_grid(Hp, GG, Mdq, dqn.shape, rq)), qmask, 1, 99.5)
            tiles = tiles_absolute(tile_residuals(dqn, hq, qmask, rq, a.tile_um), Mdq, rq, px,
                                   p.get("core_dapi0") if p.get("kind") == "block" else None)
            res[name] = dict(summary=summarize([t for t in tiles if t["in_core"]]), tiles=tiles)
            if name == "refined":
                zx0, zy0 = dqn.shape[1] // 2 - 200, dqn.shape[0] // 2 - 200
                z = overlay(dqn[zy0:zy0 + 400, zx0:zx0 + 400], hq[zy0:zy0 + 400, zx0:zx0 + 400])
                cv2.imwrite(str(a.out / "pieces" / f"{p['id']}_zoom_refined_{int(400 * rq)}um.png"),
                            cv2.cvtColor(cv2.resize(z, None, fx=2, fy=2, interpolation=cv2.INTER_NEAREST), cv2.COLOR_RGB2BGR))
        out = dict(id=p["id"], bbox_dapi0=bb, G_refined=Gp.tolist(), accepted=bool(accepted),
                   transform_used="refined" if accepted else "global", G_used=Guse.tolist(),
                   refine_vs_global=dec, piece_moved_um=move_um, mi_before=mi0, mi_after=mi1,
                   ncc_global=ncc_g, ncc_refined=ncc_r, optimizer_stop=stop,
                   residual_global=res["global"]["summary"], residual_refined=res["refined"]["summary"],
                   tiles_global=res["global"]["tiles"], tiles_refined=res["refined"]["tiles"],
                   seconds=round(time.time() - t0, 1))
        fout.write_text(json.dumps(out, indent=1))
        rg, rr = res["global"]["summary"], res["refined"]["summary"]
        f = lambda s: "n/a" if s["median_um"] is None else f"{s['median_um']:.1f} um (p90 {s['p90_um']:.1f}, {s['n_tiles']} tiles)"
        log(f"{p['id']}: moved {move_um:5.1f} um, MI {mi0:.3f}->{mi1:.3f}, {'ACCEPT' if accepted else 'keep global'} | "
            f"residual global {f(rg)} -> refined {f(rr)} [{time.time() - t0:.1f}s]")


def stage_report(a):
    import pandas as pd
    cache = a.out / "cache"
    g = json.loads((cache / "global.json").read_text())
    layout = g.get("layout", "tma")
    section = layout == "section"
    G_sift = np.array(g["G"])
    px, mpp = g["xenium_px_um"], g["he_mpp"]
    pieces = {p["id"]: p for p in g["pieces"]}
    per = load_refined(g, cache)
    unit = "block" if section else "piece"

    # ---- ONE affine refit to all accepted regions. This removes the small systematic error of the
    #      low-res SIFT fit; what remains per region is local deformation (blocks) or piece drift. ----
    G, drift, excluded = fit_global_from_regions(per, pieces, g, cache, robust=section, outlier_um=a.drift_outlier_um)
    if G is None:
        log("fewer than 3 accepted regions: exporting the SIFT affine")
        G = G_sift.copy()
    else:
        log(f"refit global affine on {len(drift)} accepted {unit}s"
            + (f" ({len(excluded)} excluded as outliers: {', '.join(excluded)})" if excluded else "")
            + f"; {unit} drift from it: median {np.median([v[0] for v in drift.values()]):.2f} um, "
              f"max {max(v[1] for v in drift.values()):.2f} um")

    rows = []
    for pid, d in per.items():
        p = pieces[pid]
        rows.append(dict(region=pid, kind=p.get("kind", "piece"), area_mm2=round(p["area_mm2"], 3),
                         x_um=round(np.mean(p["core_dapi0"][0::2]) * px, 1) if "core_dapi0" in p else None,
                         y_um=round(np.mean(p["core_dapi0"][1::2]) * px, 1) if "core_dapi0" in p else None,
                         transform_used=d["transform_used"], excluded_from_global_fit=pid in excluded,
                         moved_vs_sift_um=round(d["piece_moved_um"], 2),
                         rotation_change_deg=round(d["refine_vs_global"]["rotation_deg"], 3),
                         scale_change=round(np.sqrt(abs(d["refine_vs_global"]["det"])), 4),
                         drift_from_global_mean_um=round(drift[pid][0], 2) if pid in drift else None,
                         drift_from_global_max_um=round(drift[pid][1], 2) if pid in drift else None,
                         residual_sift_median_um=d["residual_global"]["median_um"],
                         residual_sift_p90_um=d["residual_global"]["p90_um"],
                         residual_refined_median_um=d["residual_refined"]["median_um"],
                         residual_refined_p90_um=d["residual_refined"]["p90_um"],
                         qc_tiles=d["residual_refined"]["n_tiles"]))
    df = pd.DataFrame(rows)
    df.to_csv(a.out / "per_region_alignment_qc.csv", index=False)
    lookup = region_lookup(g, cache)

    # ---- exports ----
    np.savetxt(a.out / "he_to_xenium_global_alignment.csv", G, delimiter=",", fmt="%.10g")
    # Xenium Explorer's loader is strict: same layout it exports itself (CRLF line endings,
    # no newline after the last row). H&E level-0 px -> Xenium morphology level-0 px.
    write_explorer_csv(G, a.out / f"{image_stem(g['he'])}_imagealignment.csv")
    np.savetxt(a.out / "he_to_xenium_global_alignment_um.csv", S(px) @ G, delimiter=",", fmt="%.10g")
    np.savetxt(cache / "he_to_xenium_global_alignment_sift_initial.csv", G_sift, delimiter=",", fmt="%.10g")
    json.dump({k: dict(kind=pieces[k].get("kind", "piece"), transform_used=v["transform_used"],
                       he0_to_xenium_px=v["G_used"], he0_to_xenium_um=(S(px) @ np.array(v["G_used"])).tolist())
               for k, v in per.items()}, open(a.out / "per_piece_transforms.json", "w"), indent=1)

    # ---- benchmark against an independent alignment (e.g. 10x's *_he_imagealignment.csv) ----
    ref_lines, bench = [], None
    X, what = benchmark_points(g, cache)
    region = lookup(X)
    if a.reference_alignment:
        R = np.loadtxt(a.reference_alignment, delimiter=",")
        err = lambda he: np.linalg.norm(apply(R, he) - X, axis=1) * px  # Xenium um
        bench = pd.DataFrame([dict(transform="initial SIFT affine", **pct(err(apply(np.linalg.inv(G_sift), X)))),
                              dict(transform="refit global affine (exported)", **pct(err(apply(np.linalg.inv(G), X)))),
                              dict(transform=f"per-{unit} transforms", **pct(err(piecewise_to_he(X, G, per, region))))])
        bench.insert(1, "points", what)
        bench.to_csv(a.out / "benchmark_vs_reference.csv", index=False)
        ref_lines = ["", f"Agreement with `{Path(a.reference_alignment).name}` at {what} (um):", "",
                     "| transform | median | p90 | p99 | max |", "|---|---|---|---|---|"] + \
                    [f"| {r.transform} | {r.median_um:.2f} | {r.p90_um:.2f} | {r.p99_um:.2f} | {r.max_um:.2f} |"
                     for r in bench.itertuples()]
        for r in bench.itertuples():
            log(f"vs reference, {r.transform}: median {r.median_um:.2f}, p90 {r.p90_um:.2f}, max {r.max_um:.2f} um")

    # ---- every Xenium cell centroid in H&E level-0 px ----
    cells_pq = Path(g["xenium"]) / "cells.parquet"
    n_cells = None
    if cells_pq.exists():
        cells = pd.read_parquet(cells_pq, columns=["cell_id", "x_centroid", "y_centroid"])
        xy0 = cells[["x_centroid", "y_centroid"]].to_numpy() / px
        reg = lookup(xy0)
        he_g = apply(np.linalg.inv(G), xy0)
        he_r = piecewise_to_he(xy0, G, per, reg)
        rec = he_g if section else he_r  # one section: the single affine; TMA: each core's own transform
        out = pd.DataFrame(dict(cell_id=cells["cell_id"].astype(str), region=reg,
                                he_x_px=rec[:, 0], he_y_px=rec[:, 1],
                                he_x_px_global=he_g[:, 0], he_y_px_global=he_g[:, 1],
                                he_x_px_region=he_r[:, 0], he_y_px_region=he_r[:, 1]))
        out.to_parquet(a.out / "cells_in_he_coordinates.parquet", index=False)
        n_cells = len(out)
        log(f"mapped {n_cells:,} cells into H&E pixel coordinates")

    # ---- figures ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if len(df):
        d2 = df.dropna(subset=["residual_refined_median_um"])
        fig, ax = plt.subplots(figsize=(max(8, 0.32 * len(d2) + 2), 4.2))
        x = np.arange(len(d2))
        ax.bar(x - 0.2, d2["residual_sift_median_um"], 0.4, label="initial global affine (SIFT)", color="#9aa5b1")
        ax.bar(x + 0.2, d2["residual_refined_median_um"], 0.4, label=f"per-{unit} refined", color="#2b6cb0")
        ax.errorbar(x + 0.2, d2["residual_refined_median_um"],
                    yerr=[np.zeros(len(d2)), d2["residual_refined_p90_um"] - d2["residual_refined_median_um"]],
                    fmt="none", ecolor="#1a365d", capsize=2, lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(d2["region"].str.replace("blk_", ""), rotation=90, fontsize=7)
        ax.set_ylabel("residual misalignment (um)\nmedian, whisker = p90")
        ax.set_title(f"H&E vs DAPI residual per {unit} (tile phase correlation, ~{a.qc_um:g} um/px)")
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(a.out / "02_residual_per_region.png", dpi=150)
        plt.close(fig)
    if section and drift:
        lab = np.load(cache / "pieces_labels_coarse.npy")
        Md = np.array(g["M_dapi0_to_coarse"])
        ext = np.r_[apply(np.linalg.inv(Md), [[0, 0]])[0], apply(np.linalg.inv(Md), [[lab.shape[1], lab.shape[0]]])[0]] * px / 1000
        fig, ax = plt.subplots(figsize=(10, 10 * lab.shape[0] / lab.shape[1] + 1))
        ax.imshow(lab > 0, cmap="Greys", alpha=0.2, extent=[ext[0], ext[2], ext[3], ext[1]])
        vmax = max(1.0, max(v[0] for v in drift.values()))
        cm = plt.get_cmap("magma_r")
        for pid, p in pieces.items():
            x0, y0, x1, y1 = np.array(p["core_dapi0"]) * px / 1000
            if pid in drift:
                ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, color=cm(drift[pid][0] / vmax), alpha=0.75))
                ax.text((x0 + x1) / 2, (y0 + y1) / 2, f"{drift[pid][0]:.1f}" + ("*" if pid in excluded else ""),
                        ha="center", va="center", fontsize=8)
            else:
                ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, hatch="//", ec="0.5"))
        sm = plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(0, vmax))
        plt.colorbar(sm, ax=ax, shrink=0.6, label="block drift from the single affine (mean um)")
        ax.set_xlim(ext[0], ext[2]); ax.set_ylim(ext[3], ext[1]); ax.set_aspect("equal")
        ax.set_xlabel("Xenium x (mm)"); ax.set_ylabel("Xenium y (mm)")
        ax.set_title("Local deformation: how far each block's own best affine departs from the exported one\n"
                     "(hatched = block refinement rejected; * = excluded from the fit as an outlier)", fontsize=10)
        fig.tight_layout()
        fig.savefig(a.out / "03_block_drift_map.png", dpi=150)
        plt.close(fig)

    # ---- final-transform QC (written by the `qc` stage) ----
    qc_rows = []
    qf = a.out / "qc_final.json"
    if qf.exists():
        q = json.loads(qf.read_text())
        if np.allclose(np.array(q["G"]), G, atol=1e-6):
            qc_rows.append(("**exported affine (this run)**", q["summary"]))
        else:
            log("qc_final.json is from a different transform; rerun the `qc` stage")
    for qe in sorted(a.out.glob("qc_ext_*.json")):
        q = json.loads(qe.read_text())
        qc_rows.append((f"`{Path(q['alignment']).name}`", q["summary"]))
    qc_lines = []
    if qc_rows:
        nuc_ = qc_rows[0][1].get("method") == "nuclei"
        qc_lines = ["", ("Independent check by nucleus matching: Xenium nucleus centroids vs hematoxylin nuclei detected "
                         f"in the H&E, mutual nearest neighbours within 3 um, median offset per {a.tile_um:.0f} um square (um):"
                         if nuc_ else f"Tile phase-correlation residual, {a.tile_um:.0f} um tiles (um):"), "",
                    "| transform | median | p90 | p99 | systematic offset (dx, dy) | nucleus pairs |",
                    "|---|---|---|---|---|---|"]
        for name, sm in qc_rows:
            off_ = (f"{sm['systematic_offset_um']:.2f} ({sm['mean_dx_um']:+.2f}, {sm['mean_dy_um']:+.2f})"
                    if "systematic_offset_um" in sm else "")
            qc_lines.append(f"| {name} | {sm['median_um']:.2f} | {sm['p90_um']:.2f} | {sm.get('p99_um', float('nan')):.2f} | "
                            f"{off_} | {sm.get('n_matched', '')} |")

    # ---- report ----
    def med(col):
        v = df[col].dropna() if len(df) else []
        return f"{v.median():.2f}" if len(v) else "n/a"
    o = g["orientation"]
    dG = decompose(G)
    shift = np.linalg.norm(apply(np.linalg.inv(G), X) - apply(np.linalg.inv(G_sift), X), axis=1) * mpp
    acc = int((df["transform_used"] == "refined").sum()) if len(df) else 0
    if section:
        layout_line = (f"- Layout: **one tissue section** ({g.get('layout_reason', '')}). Tissue mask "
                       f"{g['tissue_mask_mm2']:.1f} mm2 (DAPI-only mask {g['dapi_mask_mm2']:.1f} mm2), tiled into "
                       f"{len(pieces)} blocks of {g['block_um'] / 1000:g} mm")
        drift_line = ("- Local deformation (block drift from the single affine): "
                      + (f"median {np.median([v[0] for v in drift.values()]):.2f} um, "
                         f"max {max(v[1] for v in drift.values()):.2f} um" if drift else "n/a")
                      + (f"; {len(excluded)} outlier blocks excluded from the fit ({', '.join(excluded)})" if excluded else ""))
    else:
        layout_line = f"- Layout: **separate tissue pieces** ({g.get('layout_reason', 'tma')}); {len(pieces)} pieces"
        drift_line = ("- Piece drift from the single affine (what a global-only alignment would get wrong): "
                      + (f"median {np.median([v[0] for v in drift.values()]):.1f} um, "
                         f"max {max(v[1] for v in drift.values()):.1f} um" if drift else "n/a"))
    lines = [
        "# H&E -> Xenium registration report", "",
        f"- Xenium: `{g['xenium']}` (pixel {px} um)",
        f"- H&E: `{g['he']}` (pixel {mpp:.4f} um)",
        layout_line,
        f"- Orientation: mirror flip={o['flip']}, rotation {g['G_summary']['rotation_deg']:.2f} deg "
        f"({o['inliers']} SIFT inliers; the other flip gave {g['orientation_scores'][1]['inliers']})",
        f"- Initial affine (SIFT + RANSAC at {g['coarse_um']} um/px): {g['global_inliers']} inlier keypoints, "
        f"RMS {g['global_keypoint_rms_um']:.1f} um",
        f"- {acc} of {len(df)} {unit}s accepted the mutual-information refinement",
        f"- Exported affine (refit on the accepted {unit}s): rotation {dG['rotation_deg']:.3f} deg, scale "
        f"{dG['scale_x']:.5f}/{dG['scale_y']:.5f}; moves points by median {np.median(shift):.2f} um vs the SIFT fit",
        drift_line,
        f"- Tile residual, median over {unit}s: initial SIFT {med('residual_sift_median_um')} um -> "
        f"{unit} refined {med('residual_refined_median_um')} um (self-consistency only: this measure pixel-locks below ~1 um)",
        *ref_lines, *qc_lines,
        "", "Outputs:",
        f"- `{image_stem(g['he'])}_imagealignment.csv`: the exported affine in Xenium Explorer's format "
        "(H&E level-0 px -> Xenium morphology level-0 px); same matrix in `he_to_xenium_global_alignment.csv` "
        "and in microns in `he_to_xenium_global_alignment_um.csv`.",
        f"- `per_piece_transforms.json`: each {unit}'s own transform (px and um).",
        "- `cells_in_he_coordinates.parquet`: every Xenium cell centroid in H&E level-0 px"
        + (f" ({n_cells:,} cells); `he_x_px/he_y_px` use "
           + ("the single affine (one section)" if section else "each piece's own transform")
           + ", both variants are also included." if n_cells else "."),
        "- `per_region_alignment_qc.csv`, `02_residual_per_region.png`, `01_global_overlay.png`, `pieces/`"
        + (", `03_block_drift_map.png`" if section else "") + (", `benchmark_vs_reference.csv`" if bench is not None else ""),
        "", "Notes: the per-block tile residual compares the same two images the refinement optimised, so it "
        "only shows the fit converged. The nucleus-matching check is the accuracy estimate: it compares "
        "Xenium's own nucleus segmentation with nuclei found independently in the H&E, so it needs no "
        "hand-placed landmarks; per-nucleus scatter (~0.5 um) averages out over hundreds of nuclei per square. "
        "Neither replaces a look at the overlays.",
    ]
    (a.out / "REPORT.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def stage_qc(a):
    import pandas as pd
    cache = a.out / "cache"
    g = json.loads((cache / "global.json").read_text())
    if a.qc_alignment:  # QC someone else's matrix (manual, vendor) with exactly the same measurement
        G = np.loadtxt(a.qc_alignment, delimiter=",")
        tag = "ext_" + re.sub(r"[^A-Za-z0-9]+", "_", a.qc_name or Path(a.qc_alignment).stem)
    else:
        G = np.loadtxt(a.out / "he_to_xenium_global_alignment.csv", delimiter=",")
        tag = "final"
    px, mpp = g["xenium_px_um"], g["he_mpp"]
    nuc = xenium_nuclei(g, cache) if a.qc_method in ("auto", "nuclei") else None
    if a.qc_method == "nuclei" and nuc is None:
        raise HistoalignError("--qc-method nuclei needs nucleus_boundaries.parquet in the Xenium folder")
    method = "nuclei" if nuc is not None else "tiles"
    D = Pyramid(g["dapi_path"], px, is_ome=False, align=a.pyramid_align, name="DAPI") if method == "tiles" else None
    Hp = Pyramid(g["he"], mpp, align=a.pyramid_align, name="H&E")
    pieces = g["pieces"]
    sel = range(len(pieces))
    if a.pieces:
        lo, _, hi = a.pieces.partition(":")
        sel = range(int(lo), int(hi) if hi else int(lo) + 1)
    rq = a.qc_um
    fresh = lambda d: np.allclose(d["G"], G, atol=1e-6) and d.get("method") == method
    for i in sel:
        if i >= len(pieces):
            break
        p = pieces[i]
        f = cache / f"qc_{tag}_{p['id']}.json"
        if f.exists() and not a.force and fresh(json.loads(f.read_text())):
            continue
        t0 = time.time()
        if method == "nuclei":
            tiles, s = qc_region_nuclei(Hp, G, p, nuc, px, a.tile_um)
            msg = (f"{s['n_matched']:,} nucleus pairs, offset ({s['median_dx_um']:+.2f}, {s['median_dy_um']:+.2f}) um"
                   if s["n_matched"] else "too few nuclei")
        else:
            dq, Mdq = D.grid(rq, p["bbox_dapi0"])
            dqn = dapi_prep(dq)
            qmask = tissue_mask(dqn, sigma_px=max(1, 6 / rq), close_px=max(1, int(20 / rq)), min_area_px=2000 / rq ** 2)
            hq = norm01(hematoxylin(he_onto_grid(Hp, G, Mdq, dqn.shape, rq)), qmask, 1, 99.5)
            tiles = tiles_absolute(tile_residuals(dqn, hq, qmask, rq, a.tile_um), Mdq, rq, px,
                                   p.get("core_dapi0") if p.get("kind") == "block" else None)
            s = summarize([t for t in tiles if t["in_core"]])
            msg = f"{s['n_tiles']} tiles, median " + ("n/a" if s["median_um"] is None else f"{s['median_um']:.2f} um")
        f.write_text(json.dumps(dict(id=p["id"], G=G.tolist(), method=method, summary=s, tiles=tiles)))
        log(f"qc {p['id']}: {msg} [{time.time() - t0:.1f}s]")
    done = [json.loads((cache / f"qc_{tag}_{p['id']}.json").read_text()) for p in pieces
            if (cache / f"qc_{tag}_{p['id']}.json").exists()]
    done = [d for d in done if fresh(d)]
    if len(done) < len(pieces):
        log(f"qc: {len(done)}/{len(pieces)} regions done with the current transform; run again to finish")
        return
    tiles = pd.DataFrame([dict(region=d["id"], **t) for d in done for t in d["tiles"] if t["in_core"]])
    tiles.to_csv(a.out / f"qc_{tag}_tiles.csv", index=False)
    e = tiles["err_um"].to_numpy() if len(tiles) else np.array([np.nan])
    summ = dict(method=method, n_tiles=int(len(tiles)), median_um=float(np.median(e)), p90_um=float(np.percentile(e, 90)),
                p99_um=float(np.percentile(e, 99)), max_um=float(np.max(e)))
    if method == "nuclei":
        reg = pd.DataFrame([dict(region=d["id"], **d["summary"]) for d in done])
        reg.to_csv(a.out / f"qc_{tag}_regions.csv", index=False)
        w = reg["n_matched"].to_numpy(float)
        summ.update(n_matched=int(w.sum()), n_nuclei=int(reg["n_nuclei"].sum()),
                    mean_dx_um=float(np.nansum(reg["median_dx_um"] * w) / w.sum()),
                    mean_dy_um=float(np.nansum(reg["median_dy_um"] * w) / w.sum()))
        summ["systematic_offset_um"] = float(np.hypot(summ["mean_dx_um"], summ["mean_dy_um"]))
    (a.out / f"qc_{tag}.json").write_text(json.dumps(dict(G=G.tolist(), summary=summ,
                                                          alignment=a.qc_alignment or "exported"), indent=1))
    what = "the exported affine" if tag == "final" else Path(a.qc_alignment).name
    log(f"qc ({method}) of {what}: {summ['n_tiles']} tiles, median {summ['median_um']:.2f} um, p90 {summ['p90_um']:.2f} um"
        + (f"; systematic offset ({summ['mean_dx_um']:+.2f}, {summ['mean_dy_um']:+.2f}) um from "
           f"{summ['n_matched']:,} nucleus pairs" if method == "nuclei" else ""))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lab = np.load(cache / "pieces_labels_coarse.npy")
    Md = np.array(g["M_dapi0_to_coarse"])
    ext = np.r_[apply(np.linalg.inv(Md), [[0, 0]])[0], apply(np.linalg.inv(Md), [[lab.shape[1], lab.shape[0]]])[0]] * px / 1000
    fig, ax = plt.subplots(figsize=(10, 10 * lab.shape[0] / lab.shape[1] + 1.2))
    ax.imshow(lab > 0, cmap="Greys", alpha=0.2, extent=[ext[0], ext[2], ext[3], ext[1]])
    if len(tiles):
        vmax = 1.5 if method == "nuclei" else max(1.0, np.percentile(e, 98))
        sc = ax.scatter(tiles.x_um / 1000, tiles.y_um / 1000, c=tiles.err_um, s=14 if method == "tiles" else 22,
                        cmap="viridis", vmin=0, vmax=vmax, marker="s")
        plt.colorbar(sc, ax=ax, shrink=0.6, label="local misalignment (um)")
        if method == "nuclei":
            ax.quiver(tiles.x_um / 1000, tiles.y_um / 1000, tiles.dx_um, -tiles.dy_um, color="k", alpha=0.6,
                      scale=40, width=0.0015)
    ax.set_xlim(ext[0], ext[2]); ax.set_ylim(ext[3], ext[1]); ax.set_aspect("equal")
    ax.set_xlabel("Xenium x (mm)"); ax.set_ylabel("Xenium y (mm)")
    how = ("median offset of matched nuclei per" if method == "nuclei" else "phase correlation of DAPI vs hematoxylin,")
    ax.set_title(f"{what}: {how} {a.tile_um:.0f} um square\n(median {summ['median_um']:.2f} um, "
                 f"p90 {summ['p90_um']:.2f} um; arrows = direction, exaggerated)", fontsize=10)
    fig.tight_layout()
    fig.savefig(a.out / ("04_residual_map_final.png" if tag == "final" else f"04_residual_map_{tag}.png"), dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------------------------
STAGES = ("global", "refine", "report", "qc", "all")


def run(stage, a):
    """Run one stage (or "all": global, refine, report, qc, then report again with the QC table)."""
    if stage not in STAGES:
        raise ValueError(f"stage must be one of {STAGES}, got {stage!r}")
    if a.out is None:
        raise HistoalignError("an output folder (out) is required")
    if a.threads:
        import SimpleITK as sitk
        sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(int(a.threads))
    a.out.mkdir(parents=True, exist_ok=True)
    if stage in ("global", "all"):
        stage_global(a)
    if stage in ("refine", "all"):
        stage_refine(a)
    if stage in ("report", "all"):
        stage_report(a)
    if stage in ("qc", "all"):
        stage_qc(a)
        if stage == "all":
            stage_report(a)
