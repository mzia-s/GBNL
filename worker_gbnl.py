#!/usr/bin/env python3
"""
main.py
Per-nucleotide worker for graded Betti computation from 1D nucleotide positions.

What this script is for
-----------------------
Given a single sequence and a chosen nucleotide (A, C, G, or T), it:
  1) builds a Vietoris–Rips complex in 1D from that nucleotide's positions,
  2) extracts plateau endpoints in the filtration,
  3) converts the maximal faces at each endpoint to a Stanley–Reisner ideal,
  4) calls Macaulay2 to compute graded Betti numbers, and
  5) writes a compact event log to a text file.

Notes for readers
-----------------
- Defaults are repository-relative. No user-specific directories.
- If Macaulay2 is not available on PATH, you can run it through a container:
    --use-container --m2-image ~/m2-1.24.05.sif
- If you want to run all four nucleotides for one sequence, see run_gbnl.sh.
"""


import argparse, os, re, shutil, subprocess, tempfile, time, sys
from datetime import datetime, timezone
import numpy as np

# ---------- Macaulay2 helpers ----------
M2_HELPERS = r'''
emitBetti = C -> (
    stdio << "BTABLE_BEGIN" << endl;
    stdio << toString net betti C << endl;
    stdio << "BTABLE_END" << endl;
    stdio << flush;
);
'''

def _run(cmd):
    return subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)

def _have_container_cli():
    return shutil.which("singularity") or shutil.which("apptainer")

def _check_image(path):
    return bool(path) and os.path.exists(path)

def _try_local_m2_probe(m2bin):
    """
    Return True if we can execute the local M2 binary; False if not.
    Specifically guard against OSError: [Errno 8] Exec format error.
    """
    if not m2bin:
        return False
    try:
        _ = subprocess.run([m2bin, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        return True
    except OSError as e:
        if getattr(e, "errno", None) == 8:
            return False
        raise  # other OS errors bubble up

def run_m2_program(m2_code: str, use_container: bool, image_path: str) -> str:
    """
    Execute a tiny M2 script and return stdout. Strategy:
      1) If use_container: run inside singularity/apptainer with the provided image.
      2) Else try local M2 (env var M2_BIN, then PATH search).
         If we hit Exec-format error, fallback to container automatically if available.
    """
    with tempfile.NamedTemporaryFile("w", suffix=".m2", delete=False) as f:
        f.write(m2_code)
        script_path = f.name
    try:
        # Resolve candidates
        container_cli = _have_container_cli()
        have_image = _check_image(image_path)
        local_m2 = os.environ.get("M2_BIN") or shutil.which("M2") or shutil.which("M2-binary")

        # Helper: run inside container
        def _run_in_container():
            if not container_cli:
                raise RuntimeError("Container runtime (singularity/apptainer) not found on PATH.")
            if not have_image:
                raise FileNotFoundError(f"M2 image not found: {image_path}")
            cli = container_cli
            out = ""
            # Many images expose either M2 or M2-binary
            for m2bin in ("M2", "M2-binary"):
                r = _run([cli, "exec", image_path, m2bin, "--script", script_path])
                out = r.stdout
                if "BTABLE_BEGIN" in out:
                    return out
            # If net betti not emitted, still return what we got (caller will parse and see empty)
            return out

        # Helper: run local binary
        def _run_local():
            if not local_m2:
                raise RuntimeError("Macaulay2 not found on PATH and M2_BIN not set.")
            # Try to execute; catch Exec format
            try:
                return _run([local_m2, "--script", script_path]).stdout
            except OSError as e:
                if getattr(e, "errno", None) == 8:
                    # Exec format error: wrong-arch binary; fall back to container if available
                    if container_cli and have_image:
                        return _run_in_container()
                    raise RuntimeError(
                        f"Local M2 at '{local_m2}' raised Exec format error (wrong architecture). "
                        f"Either load a proper Macaulay2 module or use --use-container with a valid SIF."
                    ) from e

        # Preferred path based on flag
        if use_container:
            return _run_in_container()
        else:
            # Use local if it probes OK; otherwise container if available
            if _try_local_m2_probe(local_m2):
                return _run_local()
            if container_cli and have_image:
                return _run_in_container()
            # Nothing workable
            raise RuntimeError(
                "No usable Macaulay2 backend found.\n"
                f"- Tried local M2: {local_m2 or 'None'} (not runnable)\n"
                f"- Container CLI: {container_cli or 'None'}\n"
                f"- Image exists: {have_image} ({image_path})"
            )
    finally:
        try:
            os.remove(script_path)
        except Exception:
            pass

def m2_code_for_sr_from_facets(n_vertices: int, facets) -> str:
    """Maximal faces as monomials → SR ideal → minimal free resolution."""
    def mono(F): return "*".join(f"a_{int(v)}" for v in sorted(F))
    facet_list = ", ".join(mono(F) for F in facets) if facets else ""
    return M2_HELPERS + rf'''
needsPackage "SimplicialComplexes";
S = QQ[a_0..a_{n_vertices-1}];
D = simplicialComplex {{ {facet_list} }};
J = ideal D;                 -- Stanley–Reisner ideal I_Δ
C = res (S^1 / J);           -- minimal free resolution of S/J
emitBetti C;
exit 0
'''

def parse_betti_stdout(text: str):
    """Parse 'net betti' emitted between BTABLE_BEGIN/END to [(i,j,beta), ...]."""
    lines = text.splitlines()
    try:
        b = lines.index("BTABLE_BEGIN"); e = lines.index("BTABLE_END")
    except ValueError:
        return []
    header_nums = [int(x) for x in re.findall(r"\d+", lines[b+1])]
    triples = []
    for ln in lines[b+3:e]:
        m = re.match(r"\s*(\d+):\s*(.*)$", ln)
        if not m: continue
        r = int(m.group(1)); toks = m.group(2).split()
        for idx, tok in enumerate(toks):
            if tok == ".": continue
            beta = int(tok); i = header_nums[idx]; j = i + r
            triples.append((i, j, beta))
    return sorted(triples, key=lambda t: (t[0], t[1]))

def format_betti(triples):
    return "BETTI: -" if not triples else "BETTI: " + ", ".join(f"beta({i},{j})={b}" for (i,j,b) in triples)

# ---------- Gudhi VR helpers (from 1D points) ----------
def build_vr_simplex_tree(points: np.ndarray, max_e: float, vr_dim: int):
    try:
        import gudhi
    except Exception as e:
        raise RuntimeError("Gudhi is required (pip install gudhi).") from e
    rc = gudhi.RipsComplex(points=points, max_edge_length=max_e)
    return rc.create_simplex_tree(max_dimension=vr_dim)

def one_skeleton_counts(st, eps: float):
    V = set(); E = set(); T = set()
    for sig, f in st.get_simplices():
        if f > eps: continue
        k = len(sig)
        if k == 1: V.add(sig[0])
        elif k == 2: E.add(tuple(sorted(sig)))
        elif k == 3: T.add(tuple(sorted(sig)))
    return len(V), len(E), len(T)

def facets_at_eps(st, eps: float, vr_dim: int):
    """Return maximal simplices present at filtration ≤ eps (up to vr_dim)."""
    simplices = [(tuple(sig), filt) for sig, filt in st.get_simplices()
                 if filt <= eps and (len(sig)-1) <= vr_dim]
    # Prefer larger faces first so supersets are checked before subsets
    simplices.sort(key=lambda x: len(x[0]), reverse=True)
    facets, supersets = [], []
    for sig, _ in simplices:
        S = set(sig)
        if not any(S.issubset(T) for T in supersets):
            facets.append(sig); supersets.append(S)
    if not facets:
        verts = [tuple(sig) for sig, f in st.get_simplices() if len(sig)==1 and f <= eps]
        facets = verts
    return facets

def canonical_facet_key(facets):
    """Hashable canonical key for caching."""
    return tuple(sorted(tuple(sorted(f)) for f in facets))

def event_epsilons(st, min_e: float, max_e: float, vr_dim: int):
    """Unique filtration values where edges or higher appear within [min_e, max_e]."""
    eps = set()
    has_vertex = False
    for sig, f in st.get_simplices():
        k = len(sig) - 1
        if k == 0:
            has_vertex = True
            continue
        if k <= vr_dim and (min_e - 1e-12) <= f <= (max_e + 1e-12):
            eps.add(float(f))
    ev = sorted(eps)
    if not ev:
        fallback = min(max_e, max(min_e, 0.0 if has_vertex else min_e))
        ev = [fallback]
    return ev

# ---------- nucleotide -> 1D points ----------
def nuc_points_from_raw(raw_seq: str, nuc: str):
    """Positions are 1-based indices of the chosen nucleotide (U→T); points are those indices in R^1."""
    s = raw_seq.upper().replace("U","T")
    idx = [i+1 for i, c in enumerate(s) if c == nuc]
    pts = np.asarray(idx, dtype=float).reshape(-1, 1)
    return s, idx, pts

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seq-id", required=True)
    ap.add_argument("--seq-ord", type=int, default=None)
    ap.add_argument("--sequence", required=True)
    ap.add_argument("--nuc", required=True, choices=list("ACGT"))
    ap.add_argument("--out-root", default="./outputs")
    ap.add_argument("--vr-dim", type=int, default=1, choices=[1,2])
    ap.add_argument("--min-e", type=float, default=0.0)
    ap.add_argument("--max-e", type=float, default=50.0)
    ap.add_argument("--m2-image", default=os.path.expanduser("~/m2-1.24.05.sif"))
    ap.add_argument("--use-container", action="store_true")
    args = ap.parse_args()

    # Normalize and extract 1D points
    seq_used, pos_1based, pts = nuc_points_from_raw(args.sequence, args.nuc)

    # Output paths
    prefix = f"{args.seq_ord:03d}_" if args.seq_ord is not None else ""
    out_dir = os.path.join(args.out_root, f"{prefix}{args.seq_id}")
    os.makedirs(out_dir, exist_ok=True)
    filename = f"{prefix}{args.seq_id}_{args.nuc}_vr{args.vr_dim}_k1.txt"
    out_path = os.path.join(out_dir, filename)

    # Absent nucleotide stub
    if pts.size == 0:
        with open(out_path, "w") as out:
            out.write(f"no {args.nuc} is present IN {seq_used}\n")
        print(f"[{args.seq_id}:{args.nuc}] Wrote absence stub: {filename}")
        return

    # Worker timing
    worker_start_iso = datetime.now(timezone.utc).isoformat()
    t0 = time.time()

    # Build Rips from points, collect event epsilons, build plateau intervals
    st = build_vr_simplex_tree(pts, max_e=args.max_e, vr_dim=args.vr_dim)
    ev_eps = [e for e in event_epsilons(st, args.min_e, args.max_e, args.vr_dim)
              if args.min_e - 1e-12 <= e <= args.max_e + 1e-12]
    ev_eps = sorted(set(ev_eps))

    intervals = []
    left = float(args.min_e)
    for curr in ev_eps:
        if curr < left + 1e-12:
            left = curr
            continue
        if curr > args.max_e + 1e-12:
            break
        intervals.append((left, curr))
        left = curr
    if left < args.max_e - 1e-12:
        intervals.append((left, float(args.max_e)))

    # M2 cache keyed by canonical facets, last-key to suppress duplicates
    betti_cache = {}
    last_key = None

    with open(out_path, "w") as out:
        # ---- Header (kept identical to no-M2 worker) ----
        print(f"# SEQ: {args.seq_id}", file=out)
        if args.seq_ord is not None:
            print(f"# ORDER_INDEX: {args.seq_ord:03d}", file=out)
        print(f"# WORKER_START_ISO_UTC: {worker_start_iso}", file=out)
        print(f"# COMPLEX: VR(dim={args.vr_dim})", file=out)
        print("# INPUT_MODE: distance_matrix(|pos_i - pos_j|)", file=out)  # keep same string for parity
        if args.vr_dim >= 2:
            print("# NOTATION: VTX=vertices, EDG=edges, TRI=filled triangles (2-simplices)", file=out)
        else:
            print("# NOTATION: VTX=vertices, EDG=edges (graph-only; no triangles in VR(dim=1))", file=out)
        print("# KMER: 1", file=out)
        print(f"# NUCLEOTIDE: {args.nuc}", file=out)
        print(f"# SEQUENCE_USED (upper, U→T): {seq_used}  (len={len(seq_used)})", file=out)
        print(f"# POSITIONS(1-based) of '{args.nuc}': {pos_1based}", file=out)
        print(f"# EPS_EVENT_COUNT: {len(ev_eps)}", file=out)
        print(f"# EPS_EVENTS: {ev_eps}", file=out)
        print("# PLATEAU_INTERVALS: (left, right], evaluated at 'right' only; plus ε=min_e snapshot.", file=out)
        print("# BETTI MAP: M2 table cell (col=i, row=r) is beta(i,i+r); we print as beta(i,j) with j=i+r.", file=out)
        print(f"[NUC={args.nuc}]", file=out)

        # ---------- ε = min_e snapshot ----------
        eps0 = float(args.min_e)
        facets0 = facets_at_eps(st, eps0, vr_dim=args.vr_dim)
        key0 = canonical_facet_key(facets0)
        vtx0, edg0, tri0 = one_skeleton_counts(st, eps0)
        triples0 = betti_cache.get(key0)
        if triples0 is None:
            m2_code0 = m2_code_for_sr_from_facets(n_vertices=pts.shape[0], facets=facets0)
            stdout0 = run_m2_program(m2_code0, use_container=args.use_container, image_path=args.m2_image)
            triples0 = parse_betti_stdout(stdout0)
            betti_cache[key0] = triples0

        if args.vr_dim >= 2:
            print(f"INT=[{eps0:.6g},{eps0:.6g}]  EPS_AT_RIGHT={eps0:.6g}  | VTX={vtx0} | EDG={edg0} | TRI={tri0}  | {format_betti(triples0)}", file=out)
        else:
            print(f"INT=[{eps0:.6g},{eps0:.6g}]  EPS_AT_RIGHT={eps0:.6g}  | VTX={vtx0} | EDG={edg0}  | {format_betti(triples0)}", file=out)
        last_key = key0
        # ---------- end ε = min_e snapshot ----------

        # ---------- Plateau endpoints only (print if Betti changed) ----------
        for (a, b) in intervals:
            eps = float(b)
            facets = facets_at_eps(st, eps, vr_dim=args.vr_dim)
            key = canonical_facet_key(facets)
            vtx, edg, tri = one_skeleton_counts(st, eps)

            if key == last_key:
                continue
            last_key = key

            triples = betti_cache.get(key)
            if triples is None:
                m2_code = m2_code_for_sr_from_facets(n_vertices=pts.shape[0], facets=facets)
                stdout = run_m2_program(m2_code, use_container=args.use_container, image_path=args.m2_image)
                triples = parse_betti_stdout(stdout)
                betti_cache[key] = triples

            if args.vr_dim >= 2:
                line = (
                    f"INT=({a:.6g},{b:.6g}]  EPS_AT_RIGHT={eps:.6g}  "
                    f"| VTX={vtx} | EDG={edg} | TRI={tri}  | {format_betti(triples)}"
                )
            else:
                line = (
                    f"INT=({a:.6g},{b:.6g}]  EPS_AT_RIGHT={eps:.6g}  "
                    f"| VTX={vtx} | EDG={edg}  | {format_betti(triples)}"
                )
            print(line, file=out)

        worker_end_iso = datetime.now(timezone.utc).isoformat()
        runtime_sec = time.time() - t0
        print(f"# WORKER_END_ISO_UTC: {worker_end_iso}", file=out)
        print(f"# WORKER_RUNTIME_SEC: {runtime_sec:.3f}", file=out)

    print(f"[{args.seq_id}:{args.nuc}] Wrote {filename} in {out_dir}")

if __name__ == "__main__":
    main()
