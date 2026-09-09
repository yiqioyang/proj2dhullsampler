"""Visualize how history-matching constraints carve up the parameter space.

Two complementary views, both built on top of an already-prepared
``HistoryMatching`` case (i.e. after ``prepare_for_sampling`` / ``orchestrate``
have run, so that ``hm.paras_vars`` holds the final diagnostic groups and
``hm.results.para_l`` holds the order in which the parameter pairs were
applied):

1. :func:`animate_pair_constraints` - for one parameter pair, an animation that
   starts with the emulated samples filling the whole (p1, p2) square and then
   switches on that pair's diagnostics one at a time, updating the surviving
   points and growing a list of the applied diagnostics down the side of the
   figure.

   In notebook mode the animation is displayed inline with play/pause/scrub
   controls for a single (specified) pair; in python mode an animation file is
   written for every parameter pair.

2. :func:`plot_constraint_cascade` - a matrix of scatter panels saved as a PDF.
   Columns are parameter pairs; row 0 shows each pair constrained only by its
   own diagnostics, row 1 shows each pair after the first pair's constraint is
   also imposed, row 2 after the first two, and so on. Reading down a column
   shows how much of that pair's own feasible region is eaten away by
   constraints living in *other* parameter subspaces - i.e. how interlocked the
   parameters are.

Nothing in here mutates the case: every function only reads ``hm``.

Typical use::

    from proj2dhullsampler.history_matching_animation import (
        list_pairs, animate_pair_constraints, plot_constraint_cascade,
    )

    list_pairs(hm)                             # index -> pair -> n diagnostics
    ani = animate_pair_constraints(hm, pair=0) # notebook: inline player
    plot_constraint_cascade(hm)                # -> <case>/diagnostics/*.pdf
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

# ---------------------------------------------------------------------------
# Appearance
# ---------------------------------------------------------------------------

# The three point layers are disjoint sets, not stacked translucent overlays:
# each is thinned to `max_points` on its own, so a category stays visible even
# when it holds only a sliver of the samples.
COLOR_ALL = "#e0e0e0"       # already ruled out before this step (context field)
COLOR_CUT = "#fc9272"       # ruled out *by* this step
COLOR_NOW = "#08306b"       # still surviving
COLOR_HULL = "#000000"      # alpha-shape outline (dashed)
COLOR_PENDING = "#b0b0b0"   # not-yet-applied entries in the side list
COLOR_APPLIED = "#238b45"   # spine colour for pairs already folded into the cascade


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def _in_notebook():
    """True if we appear to be running inside a Jupyter/IPython kernel."""
    try:
        from IPython import get_ipython
    except ImportError:
        return False

    ip = get_ipython()
    if ip is None:
        return False
    # ZMQInteractiveShell -> notebook/lab/qtconsole;
    # TerminalInteractiveShell -> plain ipython REPL
    return type(ip).__name__ == "ZMQInteractiveShell"


def _resolve_mode(hm, mode=None):
    """Decide between 'notebook' and 'python'.

    Explicit `mode` wins; otherwise the case's own `hm.mode` is used (so this
    module follows the same switch as `visualize_check` / `compare_with_original`);
    otherwise we fall back to sniffing the runtime. If the case's mode disagrees
    with the detected runtime we say so rather than silently overriding it.
    """
    if mode is not None:
        if mode not in ("notebook", "python"):
            raise ValueError(f"mode must be 'notebook' or 'python', got {mode!r}")
        return mode

    case_mode = getattr(hm, "mode", None)
    detected = "notebook" if _in_notebook() else "python"

    if case_mode in ("notebook", "python"):
        if case_mode != detected:
            print(
                f"Note: the case was built with mode={case_mode!r} but this looks "
                f"like a {detected!r} session. Using {case_mode!r}; pass "
                f"mode={detected!r} to override."
            )
        return case_mode

    return detected


def _require_prepared(hm):
    """Fail early (and readably) if the case has not been taken far enough."""
    if not hasattr(hm, "tf_masks"):
        raise AttributeError(
            "hm has no tf_masks; call hm.load_case() and hm.load_mask(...) first."
        )
    if not hasattr(hm, "paras_vars"):
        raise AttributeError(
            "hm has no paras_vars; the diagnostics have not been grouped by "
            "parameter pair yet. Run remove_var2d_auto()/group_para_climatology() "
            "and prepare_for_sampling() first."
        )
    if not hasattr(hm, "p_emu"):
        raise AttributeError("hm has no p_emu; call hm.load_case() first.")


def pair_sequence(hm):
    """The parameter pairs in the order the sampler applies them.

    Prefers `hm.results.para_l` (the order `orchestrate` actually walked, minus
    any pair it had to skip) and falls back to the `hm.paras_vars` ordering,
    which is how the pairs are queued up before orchestration.
    """
    _require_prepared(hm)

    para_l = getattr(getattr(hm, "results", None), "para_l", None)
    if para_l:
        pairs = [tuple(p) for p in para_l if tuple(p) in hm.paras_vars]
        if pairs:
            return pairs
        warnings.warn(
            "hm.results.para_l does not line up with hm.paras_vars; falling back "
            "to the hm.paras_vars ordering.",
            stacklevel=2,
        )

    return [tuple(p) for p in hm.paras_vars]


def _pair_vars(hm, pair, var_order="given"):
    """Diagnostics attached to `pair`, optionally reordered.

    `var_order`:
        "given"            - the order stored in hm.paras_vars (default)
        "restrictive_last" - loosest diagnostic first, tightest last
        "restrictive_first"- tightest diagnostic first
    """
    vars_ = list(hm.paras_vars[pair])

    missing = [v for v in vars_ if v not in hm.tf_masks.columns]
    if missing:
        warnings.warn(
            f"{pair}: diagnostics {missing} are no longer columns of tf_masks; "
            "skipping them.",
            stacklevel=2,
        )
        vars_ = [v for v in vars_ if v not in missing]

    if var_order == "given":
        return vars_
    if var_order not in ("restrictive_last", "restrictive_first"):
        raise ValueError(
            "var_order must be 'given', 'restrictive_last' or 'restrictive_first', "
            f"got {var_order!r}"
        )

    counts = {v: int(hm.tf_masks[v].to_numpy().sum()) for v in vars_}
    return sorted(vars_, key=counts.get, reverse=(var_order == "restrictive_last"))


def resolve_pair(hm, pair):
    """Turn a user-supplied pair reference into a key of `hm.paras_vars`.

    Accepts an integer index into :func:`pair_sequence`, a (p1, p2) tuple/list
    in either order, or the "(p1, p2)" string form used in the saved
    specifications JSON.
    """
    pairs = pair_sequence(hm)

    if isinstance(pair, (int, np.integer)):
        return pairs[int(pair)]

    if isinstance(pair, str):
        parts = [p.strip() for p in pair.strip("() ").split(",")]
        pair = tuple(parts)

    pair = tuple(pair)
    if pair in hm.paras_vars:
        return pair
    if pair[::-1] in hm.paras_vars:
        return pair[::-1]

    raise KeyError(
        f"{pair} is not one of the case's parameter pairs. "
        f"Available pairs:\n" + "\n".join(f"  {i}: {p}" for i, p in enumerate(pairs))
    )


def list_pairs(hm, verbose=True):
    """Print (and return) the parameter pairs with their diagnostic counts."""
    pairs = pair_sequence(hm)
    rows = [(i, p, len(hm.paras_vars[p])) for i, p in enumerate(pairs)]

    if verbose:
        print(f"{len(rows)} parameter pairs (in the order the sampler applies them):")
        for i, p, n in rows:
            print(f"  [{i:>2}] {p[0]:<28} x {p[1]:<28} {n:>3} diagnostics")

    return rows


def _short(name, width=18):
    """Shorten a long parameter/diagnostic name for use in a title."""
    if len(name) <= width:
        return name
    return name[: width - 1] + "…"


def _xy_and_masks(hm, pair, vars_):
    """Point coordinates and per-diagnostic boolean masks, in one shuffled order.

    Everything is returned as positional numpy arrays (never index-aligned
    pandas objects) and pre-permuted with a fixed shuffle, so that "take the
    first `max_points` survivors" is a stable, reproducible thinning that does
    not make points jump around between frames.
    """
    xy = hm.p_emu[list(pair)].to_numpy(dtype=float)
    if vars_:
        masks = hm.tf_masks[vars_].to_numpy(dtype=bool)
    else:
        masks = np.empty((len(xy), 0), dtype=bool)

    if masks.shape[0] != xy.shape[0]:
        raise ValueError(
            f"p_emu has {xy.shape[0]} rows but tf_masks has {masks.shape[0]}; "
            "they must describe the same emulated samples."
        )
    return xy, masks


def _thin(xy, keep, max_points):
    """First `max_points` rows of `xy` where `keep` is True (both pre-shuffled)."""
    sel = xy[keep]
    if max_points is not None and sel.shape[0] > max_points:
        sel = sel[:max_points]
    return sel if sel.size else np.empty((0, 2))


def _plot_hull(ax, geom, **kwargs):
    """Outline a shapely Polygon/MultiPolygon; silently skip degenerate shapes."""
    if geom is None or getattr(geom, "is_empty", True):
        return

    geoms = list(getattr(geom, "geoms", [geom]))
    for g in geoms:
        exterior = getattr(g, "exterior", None)
        if exterior is None:
            continue
        x, y = exterior.xy
        ax.plot(x, y, **kwargs)


def _hull_for(hm, pair):
    """The alpha-shape for `pair`, preferring the post-orchestrate version."""
    valid = getattr(getattr(hm, "results", None), "valid_hulls", None)
    if valid and pair in valid:
        return valid[pair]
    return getattr(hm, "grouped_hulls", {}).get(pair)


def _style_axes(ax, pair, labels=True, fontsize=8):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([0, 0.5, 1])
    ax.tick_params(labelsize=fontsize - 2)
    if labels:
        ax.set_xlabel(pair[0], fontsize=fontsize)
        ax.set_ylabel(pair[1], fontsize=fontsize)
    else:
        ax.set_xticklabels([])
        ax.set_yticklabels([])


# ---------------------------------------------------------------------------
# 1. Per-pair animation
# ---------------------------------------------------------------------------


def _legend(fig_or_ax, cut_label, keep_label, fontsize=7, **kwargs):
    """Proxy-artist legend describing the three point categories."""
    from matplotlib.lines import Line2D

    handles = [
        Line2D([], [], ls="", marker="o", ms=4, mfc=COLOR_ALL, mec="none",
               label="already ruled out"),
        Line2D([], [], ls="", marker="o", ms=4, mfc=COLOR_CUT, mec="none",
               label=cut_label),
        Line2D([], [], ls="", marker="o", ms=4, mfc=COLOR_NOW, mec="none",
               label=keep_label),
        Line2D([], [], ls="--", color=COLOR_HULL, lw=1.2, label="alpha-shape hull"),
    ]
    return fig_or_ax.legend(handles=handles, fontsize=fontsize, **kwargs)


def build_pair_animation(
    hm,
    pair,
    max_points=3000,
    var_order="given",
    show_hull=True,
    interval=900,
    hold_last=3,
    random_state=0,
    figsize=(11, 5.5),
    point_size=6,
):
    """Build (but do not display or save) the constraint animation for one pair.

    Frame 0 shows every emulated sample filling the (p1, p2) square; frame k
    shows the samples surviving the first k diagnostics of the pair, with the
    diagnostics listed down the right-hand side as they are switched on. Points
    ruled out by the diagnostic added at frame k are highlighted so each step
    shows *what it cost*. Returns ``(fig, FuncAnimation)``.

    `hold_last` repeats the final frame that many extra times so the fully
    constrained region stays on screen for a beat before the loop restarts.
    """
    _require_prepared(hm)
    pair = resolve_pair(hm, pair)
    vars_ = _pair_vars(hm, pair, var_order=var_order)

    rng = np.random.default_rng(random_state)
    xy, masks = _xy_and_masks(hm, pair, vars_)
    perm = rng.permutation(xy.shape[0])
    xy, masks = xy[perm], masks[perm]
    n_total = xy.shape[0]

    # Cumulative survivor mask after each diagnostic; frame 0 == everything.
    keeps = [np.ones(n_total, dtype=bool)]
    for k in range(len(vars_)):
        keeps.append(keeps[-1] & masks[:, k])

    # Per frame: what is still alive, and what this very step killed.
    keep_pts = [_thin(xy, k_, max_points) for k_ in keeps]
    cut_pts = [np.empty((0, 2))] + [
        _thin(xy, keeps[k - 1] & ~keeps[k], max_points)
        for k in range(1, len(keeps))
    ]
    counts = [int(k_.sum()) for k_ in keeps]

    # A fixed thinning of the full square, drawn once as the grey context field
    # so the plot stays as "full" as it started and nothing jitters between frames.
    background = _thin(xy, np.ones(n_total, bool), max_points)

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 2, width_ratios=[3.0, 2.0], wspace=0.08)
    ax = fig.add_subplot(gs[0, 0])
    ax_txt = fig.add_subplot(gs[0, 1])
    ax_txt.axis("off")

    _style_axes(ax, pair, labels=True, fontsize=10)
    ax.scatter(
        background[:, 0], background[:, 1],
        s=point_size, c=COLOR_ALL, linewidths=0, zorder=1,
    )
    sc_cut = ax.scatter([], [], s=point_size, c=COLOR_CUT, linewidths=0, zorder=2)
    sc_now = ax.scatter([], [], s=point_size, c=COLOR_NOW, linewidths=0, zorder=3)

    hull_artists = []
    if show_hull:
        # Only drawn on the final frame: the hull is fitted to the *combined*
        # mask, so it is only meaningful once every diagnostic is switched on.
        before = list(ax.lines)
        _plot_hull(ax, _hull_for(hm, pair), color=COLOR_HULL, lw=1.3, ls="--", zorder=4)
        hull_artists = [ln for ln in ax.lines if ln not in before]
        for ln in hull_artists:
            ln.set_alpha(0.0)

    title = ax.set_title("", fontsize=11, loc="left")

    # Side panel: one Text artist per diagnostic, recoloured as frames advance.
    fontsize = float(np.clip(150.0 / max(len(vars_), 1), 5.5, 9.5))
    ax_txt.text(
        0.0, 1.0,
        f"{pair[0]}\n  × {pair[1]}\n{len(vars_)} diagnostics",
        transform=ax_txt.transAxes, va="top", ha="left",
        fontsize=fontsize + 1.5, fontweight="bold",
    )
    top, bottom = 0.86, 0.14
    step = min(0.036, (top - bottom) / max(len(vars_), 1))
    var_texts = []
    for i, v in enumerate(vars_):
        var_texts.append(
            ax_txt.text(
                0.0, top - i * step, f"{i + 1:>2}. {v}",
                transform=ax_txt.transAxes, va="top", ha="left",
                fontsize=fontsize, color=COLOR_PENDING, family="monospace",
            )
        )
    _legend(
        ax_txt, "ruled out by this step", "still surviving",
        fontsize=7, loc="lower left", frameon=False,
        bbox_to_anchor=(0.0, -0.02),
    )

    def update(frame):
        k = min(frame, len(vars_))  # frames past the end are the "hold" repeats

        sc_cut.set_offsets(cut_pts[k])
        sc_now.set_offsets(keep_pts[k])

        pct = 100.0 * counts[k] / n_total
        if k == 0:
            head = "step 0 - no constraint yet"
        else:
            lost = counts[k - 1] - counts[k]
            head = f"step {k}/{len(vars_)} - added {vars_[k - 1]}  (-{lost:,})"
        title.set_text(
            f"{head}\nsurviving: {counts[k]:,} / {n_total:,}  ({pct:.2f}%)"
        )

        for i, t in enumerate(var_texts):
            if i < k - 1:
                t.set_color("black")
                t.set_fontweight("normal")
            elif i == k - 1:
                t.set_color(COLOR_NOW)
                t.set_fontweight("bold")
            else:
                t.set_color(COLOR_PENDING)
                t.set_fontweight("normal")

        for ln in hull_artists:
            ln.set_alpha(0.9 if k == len(vars_) else 0.0)

        return [sc_cut, sc_now, title] + var_texts + list(hull_artists)

    frames = list(range(len(vars_) + 1)) + [len(vars_)] * max(hold_last, 0)
    ani = animation.FuncAnimation(
        fig, update, frames=frames, interval=interval, blit=False, repeat=True,
    )
    return fig, ani


def animate_pair_constraints(
    hm,
    pair=None,
    mode=None,
    outdir=None,
    fmt="gif",
    fps=1.2,
    dpi=110,
    display_inline=True,
    **kwargs,
):
    """Animate the per-diagnostic shrinking of one (or every) parameter pair.

    Notebook mode: builds the animation for a single `pair` (defaulting to the
    first one, with a note saying so) and displays it inline as an HTML/JS
    player with play, pause, loop and frame-step controls. Returns the
    ``FuncAnimation``; keep a reference to it or it will be garbage collected.

    Python mode: writes one animation file per parameter pair into `outdir`
    (default ``<case>/diagnostics/animations/``) and returns the list of paths.
    Passing an explicit `pair` in python mode restricts it to that pair.

    Extra keyword arguments are forwarded to :func:`build_pair_animation`
    (`max_points`, `var_order`, `show_hull`, `interval`, `figsize`, ...).
    """
    _require_prepared(hm)
    mode = _resolve_mode(hm, mode)
    pairs = pair_sequence(hm)

    if not pairs:
        raise ValueError("No parameter pairs left to animate.")

    if mode == "notebook":
        if pair is None:
            pair = pairs[0]
            print(
                f"No pair given; animating {pair}. "
                "Call list_pairs(hm) to see the others, then pass e.g. pair=3."
            )
        pair = resolve_pair(hm, pair)
        fig, ani = build_pair_animation(hm, pair, **kwargs)

        if display_inline:
            try:
                from IPython.display import HTML, display

                html = ani.to_jshtml(fps=fps)
                plt.close(fig)  # otherwise the last frame also renders as a still
                display(HTML(html))
            except ImportError:
                warnings.warn(
                    "IPython is not available; returning the animation without "
                    "displaying it.",
                    stacklevel=2,
                )
        return ani

    # python mode -> write files
    if outdir is None:
        outdir = Path(hm.diagnostics_dir) / "animations"
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    todo = [resolve_pair(hm, pair)] if pair is not None else pairs

    paths = []
    for i, p in enumerate(todo):
        fig, ani = build_pair_animation(hm, p, **kwargs)
        path = outdir / f"pair_{i:02d}_{p[0]}__{p[1]}.{fmt}"
        writer = "pillow" if fmt == "gif" else "ffmpeg"
        ani.save(str(path), writer=writer, fps=fps, dpi=dpi)
        plt.close(fig)
        paths.append(path)
        print(f"[{i + 1}/{len(todo)}] wrote {path}")

    return paths


# ---------------------------------------------------------------------------
# 2. Cascade matrix
# ---------------------------------------------------------------------------


def plot_constraint_cascade(
    hm,
    pairs=None,
    n_rows=None,
    max_points=800,
    panel_size=2.0,
    point_size=2.0,
    show_hull=True,
    random_state=0,
    outpath=None,
    save=True,
    dpi=150,
):
    """Matrix of scatter panels showing how the parameter pairs interlock.

    Columns are parameter pairs (in sampling order), each panel plotting the
    pair's own 2-D subspace with the *first* parameter on x and the second on y.
    Row 0 shows each pair constrained by its own diagnostics only. Row *i*
    (>= 1) shows each pair constrained by its own diagnostics **and** by
    everything belonging to the first *i* pairs, so reading down a column shows
    how much of that pair's feasible region is eaten away by constraints that
    live in other subspaces.

    Each panel splits the samples into three disjoint sets: outside the pair's
    own feasible region (grey), inside its own region but killed by the other
    pairs in the cascade (salmon), and surviving everything (dark blue). The
    three sets are thinned to `max_points` *independently*, so that a category
    holding only a sliver of the samples still shows up; the trade-off is that
    point density is not comparable between colours (the per-panel counts are).

    A column gets a green frame once its own pair has entered the cascade. For
    those panels the dark-blue set is exactly the cumulative surviving set, so
    the count printed in the corner is the same across every green-framed column
    in that row. Their salmon set is generally *not* empty - it is the points
    this pair accepts but some other applied pair rejects.

    With every pair as a column the grid is n_pairs x (n_pairs + 1) panels, so
    the PDF is physically large; `max_points` keeps it to a sane file size (at
    the default, ~11 MB for 18 pairs). Use `pairs` and/or `n_rows` to cut it
    down. The scatters are deliberately *not* rasterized: the PDF backend routes
    rasterized artists through the mixed-mode renderer, which allocates a
    full-figure-sized buffer per artist and blows up to many GB on a grid this
    size.

    Returns the figure; also writes a PDF unless `save=False`.
    """
    _require_prepared(hm)

    all_pairs = pair_sequence(hm)
    if pairs is None:
        pairs = all_pairs
    else:
        pairs = [resolve_pair(hm, p) for p in pairs]
    n_pairs = len(pairs)
    if n_pairs == 0:
        raise ValueError("No parameter pairs to plot.")

    # row 0 = own only; rows 1..n_pairs = cumulative through pairs[0..i-1]
    if n_rows is None:
        n_rows = len(all_pairs) + 1
    n_rows = int(min(n_rows, len(all_pairs) + 1))

    rng = np.random.default_rng(random_state)
    n_total = hm.p_emu.shape[0]
    perm = rng.permutation(n_total)

    # Per-pair coordinates and own-constraint masks, all in the same shuffled order.
    xy_by_pair, own_by_pair = {}, {}
    for p in pairs:
        xy, masks = _xy_and_masks(hm, p, _pair_vars(hm, p))
        xy_by_pair[p] = xy[perm]
        own_by_pair[p] = (
            masks[perm].all(axis=1) if masks.shape[1] else np.ones(n_total, bool)
        )

    # The cascade follows the *sampling* order, which may include pairs the
    # caller left out of `pairs`; those still constrain, they just get no column.
    cascade_pairs = all_pairs[: n_rows - 1]
    cumulative = [np.ones(n_total, dtype=bool)]
    for p in cascade_pairs:
        _, masks = _xy_and_masks(hm, p, _pair_vars(hm, p))
        step = masks[perm].all(axis=1) if masks.shape[1] else np.ones(n_total, bool)
        cumulative.append(cumulative[-1] & step)

    fig, axes = plt.subplots(
        n_rows, n_pairs,
        figsize=(panel_size * n_pairs + 2.6, panel_size * n_rows + 1.8),
        squeeze=False,
    )

    for i in range(n_rows):
        cum = cumulative[i]
        applied = set(cascade_pairs[:i])

        for j, p in enumerate(pairs):
            ax = axes[i][j]
            xy = xy_by_pair[p]
            own = own_by_pair[p]
            joint = cum & own

            _style_axes(ax, p, labels=False)
            if p in applied:
                for spine in ax.spines.values():
                    spine.set_edgecolor(COLOR_APPLIED)
                    spine.set_linewidth(1.8)

            for keep, color, z in (
                (~own, COLOR_ALL, 1),
                (own & ~joint, COLOR_CUT, 2),
                (joint, COLOR_NOW, 3),
            ):
                pts = _thin(xy, keep, max_points)
                ax.scatter(
                    pts[:, 0], pts[:, 1],
                    s=point_size, c=color, linewidths=0, zorder=z,
                )
            if show_hull:
                _plot_hull(ax, _hull_for(hm, p), color=COLOR_HULL, lw=0.8,
                           ls="--", zorder=4)

            n_own = int(own.sum())
            n_joint = int(joint.sum())
            frac = 100.0 * n_joint / n_own if n_own else 0.0
            ax.text(
                0.03, 0.97, f"{n_joint:,}\n{frac:.0f}% of own",
                transform=ax.transAxes, va="top", ha="left", fontsize=7,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.8),
            )

            if i == 0:
                ax.set_title(
                    f"x: {_short(p[0], 20)}\ny: {_short(p[1], 20)}",
                    fontsize=8, pad=5,
                )

        label = (
            "row 0\nown constraint\nonly"
            if i == 0
            else f"row {i}\n+ {_short(cascade_pairs[i - 1][0], 18)}\n"
            f"   × {_short(cascade_pairs[i - 1][1], 18)}"
        )
        axes[i][0].set_ylabel(label, fontsize=8, rotation=0, ha="right", va="center")

    # Reserve a constant-height header band for the title + legend, expressed as
    # a fraction of this figure's height (which grows with the number of rows).
    head = min(1.6 / (panel_size * n_rows + 1.8), 0.4)
    fig.suptitle(
        f"Constraint cascade - {getattr(hm, 'case_name', '')}\n"
        "columns: parameter pairs (sampling order)   |   "
        "rows: cumulative constraints from the pairs listed on the left\n"
        "axes are 0-1 normalized parameters, x = first, y = second",
        fontsize=11, y=1.0 - 0.10 * head, va="top",
    )
    _legend(
        fig, "killed by the other pairs", "surviving everything",
        fontsize=8, loc="upper center", ncol=4, frameon=False,
        bbox_to_anchor=(0.5, 1.0 - 0.72 * head),
    )
    fig.tight_layout(rect=[0, 0, 1, 1.0 - head])

    if save:
        if outpath is None:
            outpath = Path(hm.diagnostics_dir) / "constraint_cascade.pdf"
        outpath = Path(outpath)
        outpath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
        print(f"wrote {outpath}")

    return fig
