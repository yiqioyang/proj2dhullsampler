"""Visualize how history-matching constraints carve up the parameter space.

Three views, all built on top of an already-prepared
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

2. :func:`plot_constraint_interlock` - a two-column matrix of scatter panels
   saved as a PDF, one row per parameter pair. On the left, ``hm.p_emu`` split
   by that pair's own diagnostics alone - the whole region the pair itself
   allows. On the right, the identical square with those survivors as a
   backdrop and the parameter sets the sampler actually *drew* scattered on
   top. Where the drawn samples leave the pair's own region uncovered is what
   the interlocking with the *other* pairs costs it - i.e. how much of a pair's
   feasible region is really eaten away by constraints living in other
   parameter subspaces.

   This one needs the case taken one step further than the animation: its right
   column is ``hm.results.unscaled_samples``, so it has to run after
   ``hm.draw()``.

3. :func:`plot_dropped_constraints` - the diagnostics that were *dropped* along
   the way (every list in ``hm.dropped_vars``), each outlined in the (p1, p2)
   square of its own two sensitive parameters. Variables sharing a pair are
   overlaid in the same panel, so each figure answers "what would these
   diagnostics have said about this pair, had they been kept?". Inline pages in
   notebook mode; PNG pages plus one multi-page PDF in python mode.

Nothing in here mutates the case: every function only reads ``hm``.

Typical use::

    from proj2dhullsampler.history_matching_animation import (
        list_pairs, animate_pair_constraints, plot_constraint_interlock,
        plot_dropped_constraints,
    )

    list_pairs(hm)                             # index -> pair -> n diagnostics
    ani = animate_pair_constraints(hm, pair=0) # notebook: inline player
    plot_constraint_interlock(hm)              # -> <case>/diagnostics/*.pdf
                                               #    (needs hm.draw() first)
    plot_dropped_constraints(hm)               # -> <case>/diagnostics/dropped_vars/
"""
import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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


def _pct(x):
    """Percentage as a short string, readable across the many decades these
    figures span (100 down to a few surviving points in 10^7)."""
    if x >= 99.95:
        return "100"
    if x >= 10:
        return f"{x:.0f}"
    if x >= 1:
        return f"{x:.1f}"
    if x >= 0.01:
        return f"{x:.2f}"
    return f"{x:.1e}"


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


def _legend(fig_or_ax, cut_label, keep_label, fontsize=7,
            all_label="already ruled out", **kwargs):
    """Proxy-artist legend describing the three point categories."""
    from matplotlib.lines import Line2D

    handles = [
        Line2D([], [], ls="", marker="o", ms=4, mfc=COLOR_ALL, mec="none",
               label=all_label),
        Line2D([], [], ls="", marker="o", ms=4, mfc=COLOR_CUT, mec="none",
               label=cut_label),
        Line2D([], [], ls="", marker="o", ms=4, mfc=COLOR_NOW, mec="none",
               label=keep_label),
        Line2D([], [], ls="--", color=COLOR_HULL, lw=1.2, label="alpha-shape hull"),
    ]
    return fig_or_ax.legend(handles=handles, fontsize=fontsize, **kwargs)


def _build_pair_figure(
    hm,
    pair,
    max_points=3000,
    var_order="given",
    show_hull=True,
    random_state=0,
    figsize=(11, 5.5),
    point_size=6,
):
    """The figure and its per-step draw callable, with no animation attached.

    Returns ``(fig, draw_step, n_steps)`` where ``draw_step(k)`` renders step k
    (0 = unconstrained, n_steps = every diagnostic applied). Both the animation
    and the still-frame writer are built on top of this, so they cannot drift
    apart.
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

    return fig, update, len(vars_)


def build_pair_animation(
    hm,
    pair,
    interval=900,
    hold_last=3,
    **kwargs,
):
    """Build (but do not display or save) the constraint animation for one pair.

    Frame 0 shows every emulated sample filling the (p1, p2) square; frame k
    shows the samples surviving the first k diagnostics of the pair, with the
    diagnostics listed down the right-hand side as they are switched on. Points
    ruled out by the diagnostic added at frame k are highlighted so each step
    shows *what it cost*. Returns ``(fig, FuncAnimation)``.

    `hold_last` repeats the final frame that many extra times so the fully
    constrained region stays on screen for a beat before the loop restarts.
    Remaining keyword arguments go to :func:`_build_pair_figure` (`max_points`,
    `var_order`, `show_hull`, `random_state`, `figsize`, `point_size`).
    """
    fig, draw_step, n_steps = _build_pair_figure(hm, pair, **kwargs)

    frames = list(range(n_steps + 1)) + [n_steps] * max(hold_last, 0)
    ani = animation.FuncAnimation(
        fig, draw_step, frames=frames, interval=interval, blit=False, repeat=True,
    )
    # Render an arbitrary step without reaching into FuncAnimation's private ._func.
    ani.draw_step = draw_step
    ani.n_steps = n_steps
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
        # Number by the pair's real position in the sampling order, so that
        # asking for one pair gives the same filename as doing all of them.
        path = outdir / f"pair_{pairs.index(p):02d}_{p[0]}__{p[1]}.{fmt}"
        writer = "pillow" if fmt == "gif" else "ffmpeg"
        ani.save(str(path), writer=writer, fps=fps, dpi=dpi)
        plt.close(fig)
        paths.append(path)
        print(f"[{i + 1}/{len(todo)}] wrote {path}")

    return paths


def save_pair_frames(
    hm,
    pair=None,
    outdir=None,
    fmt="png",
    dpi=110,
    per_pair_subdir=True,
    **kwargs,
):
    """Write every step of the animation as a separate still image.

    The same frames the animation plays through, but as individual files:
    ``step_00`` is the unconstrained square and ``step_k`` is the state after
    the k-th diagnostic, so flipping through them is the animation frame by
    frame. Useful when you want to drop a particular step into a slide or a
    paper, or to page through the steps without a player.

    `pair` takes an index / tuple / "(p1, p2)" string; leave it out to do every
    parameter pair. `outdir` defaults to ``<case>/diagnostics/frames/``, with
    one subdirectory per pair unless `per_pair_subdir=False`. Extra keyword
    arguments go to :func:`build_pair_animation` (`max_points`, `var_order`,
    `show_hull`, `figsize`, ...).

    Returns the list of written paths.
    """
    _require_prepared(hm)

    if outdir is None:
        outdir = Path(hm.diagnostics_dir) / "frames"
    outdir = Path(outdir)

    pairs = pair_sequence(hm)
    todo = [resolve_pair(hm, pair)] if pair is not None else pairs

    paths = []
    for i, p in enumerate(todo):
        # Deliberately not build_pair_animation(): a FuncAnimation that is never
        # played or saved warns on garbage collection, and stills do not need one.
        fig, draw_step, n_steps = _build_pair_figure(hm, p, **kwargs)

        # Numbered by real position in the sampling order (see above).
        stem = f"pair_{pairs.index(p):02d}_{p[0]}__{p[1]}"
        if per_pair_subdir:
            target = outdir / stem
            name = "step_{k:02d}"
        else:
            target = outdir
            name = stem + "_step_{k:02d}"
        target.mkdir(parents=True, exist_ok=True)

        for k in range(n_steps + 1):
            draw_step(k)
            path = target / (name.format(k=k) + f".{fmt}")
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
            paths.append(path)

        plt.close(fig)
        print(f"[{i + 1}/{len(todo)}] {n_steps + 1} frames -> {target}")

    return paths


# ---------------------------------------------------------------------------
# 2. Alone vs interlocked (two columns per pair)
# ---------------------------------------------------------------------------


def _drawn_samples(hm):
    """The parameter sets the sampler drew, on the 0-1 scale the hulls live on.

    Strictly ``hm.results.unscaled_samples`` - i.e. this figure can only be made
    *after* ``hm.draw()``. The right-hand column is meant to be the run's actual
    output and nothing else, so there is deliberately no fallback that would
    quietly plot something adjacent to it (a saved CSV from another run, a
    re-draw with different settings) under the same colours.
    """
    drawn = getattr(getattr(hm, "results", None), "unscaled_samples", None)
    if drawn is None or not len(drawn):
        raise AttributeError(
            "hm.results.unscaled_samples is empty; plot_constraint_interlock "
            "plots the parameter sets the sampler drew, so it has to run after "
            "hm.draw(). To use a finished run's saved output instead, put it "
            "back on the normalized scale and assign it: "
            "hm.results.unscaled_samples = (real - ppe.min()) / (ppe.max() - ppe.min())."
        )

    missing = [p for p in hm.para_nm if p not in drawn.columns]
    if missing:
        raise KeyError(f"hm.results.unscaled_samples is missing parameters {missing}.")
    return drawn[list(hm.para_nm)].astype(float)


def _cloud(pts):
    """Convex hull of the drawn samples, as a shapely polygon (or None).

    A convex hull is a deliberately crude stand-in for "the region these samples
    occupy": it is defined for as few as three points (an alpha-shape of a
    couple of hundred scattered draws is not), and it never collapses on
    sampling noise. It rounds *up* - a cloud with a hole in it still reports the
    area of its outline - so anything measured against it is an upper bound on
    how much of its own region a pair actually gets to use.
    """
    from shapely import MultiPoint, prepare

    if pts.shape[0] < 3:
        return None
    cloud = MultiPoint([tuple(xy) for xy in pts]).convex_hull
    if cloud.is_empty or cloud.geom_type not in ("Polygon", "MultiPolygon"):
        return None               # every sample collinear or coincident

    prepare(cloud)  # builds an internal index for the repeated contains() below
    return cloud


def _interlock_legend(fig, fontsize=8, **kwargs):
    """Legend for the two-column alone-vs-interlocked figure."""
    from matplotlib.lines import Line2D

    # The drawn samples are marked exactly like the left panel's survivors, so
    # they share one entry rather than getting an identical-looking second one.
    handles = [
        Line2D([], [], ls="", marker="o", ms=4, mfc=COLOR_ALL, mec="none",
               label="ruled out by this pair's own diagnostics"),
        Line2D([], [], ls="", marker="o", ms=4, mfc=COLOR_NOW, mec="none",
               label="left: survivors of those diagnostics  |  "
                     "right: the parameter sets drawn"),
        Line2D([], [], ls="", marker="o", ms=4, mfc=COLOR_CUT, mec="none",
               label="right: those same survivors, as the backdrop"),
        Line2D([], [], ls="--", color=COLOR_HULL, lw=1.2, label="alpha-shape hull"),
    ]
    return fig.legend(handles=handles, fontsize=fontsize, **kwargs)


def plot_constraint_interlock(
    hm,
    pairs=None,
    n_rows=None,
    sort_by="order",
    max_points=1500,
    panel_size=2.4,
    point_size=2.5,
    sample_size=None,
    show_hull=True,
    random_state=0,
    outpath=None,
    save=True,
    dpi=150,
    verbose=True,
):
    """Two panels per parameter pair: constrained alone, then constrained jointly.

    One row per parameter pair, both panels drawn on the same axes (the pair's
    first parameter on x, its second on y, both 0-1 normalized):

    * **left** - the pair *on its own*. The emulated samples ``hm.p_emu``, split
      by that pair's own diagnostics alone: dark blue where every diagnostic in
      ``hm.paras_vars[pair]`` passes, grey where at least one fails. This is the
      whole region the pair's own constraints allow, ignoring every other pair.
    * **right** - the identical square, with those same survivors now drawn in
      salmon as a backdrop and the parameter sets the sampler actually *drew*
      scattered on top in dark blue. Those had to clear every pair at once, so
      salmon left uncovered is region this pair permits but the rest of the
      constraints do not.

    The left panel is therefore the constraint the pair asserts, the right panel
    is the constraint it experiences, and the gap between them is what the
    interlocking with the other pairs costs it.

    **Where each column's points come from.** The left column is
    ``hm.p_emu`` masked with ``hm.tf_masks`` - the emulated design and the raw
    per-diagnostic verdicts, not the smoothed alpha-shape fitted to them, so a
    spur or a gap that the hull rounds off is visible here. (The hull is still
    outlined, dashed, for comparison: that outline is what the sampler actually
    draws from.) The right column is ``hm.results.unscaled_samples``, so this
    figure only works *after* ``hm.draw()`` - see :func:`_drawn_samples`.

    A consequence worth remembering: the right column carries however many
    samples the run produced (often a few hundred), and they are a *sample* of
    the joint feasible region, not a picture of it. The printed "% of own" is
    the share of this pair's own survivors that fall inside the convex hull of
    those drawn samples - a coarse upper bound on how much of its own region the
    pair actually gets to use, and the number ``sort_by="shrinkage"`` orders on.
    A pair at 80% is essentially independent; a pair at 5% is being carved
    almost entirely by constraints living in other parameter subspaces.

    `pairs` / `n_rows` select which pairs get a row (the drawn samples are
    whatever the run produced either way - they do not change with this list).
    `sort_by` is "order" (the order the sampler applies the pairs) or
    "shrinkage" (most interlocked first). The two ``p_emu`` categories are
    thinned to `max_points` *independently*, so a category holding only a sliver
    of the design still shows up; the trade-off is that point density is not
    comparable between colours - the printed counts are. The drawn samples are
    never thinned. `sample_size` defaults to `point_size`, so the drawn samples
    are marked exactly like every other point in the figure.

    Returns the figure; also writes a PDF unless `save=False`.
    """
    _require_prepared(hm)

    if sort_by not in ("order", "shrinkage"):
        raise ValueError(f"sort_by must be 'order' or 'shrinkage', got {sort_by!r}")

    seq = list(pair_sequence(hm))
    pairs = list(seq) if pairs is None else [resolve_pair(hm, p) for p in pairs]
    if n_rows is not None:
        pairs = pairs[: int(n_rows)]
    n = len(pairs)
    if n == 0:
        raise ValueError("No parameter pairs to plot.")

    drawn = _drawn_samples(hm)
    if verbose:
        print(f"{len(drawn):,} drawn samples from hm.results.unscaled_samples")
        print(f"{len(hm.p_emu):,} emulated samples split by {n} pairs' diagnostics")

    # One fixed shuffle shared by every pair, so "the first max_points
    # survivors" is a reproducible unbiased thinning rather than whatever order
    # the design happens to be stored in.
    rng = np.random.default_rng(random_state)
    perm = rng.permutation(len(hm.p_emu))

    from shapely import contains, points as _points

    rows = []
    for p in pairs:
        vars_ = _pair_vars(hm, p)
        xy, masks = _xy_and_masks(hm, p, vars_)
        xy, masks = xy[perm], masks[perm]
        keep = masks.all(axis=1) if masks.shape[1] else np.ones(xy.shape[0], bool)

        pts = drawn[[p[0], p[1]]].to_numpy(dtype=float)
        cloud = _cloud(pts)
        # Measured against the same population the left panel shows: of the
        # emulated samples this pair allows, how many lie where the sampler
        # actually ended up drawing?
        inside = (
            float(contains(cloud, _points(xy[keep])).mean())
            if cloud is not None and keep.any() else 0.0
        )

        rows.append(dict(
            pair=p,
            idx=seq.index(p) if p in seq else None,
            n_vars=len(vars_),
            cut=_thin(xy, ~keep, max_points),
            kept=_thin(xy, keep, max_points),
            drawn=pts,
            n_kept=int(keep.sum()),
            frac_emu=100.0 * float(keep.mean()),
            frac_own=100.0 * inside,
        ))

    order = list(range(n))
    if sort_by == "shrinkage":
        order.sort(key=lambda j: (rows[j]["frac_own"], j))

    fig, axes = plt.subplots(
        n, 2,
        figsize=(2 * panel_size + 2.0, panel_size * n + 1.6),
        squeeze=False,
    )
    sample_size = point_size if sample_size is None else sample_size

    for row, j in enumerate(order):
        r = rows[j]
        p = r["pair"]

        # The left panel's blue survivors are the right panel's salmon backdrop:
        # the same points, so the eye can carry the region across and see what
        # the drawn samples do and do not cover.
        layers = (
            ((r["cut"], COLOR_ALL, 1), (r["kept"], COLOR_NOW, 3)),
            ((r["cut"], COLOR_ALL, 1), (r["kept"], COLOR_CUT, 2)),
        )
        labels = (
            f"{r['n_kept']:,} survive\n{_pct(r['frac_emu'])}% of p_emu",
            f"{len(r['drawn']):,} drawn\n{_pct(r['frac_own'])}% of own",
        )

        for k in (0, 1):
            ax = axes[row][k]
            _style_axes(ax, p, labels=False)
            for pts, color, z in layers[k]:
                if pts.shape[0]:
                    ax.scatter(
                        pts[:, 0], pts[:, 1],
                        s=point_size, c=color, linewidths=0, zorder=z,
                    )
            if k == 1 and r["drawn"].shape[0]:
                ax.scatter(
                    r["drawn"][:, 0], r["drawn"][:, 1],
                    s=sample_size, c=COLOR_NOW, linewidths=0, zorder=5,
                )
            if show_hull:
                _plot_hull(ax, _hull_for(hm, p), color=COLOR_HULL, lw=0.9,
                           ls="--", zorder=4)
            ax.set_xlabel(_short(p[0], 22), fontsize=7.5, labelpad=2)
            ax.text(
                0.03, 0.97, labels[k],
                transform=ax.transAxes, va="top", ha="left", fontsize=7,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85),
            )
            if row == 0:
                ax.set_title(
                    ("constrained by its own diagnostics alone" if k == 0
                     else "what the joint constraint actually leaves"),
                    fontsize=9, pad=6,
                )

        axes[row][0].set_ylabel(_short(p[1], 22), fontsize=7.5, labelpad=2)
        # Row stamp: where the pair sits in the sampling order and how many
        # diagnostics it carries - the two things not readable off the panels.
        axes[row][1].text(
            1.05, 0.5,
            (f"pair {r['idx']}\n" if r["idx"] is not None else "")
            + f"{r['n_vars']} diag",
            transform=axes[row][1].transAxes, va="center", ha="left",
            fontsize=8, color="#555555",
        )

    head = min(1.0 / (panel_size * n + 1.6), 0.4)
    fig.suptitle(
        f"Constrained alone vs interlocked - {getattr(hm, 'case_name', '')}\n"
        "left: hm.p_emu under this pair's own diagnostics   |   "
        "right: the same survivors under the parameter sets the sampler drew\n"
        f"{len(drawn):,} drawn samples against {len(hm.p_emu):,} emulated ones; "
        "axes are 0-1 normalized parameters, x = first, y = second",
        fontsize=11, y=1.0 - 0.06 * head, va="top",
    )
    _interlock_legend(
        fig, fontsize=8, loc="upper center", ncol=2, frameon=False,
        bbox_to_anchor=(0.5, 1.0 - 0.55 * head),
    )
    fig.tight_layout(rect=[0, 0, 1, 1.0 - head])

    if save:
        if outpath is None:
            outpath = Path(hm.diagnostics_dir) / "constraint_interlock.pdf"
        outpath = Path(outpath)
        outpath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
        print(f"wrote {outpath}")

    return fig


# ---------------------------------------------------------------------------
# 3. Dropped diagnostics, in the square of their own sensitive parameters
# ---------------------------------------------------------------------------

# Categorical slots in fixed order, paired with a line style so an outline
# stays identifiable without relying on colour alone.
DROPPED_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#4a3aa7", "#e34948"]
DROPPED_STYLES = ["-", "--", "-", "--", "-", "--"]

# Human-readable names for the lists HistoryMatching stores in hm.dropped_vars.
DROP_REASONS = {
    "by_name": "dropped by name",
    "by_emulator_performance": "emulator error",
    "useless": "no constraint",
    "tight": "too tight",
    "nooverlap2d": "no 2-D overlap",
    "too_few_vars": "too few in pair",
    "during_iteration": "dropped in orchestrate",
}


def _dropped_source(hm, dropped):
    """The ``{reason: [vars] | {pair: [vars]}}`` mapping to plot.

    `dropped` may be that mapping itself, or a path to a saved
    ``<result>_dropped_vars.json``. Left as None, the live ``hm.dropped_vars``
    is used if it holds anything; otherwise (e.g. a case reloaded from disk
    rather than re-run) the run's JSON under ``<case>/output/`` is read.
    """
    if isinstance(dropped, dict):
        return dropped
    if dropped is not None:
        return json.loads(Path(dropped).read_text())

    live = vars(getattr(hm, "dropped_vars", None) or type("_", (), {})())
    if any(len(v) for v in live.values()):
        return live

    out = Path(hm.root) / "output"
    result_name = getattr(hm, "result_name", None)
    if result_name and (out / f"{result_name}_dropped_vars.json").exists():
        path = out / f"{result_name}_dropped_vars.json"
    else:
        found = sorted(out.glob("*_dropped_vars.json"))
        if len(found) != 1:
            raise FileNotFoundError(
                "hm.dropped_vars is empty and no single *_dropped_vars.json was "
                f"found under {out} (found {[p.name for p in found]}). Pass "
                "dropped=<path to the JSON> or dropped=<dict>."
            )
        path = found[0]
    print(f"hm.dropped_vars is empty; reading {path}")
    return json.loads(path.read_text())


def _full_masks(hm):
    """Pass/fail masks for *every* diagnostic, dropped ones included.

    ``hm.tf_masks`` loses a column at each drop; ``hm.tf_masks_raw`` is the copy
    ``load_mask`` took before any of that.
    """
    raw = getattr(hm, "tf_masks_raw", None)
    if raw is not None:
        return raw
    level = getattr(getattr(hm, "specifications", None), "uncertainty_threshold", None)
    if level is None:
        raise AttributeError("hm has no tf_masks_raw; call hm.load_mask(...) first.")
    return pd.read_csv(Path(hm.root) / f"tf_masks_level_{level}.csv", index_col=0)


def collect_dropped_vars(hm, dropped=None, reasons=None, verbose=True):
    """One row per dropped diagnostic: why it went, its pair, how much passes.

    The sensitive parameters come from ``<case>/meta.csv`` rather than
    ``hm.meta``, because every drop step also removes the variable from
    ``hm.meta``. The pair is ordered the way ``group_para_climatology`` orders
    it (by parameter index), so a dropped variable lands on the same (p1, p2)
    key as the kept diagnostics of that pair.

    `reasons` restricts to some of the ``hm.dropped_vars`` lists (e.g.
    ``["tight", "useless"]``). Returns a DataFrame with columns
    ``var, reason, p1, p2, n_pass, pct_pass``.
    """
    source = _dropped_source(hm, dropped)
    meta = pd.read_csv(Path(hm.root) / "meta.csv", index_col=0)
    masks = _full_masks(hm)

    rows, seen, skipped = [], set(), []
    for reason, entry in source.items():
        if reasons is not None and reason not in reasons:
            continue
        # during_iteration is {pair: [vars]}; everything else is a flat list
        names = [v for vs in entry.values() for v in vs] if isinstance(entry, dict) else entry
        for v in names:
            if v in seen:
                continue
            seen.add(v)
            if v not in meta.columns or v not in masks.columns:
                skipped.append(v)
                continue
            inds = sorted(int(i) for i in meta[v].to_numpy())
            if len(inds) != 2:
                skipped.append(v)
                continue
            m = masks[v].to_numpy(dtype=bool)
            rows.append(dict(
                var=v, reason=reason,
                p1=hm.para_nm[inds[0]], p2=hm.para_nm[inds[1]],
                n_pass=int(m.sum()), pct_pass=100.0 * float(m.mean()),
            ))

    if skipped:
        warnings.warn(
            f"{len(skipped)} dropped variables have no 2-parameter entry in "
            f"meta.csv or no mask column; skipping them: {skipped}",
            stacklevel=2,
        )

    table = pd.DataFrame(rows, columns=["var", "reason", "p1", "p2", "n_pass", "pct_pass"])
    if verbose:
        n_pairs = table.groupby(["p1", "p2"]).ngroups if len(table) else 0
        print(f"{len(table)} dropped variables over {n_pairs} parameter pairs")
        if len(table):
            print(table["reason"].value_counts().to_string())
    return table


def _pass_fraction(xy, mask, grid):
    """Share of emulated samples passing, per cell of a `grid` x `grid` binning.

    Each emulator only sees its own two parameters, so within a cell the mask is
    essentially constant and the 0.5 contour of this field is the variable's
    pass region.
    """
    edges = np.linspace(0.0, 1.0, grid + 1)
    total, _, _ = np.histogram2d(xy[:, 0], xy[:, 1], bins=[edges, edges])
    passed, _, _ = np.histogram2d(xy[mask, 0], xy[mask, 1], bins=[edges, edges])
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(total > 0, passed / total, 0.0)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, frac.T  # .T: histogram2d is (x, y); contour wants (row=y, col=x)


def _dropped_panels(hm, table, max_vars_per_panel):
    """Panels to draw: pairs with the most dropped variables first, each pair
    split into chunks small enough to tell the outlines apart."""
    kept = {tuple(p) for p in getattr(hm, "paras_vars", {})}
    groups = sorted(
        table.groupby(["p1", "p2"], sort=False),
        key=lambda kv: (-len(kv[1]), kv[0]),
    )
    panels = []
    for pair, g in groups:
        g = g.sort_values("pct_pass", ascending=False)  # loosest drawn first, underneath
        chunks = [g.iloc[i:i + max_vars_per_panel] for i in range(0, len(g), max_vars_per_panel)]
        for c, chunk in enumerate(chunks):
            panels.append(dict(
                pair=tuple(pair), rows=chunk, n_total=len(g),
                chunk=c, n_chunks=len(chunks), kept=tuple(pair) in kept,
            ))
    return panels


def _panel_label_lines(panel, show_hull):
    """Height of a panel's legend in text lines (for sizing the legend row)."""
    return 2 * len(panel["rows"]) + (1 if show_hull and panel["kept"] else 0)


def _draw_dropped_panel(ax, ax_leg, hm, panel, xy_cache, masks, grid,
                        show_hull, fontsize):
    from matplotlib.colors import to_rgba
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    pair = panel["pair"]
    if pair not in xy_cache:
        xy_cache[pair] = hm.p_emu[list(pair)].to_numpy(dtype=float)
    xy = xy_cache[pair]

    _style_axes(ax, pair, labels=True, fontsize=fontsize + 1)
    handles = []

    for i, r in enumerate(panel["rows"].itertuples()):
        color, ls = DROPPED_COLORS[i], DROPPED_STYLES[i]
        mask = masks[r.var].to_numpy(dtype=bool)
        label = f"{r.var}\n   {DROP_REASONS.get(r.reason, r.reason)} · {_pct(r.pct_pass)}% pass"

        centers, frac = _pass_fraction(xy, mask, grid)
        if (frac >= 0.5).sum() >= 4:
            # light fill says which side of the line passes; the line carries identity
            ax.contourf(centers, centers, frac, levels=[0.5, 1.01],
                        colors=[color], alpha=0.07, zorder=2)
            ax.contour(centers, centers, frac, levels=[0.5],
                       colors=[color], linestyles=[ls], linewidths=1.8, zorder=3 + i)
            handles.append(Patch(fc=to_rgba(color, 0.15), ec=color, ls=ls, lw=1.5,
                                 label=label))
        else:
            # too small a region to resolve on the grid: show the survivors themselves
            pts = xy[mask]
            ax.scatter(pts[:, 0], pts[:, 1], s=10, marker="x", c=color,
                       linewidths=0.9, zorder=3 + i)
            handles.append(Line2D([], [], ls="", marker="x", color=color, label=label))

    if show_hull and panel["kept"]:
        before = len(ax.lines)
        _plot_hull(ax, _hull_for(hm, pair), color=COLOR_HULL, lw=1.1, ls=":", zorder=10)
        if len(ax.lines) > before:
            handles.append(Line2D([], [], ls=":", color=COLOR_HULL, lw=1.2,
                                  label="kept diagnostics (alpha-shape)"))

    status = "pair kept" if panel["kept"] else "pair not kept"
    part = (f"  ·  part {panel['chunk'] + 1}/{panel['n_chunks']}"
            if panel["n_chunks"] > 1 else "")
    ax.set_title(f"{panel['n_total']} dropped with a constraint  ·  {status}{part}",
                 fontsize=fontsize + 1, loc="left", pad=4)

    ax_leg.axis("off")
    ax_leg.legend(handles=handles, loc="upper left", fontsize=fontsize,
                  frameon=False, borderaxespad=0.0, handlelength=2.2,
                  labelspacing=0.5)


def plot_dropped_constraints(
    hm,
    dropped=None,
    reasons=None,
    pairs=None,
    mode=None,
    outdir=None,
    ncols=4,
    rows_per_page=3,
    max_vars_per_panel=6,
    grid=60,
    show_hull=True,
    panel_size=3.0,
    fontsize=7,
    dpi=150,
    verbose=True,
):
    """Where each *dropped* diagnostic would have constrained its own pair.

    For every variable in ``hm.dropped_vars`` (or the saved
    ``*_dropped_vars.json``, see :func:`collect_dropped_vars`), its two
    sensitive parameters are read from ``<case>/meta.csv`` and its pass region
    is drawn on that (p1, p2) square: the emulated samples ``hm.p_emu`` binned
    on a `grid` x `grid` lattice, outlined where at least half of a cell passes
    in ``hm.tf_masks_raw``. Variables sharing a pair are overlaid in one panel,
    each with its own colour and line style; the legend under the panel names
    the variable, why it was dropped and what share of ``p_emu`` it lets
    through. A region too small to outline on the grid has its surviving
    samples marked with crosses instead.

    Variables that do not constrain anything - passing for every emulated
    sample, or for none - are left out (the count is printed), and a pair left
    with no variables gets no panel.

    If the pair still carries kept diagnostics, their alpha-shape is drawn
    dotted, so you can see whether a dropped variable disagreed with the
    constraint that was actually used.

    Pairs with the most dropped variables come first. A pair with more than
    `max_vars_per_panel` of them is split over consecutive panels ("part 1/2").
    `pairs` restricts to some pairs (index / tuple / "(p1, p2)" string, either
    order); `reasons` restricts to some of the drop lists.

    Notebook mode shows the pages inline and returns the figures. Python mode
    writes ``page_XX.png`` files and one ``dropped_constraints.pdf`` into
    `outdir` (default ``<case>/diagnostics/dropped_vars/``) and returns the
    paths.
    """
    from matplotlib.backends.backend_pdf import PdfPages

    mode = _resolve_mode(hm, mode)
    table = collect_dropped_vars(hm, dropped=dropped, reasons=reasons, verbose=verbose)
    if table.empty:
        print("No dropped variables to plot.")
        return []

    if pairs is not None:
        wanted = set()
        for p in pairs:
            if isinstance(p, (int, np.integer)):
                wanted.add(tuple(pair_sequence(hm)[int(p)]))
                continue
            if isinstance(p, str):
                p = [s.strip() for s in p.strip("() ").split(",")]
            wanted.update({tuple(p), tuple(p)[::-1]})
        table = table[[(a, b) in wanted for a, b in zip(table["p1"], table["p2"])]]
        if table.empty:
            raise ValueError(f"None of the dropped variables belong to pairs {pairs}.")

    masks = _full_masks(hm)
    if len(masks) != len(hm.p_emu):
        raise ValueError(
            f"p_emu has {len(hm.p_emu)} rows but the masks have {len(masks)}; "
            "they must describe the same emulated samples."
        )

    no_constraint = (table["n_pass"] == 0) | (table["n_pass"] == len(masks))
    if verbose and no_constraint.any():
        print(f"Leaving out {int(no_constraint.sum())} dropped variables that pass "
              "everywhere or nowhere (no constraint to draw)")
    table = table[~no_constraint]
    if table.empty:
        print("No dropped variable constrains its parameter pair; nothing to plot.")
        return []

    panels = _dropped_panels(hm, table, max_vars_per_panel)
    per_page = ncols * rows_per_page
    pages = [panels[i:i + per_page] for i in range(0, len(panels), per_page)]
    case = getattr(hm, "case_name", "")

    if mode == "python":
        outdir = Path(hm.diagnostics_dir) / "dropped_vars" if outdir is None else Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        for old in outdir.glob("page_*.png"):  # a rerun may produce fewer pages
            old.unlink()
        pdf = PdfPages(outdir / "dropped_constraints.pdf")

    xy_cache, figs, paths = {}, [], []
    line_h = fontsize * 1.3 / 72.0  # inches per legend text line, spacing included
    hspace, top, bottom = 0.12, 0.98, 0.01
    try:
        for n, page in enumerate(pages):
            n_rows = -(-len(page) // ncols)
            leg_h = [
                0.2 + line_h * max(_panel_label_lines(p, show_hull)
                                   for p in page[r * ncols:(r + 1) * ncols])
                for r in range(n_rows)
            ]

            heights = [h for r in range(n_rows) for h in (panel_size, leg_h[r])]
            # gridspec spacing is a fraction of the mean row height, on top of the rows
            gaps = hspace * sum(heights) / len(heights) * (len(heights) - 1)
            fig = plt.figure(figsize=(ncols * (panel_size + 0.6),
                                      (sum(heights) + gaps) / (top - bottom)))
            # the suptitle sits above the grid and bbox_inches="tight" keeps it,
            # so the grid can use the whole figure
            gs = fig.add_gridspec(2 * n_rows, ncols, height_ratios=heights,
                                  hspace=hspace, wspace=0.45,
                                  top=top, bottom=bottom, left=0.06, right=0.98)

            for k, panel in enumerate(page):
                r, c = divmod(k, ncols)
                ax = fig.add_subplot(gs[2 * r, c])
                ax_leg = fig.add_subplot(gs[2 * r + 1, c])
                _draw_dropped_panel(ax, ax_leg, hm, panel, xy_cache, masks,
                                    grid, show_hull, fontsize)

            fig.suptitle(
                f"Dropped diagnostics in their own sensitive-parameter square - {case}"
                f"   (page {n + 1}/{len(pages)})\n"
                "shaded/outlined: where the variable passes (hm.tf_masks_raw on hm.p_emu); "
                "axes are 0-1 normalized parameters",
                fontsize=fontsize + 3, y=1.0, va="bottom",
            )

            if mode == "python":
                path = outdir / f"page_{n + 1:02d}.png"
                fig.savefig(path, dpi=dpi, bbox_inches="tight")
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
                paths.append(path)
                if verbose:
                    print(f"[{n + 1}/{len(pages)}] wrote {path}")
            else:
                plt.show()
                figs.append(fig)
    finally:
        if mode == "python":
            pdf.close()

    if mode == "python":
        paths.append(outdir / "dropped_constraints.pdf")
        if verbose:
            print(f"wrote {paths[-1]}")
        return paths
    return figs
