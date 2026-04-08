"""
Discrete-Event Simulation – Computer Processing System
=======================================================
Implements the full 12-step simulation project cycle described in
projectOutline.md.

System architecture
-------------------
  Arrivals  ──►  Buffer 1 (Type I)  ─┐
                                      ├─► Node 1 processor ──► Router ──► Node 2 (2 processors)
  Arrivals  ──►  Buffer 2 (Type II) ─┘

Node 1 alternates between Buffer 1 (50 ms window) and Buffer 2 (30 ms window).
A current packet is never interrupted – the switch waits until service completes.

Usage
-----
    python simulation.py

Outputs:
  • Console report with point estimates and 95% confidence intervals
  • Chi-square goodness-of-fit results
  • histograms.png  – six histogram panels with fitted PDF overlays
"""

import heapq
import math
import os
import random

import matplotlib
matplotlib.use("Agg")          # non-interactive backend (safe for all environments)
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# ---------------------------------------------------------------------------
# Event-type constants (integers kept small for fast comparisons)
# ---------------------------------------------------------------------------
ARRIVAL_1    = 0   # Type I packet arrives at Buffer 1
ARRIVAL_2    = 1   # Type II packet arrives at Buffer 2
DEPARTURE_N1 = 2   # Node 1 processor finishes a packet
DEPARTURE_N2 = 3   # One of Node 2's processors finishes
SWITCH_N1    = 4   # Node 1 buffer-phase timer fires

# Phase durations (ms)
PHASE_DUR = {1: 50.0, 2: 30.0}


# ---------------------------------------------------------------------------
# Random variate generators  (inverse-transform method)
# ---------------------------------------------------------------------------

def _rvg_interarrival_1(rng):
    """Type I inter-arrival: X = -4 ln(U)  [mean = 4 ms]"""
    return -4.0 * math.log(rng.random())


def _rvg_interarrival_2(rng):
    """Type II inter-arrival: X = -12 ln(U)  [mean = 12 ms]"""
    return -12.0 * math.log(rng.random())


def _rvg_service_b1(rng):
    """Node 1, Buffer 1 service: X = 1 + 2U  [Uniform(1, 3)]"""
    return 1.0 + 2.0 * rng.random()


def _rvg_service_b2(rng):
    """Node 1, Buffer 2 service: X = 2 + 4U  [Uniform(2, 6)]"""
    return 2.0 + 4.0 * rng.random()


def _rvg_service_n2(rng):
    """Node 2 processor service: X = -5 ln(U)  [Exp, mean = 5 ms]"""
    return -5.0 * math.log(rng.random())


# ---------------------------------------------------------------------------
# Simulation class  (one replication)
# ---------------------------------------------------------------------------

class Simulation:
    """
    Single replication of the discrete-event simulation.

    Parameters
    ----------
    sim_time : float
        Simulation horizon in milliseconds.
    rng : random.Random
        Seeded RNG instance (ensures reproducibility per replication).
    """

    def __init__(self, sim_time, rng):
        self.sim_time = sim_time
        self.rng      = rng

        # Future Event List: min-heap of (time, seq, event_type, data)
        self._fel  = []
        self._seq  = 0            # monotonic tie-breaker

        # ── Clock ──────────────────────────────────────────────────────────
        self.clock = 0.0

        # ── Node 1 state ───────────────────────────────────────────────────
        # Each buffer stores the clock-time at which the packet arrived
        self.buf1 = []            # Type I packets waiting
        self.buf2 = []            # Type II packets waiting

        self.n1_busy           = False   # True while processor is occupied
        self.n1_serving        = None    # 1 or 2 – which buffer is being served
        self.n1_phase          = 1       # current priority phase (1=B1, 2=B2)
        self.n1_switch_pending = False   # True when the phase timer expired
                                         # mid-service (switch deferred)

        # ── Node 2 state ───────────────────────────────────────────────────
        self.n2_buf    = []              # shared buffer: arrival times
        self.n2_busy   = [False, False]  # per-processor busy flags

        # ── Statistical accumulators ───────────────────────────────────────
        self._area_t  = 0.0      # last clock update for area integrals

        # Time-weighted area integrals (Σ queue_length × Δt)
        self.area_q1_wait  = 0.0   # Buffer 1 – waiting packets only
        self.area_q1_total = 0.0   # Buffer 1 – waiting + in service
        self.area_q2_wait  = 0.0
        self.area_q2_total = 0.0
        self.area_qn2      = 0.0   # Node 2 shared buffer (waiting only)

        # Per-packet sample lists
        self.wait_b1     = []   # waiting times in Buffer 1
        self.wait_b2     = []   # waiting times in Buffer 2
        self.wait_n2     = []   # waiting times in Node 2 buffer
        self.sojourn_b1  = []   # wait + service at Node 1 – Type I
        self.sojourn_b2  = []   # wait + service at Node 1 – Type II

        # Packet counters
        self.n_type1      = 0   # Type I arrivals
        self.n_type2      = 0   # Type II arrivals
        self.n_exit_n1    = 0   # packets that leave Node 1 (all types, pre-routing)
        self.n_discard    = 0   # Type II packets discarded (50 % rule)
        self.n_redirect   = 0   # packets redirected by Router 1
        self.n_had_wait_b1 = 0  # Type I packets that found a non-empty/busy B1
        self.n_had_wait_b2 = 0  # Type II packets that found a non-empty/busy B2

    # ── Internal helpers ───────────────────────────────────────────────────

    def _schedule(self, time, etype, data=None):
        """Push an event onto the Future Event List."""
        self._seq += 1
        heapq.heappush(self._fel, (time, self._seq, etype, data))

    def _update_areas(self):
        """
        Accumulate time-weighted queue-length integrals.
        Must be called at the START of every event handler, before any
        state modification, so that the integral reflects the state that
        persisted from the previous event up to now.
        """
        dt = self.clock - self._area_t
        if dt <= 0.0:
            return

        q1 = len(self.buf1)
        q2 = len(self.buf2)
        s1 = 1 if (self.n1_busy and self.n1_serving == 1) else 0
        s2 = 1 if (self.n1_busy and self.n1_serving == 2) else 0

        self.area_q1_wait  += q1 * dt
        self.area_q1_total += (q1 + s1) * dt
        self.area_q2_wait  += q2 * dt
        self.area_q2_total += (q2 + s2) * dt
        self.area_qn2      += len(self.n2_buf) * dt

        self._area_t = self.clock

    # ── Main loop ──────────────────────────────────────────────────────────

    def run(self):
        """Execute the simulation until sim_time is reached."""
        # Seed the FEL with first arrivals and the first phase-switch timer
        self._schedule(_rvg_interarrival_1(self.rng), ARRIVAL_1)
        self._schedule(_rvg_interarrival_2(self.rng), ARRIVAL_2)
        self._schedule(PHASE_DUR[1], SWITCH_N1)

        while self._fel:
            time, _, etype, data = heapq.heappop(self._fel)
            if time > self.sim_time:
                break

            self._update_areas()   # update integrals BEFORE changing state
            self.clock = time

            if   etype == ARRIVAL_1:    self._on_arrival_1()
            elif etype == ARRIVAL_2:    self._on_arrival_2()
            elif etype == DEPARTURE_N1: self._on_departure_n1(data)
            elif etype == DEPARTURE_N2: self._on_departure_n2(data)
            elif etype == SWITCH_N1:    self._on_switch_n1()

    # ── Event handlers ─────────────────────────────────────────────────────

    def _on_arrival_1(self):
        """Type I packet arrives at Buffer 1."""
        self.n_type1 += 1
        # Schedule the next Type I arrival
        self._schedule(self.clock + _rvg_interarrival_1(self.rng), ARRIVAL_1)

        # Packet must wait if Node 1 is already busy (serving B1 or B2)
        if self.n1_busy:
            self.n_had_wait_b1 += 1

        self.buf1.append(self.clock)   # enqueue with arrival timestamp

        # If Node 1 is idle, start service immediately
        if not self.n1_busy:
            self._start_n1()

    def _on_arrival_2(self):
        """Type II packet arrives at Buffer 2."""
        self.n_type2 += 1
        self._schedule(self.clock + _rvg_interarrival_2(self.rng), ARRIVAL_2)

        if self.n1_busy:
            self.n_had_wait_b2 += 1

        self.buf2.append(self.clock)

        if not self.n1_busy:
            self._start_n1()

    def _start_n1(self):
        """
        Pick the next packet for Node 1 service based on the current phase
        and pending-switch flag, then schedule its departure.
        """
        buf = self._select_buf()
        if buf is None:
            return   # both buffers empty – Node 1 stays idle

        self.n1_busy    = True
        self.n1_serving = buf

        if buf == 1:
            arr_time = self.buf1.pop(0)
            wait     = self.clock - arr_time
            self.wait_b1.append(wait)
            svc = _rvg_service_b1(self.rng)
            self.sojourn_b1.append(wait + svc)
        else:
            arr_time = self.buf2.pop(0)
            wait     = self.clock - arr_time
            self.wait_b2.append(wait)
            svc = _rvg_service_b2(self.rng)
            self.sojourn_b2.append(wait + svc)

        self._schedule(self.clock + svc, DEPARTURE_N1, buf)

    def _select_buf(self):
        """
        Determine which buffer to serve next.

        Rules (in priority order):
        1. If a phase-switch is pending, try to serve the NEW phase's buffer.
           If that buffer is empty, serve the other (or None if both empty).
        2. Otherwise serve the CURRENT phase's buffer; fall back to the other.
        Returns the buffer number (1 or 2) or None if both are empty.
        """
        b1 = bool(self.buf1)
        b2 = bool(self.buf2)

        if not b1 and not b2:
            return None

        if self.n1_switch_pending:
            self.n1_switch_pending = False
            target = self.n1_phase   # the phase we just switched TO
            if (target == 1 and b1) or (target == 2 and b2):
                return target
            # Target buffer is empty – serve whichever buffer has packets
            return 1 if b1 else 2

        # No pending switch: honour current phase priority
        if self.n1_phase == 1:
            return 1 if b1 else 2
        else:
            return 2 if b2 else 1

    def _on_departure_n1(self, buf):
        """
        A packet finishes service at Node 1.

        Routing logic
        -------------
        • Type II (buf == 2): discard with probability 0.5.
        • Remaining packets: forward to Node 2 unless Router 1 redirects
          (condition: more than 5 packets already waiting in Node 2 buffer).
        """
        self.n1_busy = False
        self.n_exit_n1 += 1

        # ── Type II discard check ─────────────────────────────────────────
        if buf == 2 and self.rng.random() <= 0.5:
            self.n_discard += 1
            self._start_n1()   # start next Node 1 service before returning
            return

        # ── Router 1: redirect if Node 2 buffer is over capacity ──────────
        if len(self.n2_buf) > 5:
            self.n_redirect += 1
        else:
            self.n2_buf.append(self.clock)   # add to Node 2 shared buffer
            self._start_n2()                 # assign to a free processor if any

        # Start serving the next waiting packet at Node 1
        self._start_n1()

    def _start_n2(self):
        """
        Assign waiting Node 2 packets to any free processors (greedy).
        Called whenever a new packet enters the Node 2 buffer OR a
        processor becomes free.
        """
        for proc in (0, 1):
            if not self.n2_busy[proc] and self.n2_buf:
                arr_time = self.n2_buf.pop(0)
                wait     = self.clock - arr_time
                self.wait_n2.append(wait)
                self.n2_busy[proc] = True
                self._schedule(
                    self.clock + _rvg_service_n2(self.rng),
                    DEPARTURE_N2,
                    proc,
                )

    def _on_departure_n2(self, proc):
        """A Node 2 processor finishes; free it and pull the next queued packet."""
        self.n2_busy[proc] = False
        self._start_n2()

    def _on_switch_n1(self):
        """
        Node 1 buffer-phase timer fires.

        Transition the active phase (1 → 2 or 2 → 1) and schedule the
        next phase-switch event.

        If Node 1 is currently serving the OUTGOING phase's buffer, set
        n1_switch_pending so the switch takes effect after the current
        packet's service completes.
        If Node 1 is idle, immediately start serving from the new phase.
        """
        old_phase       = self.n1_phase
        self.n1_phase   = 3 - old_phase            # toggle: 1↔2
        self._schedule(self.clock + PHASE_DUR[self.n1_phase], SWITCH_N1)

        if self.n1_busy and self.n1_serving == old_phase:
            # Mid-service on the outgoing phase → defer switch
            self.n1_switch_pending = True
        elif not self.n1_busy:
            # Idle → start serving the new phase (or whatever is available)
            self._start_n1()
        # If busy serving the NEW phase already, no action needed

    # ── Results ────────────────────────────────────────────────────────────

    def results(self):
        """
        Return a dict of scalar performance metrics and raw sample lists.
        Scalars are used across replications for confidence intervals;
        raw lists are pooled for histograms and goodness-of-fit tests.
        """
        T = self.sim_time

        def _mean(lst):
            return float(np.mean(lst)) if lst else 0.0

        total_routed = self.n_exit_n1 - self.n_discard  # reached the router

        return {
            # ── Per-packet mean times ─────────────────────────────────────
            "avg_wait_b1":    _mean(self.wait_b1),
            "avg_wait_b2":    _mean(self.wait_b2),
            "avg_wait_n2":    _mean(self.wait_n2),
            "avg_sojourn_b1": _mean(self.sojourn_b1),   # wait + service, Type I
            "avg_sojourn_b2": _mean(self.sojourn_b2),   # wait + service, Type II
            # ── Time-average queue lengths (Little's-Law compatible) ───────
            "avg_q1_wait":  self.area_q1_wait  / T,
            "avg_q1_total": self.area_q1_total / T,
            "avg_q2_wait":  self.area_q2_wait  / T,
            "avg_q2_total": self.area_q2_total / T,
            "avg_qn2":      self.area_qn2      / T,
            # ── Probabilities ─────────────────────────────────────────────
            "redirect_prob": (self.n_redirect / total_routed)
                             if total_routed > 0 else 0.0,
            "p_wait_b1": (self.n_had_wait_b1 / self.n_type1)
                         if self.n_type1 > 0 else 0.0,
            "p_wait_b2": (self.n_had_wait_b2 / self.n_type2)
                         if self.n_type2 > 0 else 0.0,
            # ── Raw sample lists (underscore prefix = not a scalar metric) ─
            "_wait_b1":    self.wait_b1,
            "_wait_b2":    self.wait_b2,
            "_wait_n2":    self.wait_n2,
            "_sojourn_b1": self.sojourn_b1,
            "_sojourn_b2": self.sojourn_b2,
        }


# ---------------------------------------------------------------------------
# Multi-replication runner
# ---------------------------------------------------------------------------

_SCALAR_KEYS = [
    "avg_wait_b1", "avg_wait_b2", "avg_wait_n2",
    "avg_sojourn_b1", "avg_sojourn_b2",
    "avg_q1_wait", "avg_q1_total",
    "avg_q2_wait", "avg_q2_total",
    "avg_qn2",
    "redirect_prob", "p_wait_b1", "p_wait_b2",
]

_POOL_KEYS = ["wait_b1", "wait_b2", "wait_n2", "sojourn_b1", "sojourn_b2"]


def run_replications(n_reps=30, sim_time=100_000.0, base_seed=42):
    """
    Run *n_reps* independent replications and collect per-replication
    scalar metrics and pooled sample lists.

    Parameters
    ----------
    n_reps    : number of independent replications
    sim_time  : simulation horizon per replication (ms)
    base_seed : replication i uses seed (base_seed + i)

    Returns
    -------
    accum : dict  key → list of per-replication scalar values
    pools : dict  key → list of all individual sample values (all reps pooled)
    """
    accum = {k: [] for k in _SCALAR_KEYS}
    pools = {k: [] for k in _POOL_KEYS}

    for i in range(n_reps):
        rng = random.Random(base_seed + i)
        sim = Simulation(sim_time, rng)
        sim.run()
        r = sim.results()

        for k in _SCALAR_KEYS:
            accum[k].append(r[k])
        for k in _POOL_KEYS:
            pools[k].extend(r[f"_{k}"])

    return accum, pools


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------

def ci_t(data, conf=0.95):
    """
    Compute a t-distribution confidence interval.

    Returns (mean, lower_bound, upper_bound).
    """
    n    = len(data)
    mean = float(np.mean(data))
    if n < 2:
        return mean, mean, mean
    se = float(stats.sem(data))
    h  = se * stats.t.ppf((1.0 + conf) / 2.0, df=n - 1)
    return mean, mean - h, mean + h


def chi_square_gof(data, dist_name):
    """
    Chi-square goodness-of-fit test.

    Parameters
    ----------
    data      : 1-D array-like of observed samples
    dist_name : 'expon' or 'uniform'

    Returns (chi2_stat, degrees_of_freedom, p_value).
    """
    arr = np.array(data, dtype=float)
    n   = len(arr)

    # Fit the distribution
    if dist_name == "expon":
        loc, scale = stats.expon.fit(arr, floc=0.0)
        dist_obj   = stats.expon(loc=loc, scale=scale)
        n_params   = 1   # only scale estimated (loc fixed at 0)
    elif dist_name == "uniform":
        loc, scale = stats.uniform.fit(arr)
        dist_obj   = stats.uniform(loc=loc, scale=scale)
        n_params   = 2   # both loc and scale estimated
    else:
        raise ValueError(f"Unsupported distribution: {dist_name!r}")

    # Bin count (Sturges' rule, minimum 10)
    n_bins = max(10, int(1 + 3.322 * math.log10(n)))
    counts, edges = np.histogram(arr, bins=n_bins)

    expected = np.array([
        n * (dist_obj.cdf(edges[i + 1]) - dist_obj.cdf(edges[i]))
        for i in range(n_bins)
    ])

    # Merge bins with expected frequency < 5 (left-to-right, remainder merged last)
    obs_m, exp_m = [], []
    obs_acc, exp_acc = 0, 0.0
    for o, e in zip(counts, expected):
        obs_acc += o
        exp_acc += e
        if exp_acc >= 5.0:
            obs_m.append(obs_acc)
            exp_m.append(exp_acc)
            obs_acc, exp_acc = 0, 0.0
    if obs_acc > 0:
        if obs_m:
            obs_m[-1] += obs_acc
            exp_m[-1] += exp_acc
        else:
            obs_m.append(obs_acc)
            exp_m.append(exp_acc)

    obs_m = np.array(obs_m, dtype=float)
    exp_m = np.array(exp_m, dtype=float)

    k   = len(obs_m)
    df  = max(k - 1 - n_params, 1)
    chi2_stat = float(np.sum((obs_m - exp_m) ** 2 / exp_m))
    p_value   = float(1.0 - stats.chi2.cdf(chi2_stat, df))
    return chi2_stat, df, p_value


# ---------------------------------------------------------------------------
# Histogram plotting
# ---------------------------------------------------------------------------

_PANEL_CFG = [
    # (pool_key,   title,                                        fit_dist,  colour)
    ("wait_b1",    "Waiting Time in Buffer 1 (ms)",              "expon",   "steelblue"),
    ("wait_b2",    "Waiting Time in Buffer 2 (ms)",              "expon",   "darkorange"),
    ("wait_n2",    "Waiting Time at Node 2 Buffer (ms)",         "expon",   "seagreen"),
    ("sojourn_b1", "Total Time at Node 1 \u2013 Type I (ms)",    "expon",   "mediumpurple"),
    ("sojourn_b2", "Total Time at Node 1 \u2013 Type II (ms)",   "uniform", "tomato"),
]


def plot_histograms(pools, output_dir="."):
    """
    Generate histogram panels with fitted-PDF overlays and save to
    *output_dir*/histograms.png.
    """
    os.makedirs(output_dir, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes_flat = axes.flatten()

    for idx, (key, title, fit_dist, colour) in enumerate(_PANEL_CFG):
        ax   = axes_flat[idx]
        data = np.array(pools[key], dtype=float)
        if len(data) == 0:
            ax.set_visible(False)
            continue

        ax.hist(data, bins=60, color=colour, edgecolor="white",
                alpha=0.8, density=True)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel("Density")
        ax.grid(True, alpha=0.3)

        # Fitted PDF overlay
        x = np.linspace(max(0, data.min()), data.max(), 300)
        if fit_dist == "expon":
            loc, scale = stats.expon.fit(data, floc=0.0)
            ax.plot(x, stats.expon.pdf(x, loc=loc, scale=scale),
                    "r-", lw=2, label=f"Exp(μ={scale:.2f})")
        elif fit_dist == "uniform":
            loc, scale = stats.uniform.fit(data)
            ax.plot(x, stats.uniform.pdf(x, loc=loc, scale=scale),
                    "r-", lw=2, label=f"Unif({loc:.2f}, {loc+scale:.2f})")
        ax.legend(fontsize=8)

    axes_flat[-1].set_visible(False)   # sixth panel unused
    fig.suptitle(
        "Computer Processing System – DES Output Distributions",
        fontsize=13, y=1.01,
    )
    plt.tight_layout()
    path = os.path.join(output_dir, "histograms.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Histograms saved → {path}")


# ---------------------------------------------------------------------------
# Console report
# ---------------------------------------------------------------------------

_METRIC_LABELS = [
    ("avg_wait_b1",    "Avg Waiting Time – Buffer 1 (ms)"),
    ("avg_wait_b2",    "Avg Waiting Time – Buffer 2 (ms)"),
    ("avg_wait_n2",    "Avg Waiting Time – Node 2 Buffer (ms)"),
    ("avg_sojourn_b1", "Avg Sojourn at Node 1 – Type I  (wait + service, ms)"),
    ("avg_sojourn_b2", "Avg Sojourn at Node 1 – Type II (wait + service, ms)"),
    ("avg_q1_wait",    "Avg # Type I Pkts Waiting in Buffer 1"),
    ("avg_q1_total",   "Avg # Type I Pkts In-System (wait + serving)"),
    ("avg_q2_wait",    "Avg # Type II Pkts Waiting in Buffer 2"),
    ("avg_q2_total",   "Avg # Type II Pkts In-System (wait + serving)"),
    ("avg_qn2",        "Avg # Packets Waiting in Node 2 Buffer"),
    ("redirect_prob",  "Packet Redirection Probability (Router 1)"),
    ("p_wait_b1",      "Fraction of Type I Pkts That Must Wait in Buffer 1"),
    ("p_wait_b2",      "Fraction of Type II Pkts That Must Wait in Buffer 2"),
]


def print_report(accum, n_reps, sim_time):
    """Print all performance metrics with 95% t-distribution CIs."""
    sep = "=" * 72
    print(f"\n{sep}")
    print("SIMULATION REPORT")
    print(f"  Replications : {n_reps}")
    print(f"  Sim horizon  : {sim_time:,.0f} ms per replication")
    print(sep)

    for key, label in _METRIC_LABELS:
        data = accum[key]
        mean, lo, hi = ci_t(data)
        print(f"\n  {label}")
        print(f"    Point estimate : {mean:10.4f}")
        print(f"    95% CI         : [{lo:.4f},  {hi:.4f}]")

    print(f"\n{sep}\n")


def print_gof(pools):
    """Print chi-square goodness-of-fit results."""
    sep = "=" * 72
    print(sep)
    print("GOODNESS-OF-FIT  (Chi-Square Tests)")
    print(sep)

    tests = [
        ("wait_b1",    "expon",   "Waiting time – Buffer 1       vs Exponential"),
        ("wait_b2",    "expon",   "Waiting time – Buffer 2       vs Exponential"),
        ("wait_n2",    "expon",   "Waiting time – Node 2 buffer  vs Exponential"),
        ("sojourn_b1", "expon",   "Sojourn at Node 1 – Type I    vs Exponential"),
        ("sojourn_b2", "uniform", "Sojourn at Node 1 – Type II   vs Uniform"),
    ]

    for key, dist, label in tests:
        data = pools[key]
        if len(data) < 20:
            print(f"\n  {label}")
            print(f"    Insufficient data ({len(data)} samples)")
            continue
        chi2, df, p = chi_square_gof(data, dist)
        verdict = "fail to reject H0 (good fit)" if p > 0.05 else "REJECT H0 (poor fit)"
        print(f"\n  {label}")
        print(f"    χ²={chi2:.3f},  df={df},  p-value={p:.4f}  →  {verdict}")

    print(f"\n{sep}\n")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    N_REPS    = 30
    SIM_TIME  = 100_000.0   # ms
    BASE_SEED = 42

    print(f"\nRunning {N_REPS} replications × {SIM_TIME:,.0f} ms …")
    accum, pools = run_replications(
        n_reps=N_REPS, sim_time=SIM_TIME, base_seed=BASE_SEED
    )

    print_report(accum, N_REPS, SIM_TIME)
    print_gof(pools)
    print("Generating histograms …")
    plot_histograms(pools)
    print("\nDone.")
