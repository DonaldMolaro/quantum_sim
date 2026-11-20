#!/usr/bin/env python3
"""Run every example through the simulator and build browser visualizations.

For each example this script drives ./quantum_sim, uses the EXPORT command to capture the
circuit and the state after every step, and writes a self-contained HTML page (the viewer
with the data embedded) plus an index page that links them all.

Sources, in the order they appear on the index page:
  1. examples/*.qsim            runnable circuit files
  2. algorithm catalog (below)   one representative run of every algorithm that has a circuit
  3. examples/*.md               best-effort replay of the commands in each tutorial's code blocks

Usage: scripts/render_examples.sh [--out DIR] [--only TEXT] [--seed N] [--binary PATH]
"""
import argparse
import hashlib
import html
import json
import os
import re
import subprocess
import sys

MIN_STEPS = 3  # INIT plus a single gate is not worth a page

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# One representative run per algorithm. (name, title, commands, extra environment, expects_circuit)
CATALOG = [
    ("grover", "Grover search (3 qubits, targets 3 and 5)", ["GROVER 3 5"], {}, True),
    ("grover_auto", "Grover with auto-tuned iterations", ["GROVER AUTO 3 6 3 5"], {}, True),
    ("quantum_counting", "Quantum counting (4 qubits, 1 target)", ["QCOUNT RUN 4 8 3"], {}, True),
    ("qpe_4bit_t_gate", "Phase estimation of the T gate", ["QPE 4 PI/4"], {}, True),
    ("qpe_4bit_s_gate", "Phase estimation of the S gate", ["QPE 4 PI/2"], {}, True),
    ("shor_15", "Shor's algorithm, N=15, a=7", ["SHOR 15"], {"QSIM_SHOR_FORCE_A": "7"}, True),
    ("deutsch_jozsa_balanced", "Deutsch-Jozsa, balanced oracle", ["DEUTSCH_JOZSA 3 BALANCED_XOR0"], {}, True),
    ("deutsch_jozsa_constant", "Deutsch-Jozsa, constant oracle", ["DEUTSCH_JOZSA 3 CONST0"], {}, True),
    ("bernstein_vazirani", "Bernstein-Vazirani, secret 11", ["BV 4 11"], {}, True),
    ("qec_bit_flip", "3-qubit bit-flip error correction", ["QEC RUN 1 2"], {}, True),
    ("qrng", "Quantum random number generator", ["QRNG 4"], {}, True),
    ("qubo_grover", "QUBO solved with Grover threshold search",
     ["QUBO GROVER 3 0 1 1 0 0 0 1 0 0 0 1"], {}, True),
    ("qaoa", "QAOA on a 3-variable QUBO (final optimized circuit)",
     ["VQA QAOA 3 1 0 10 0.25 -2 0 2 0 1 0 2 0 -3"], {}, True),
    ("vqe_zz", "VQE for Z0 Z1 (final optimized ansatz)", ["VQE RUN 2 2 10 0.3 0 1 1.0 2 Z 0 Z 1"], {}, True),
    ("latin_square", "Latin square search with Grover (12 qubits)", ["LATIN DEMO"], {}, True),
    # Classical in this simulator: there is no quantum state to show.
    ("simon", "Simon's algorithm", ["SIMON 4 10 8"], {}, False),
    ("anneal", "Simulated / quantum annealing", ["ANNEAL DEMO"], {}, False),
    ("tsp", "Travelling salesman via QUBO", ["TSP DEMO"], {}, False),
]

ALGORITHM_COMMANDS = {
    "GROVER", "BV", "BERNSTEIN_VAZIRANI", "DEUTSCH_JOZSA", "DEUTSCH", "QPE", "PHASE_ESTIMATION",
    "QCOUNT", "QUANTUM_COUNTING", "SHOR", "QEC", "QRNG", "QUBO", "VQA", "QAOA", "VQE", "LATIN",
    "SIMON", "ANNEAL", "TSP",
}
OTHER_COMMANDS = {
    "INIT", "H", "X", "Y", "Z", "S", "T", "SDG", "TDG", "P", "RX", "RY", "RZ", "RU", "RESET",
    "CX", "CNOT", "CZ", "CY", "CH", "CRX", "CRY", "CRZ", "CP", "CU", "SWAP", "ISWAP",
    "XX", "YY", "ZZ", "CCX", "TOFFOLI", "CSWAP", "FREDKIN", "MCX", "MEASURE", "IF",
    "NOISE", "SEED", "QFTMODE",
}
KNOWN = ALGORITHM_COMMANDS | OTHER_COMMANDS


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def clean_command(line):
    """Return the command text if the line looks like a simulator command, else None."""
    line = line.strip()
    if line.startswith("QSIM>"):
        line = line[5:].strip()
    if not line or line.startswith("#") or "<" in line or "[" in line:
        return None  # comments and syntax summaries such as "QPE <m> <phase>"
    word = line.split()[0].upper()
    return line if word in KNOWN else None


def qsim_commands(path):
    """(command, heading) pairs; a .qsim file has no headings."""
    cmds = []
    for line in read(path).splitlines():
        if line.strip().upper().startswith("EXPORT"):
            continue  # the script does its own exporting
        c = clean_command(line)
        if c:
            cmds.append((c, None))
    return cmds


def markdown_commands(path):
    """(command, nearest preceding heading) pairs from the fenced code blocks."""
    cmds, in_block, heading = [], False, None
    for line in read(path).splitlines():
        if line.strip().startswith("```"):
            in_block = not in_block
            continue
        if in_block:
            c = clean_command(line)
            if c:
                cmds.append((c, heading))
        elif line.startswith("#"):
            heading = re.sub(r"^#+\s*", "", line).strip()
            heading = re.sub(r"^\d+\)\s*", "", heading)  # "1) Run the demo" -> "Run the demo"
    return cmds


def with_exports(commands, tmp_prefix):
    """Insert EXPORT before every INIT / algorithm command, and one at the end.

    A new INIT or algorithm run replaces the recording, so exporting just before it
    captures the circuit built so far. Duplicates are removed afterwards.
    """
    out, titles, current = [], [], None
    for c, heading in commands:
        word = c.split()[0].upper()
        if word == "INIT" or word in ALGORITHM_COMMANDS:
            out.append("EXPORT %s_%03d.json" % (tmp_prefix, len(titles)))
            titles.append(current)
            current = None
        out.append(c)
        if current is None:
            current = heading
    out.append("EXPORT %s_%03d.json" % (tmp_prefix, len(titles)))
    titles.append(current)
    return out, titles


def run_session(binary, commands, env_extra, seed, timeout):
    env = dict(os.environ)
    env.update(env_extra)
    env["QSIM_RNG_SEED"] = str(seed)
    script = "\n".join(["SEED %d" % seed] + commands + ["QUIT"]) + "\n"
    try:
        proc = subprocess.run([binary], input=script, capture_output=True, text=True,
                              cwd=ROOT, env=env, timeout=timeout)
    except subprocess.TimeoutExpired:
        return False, "timed out after %ds" % timeout
    if proc.returncode != 0:
        return False, "simulator exited with status %d" % proc.returncode
    return True, ""


def load_export(path):
    try:
        with open(path, encoding="utf-8") as f:
            doc = json.load(f)
    except (OSError, ValueError):
        return None
    steps = doc.get("steps") or []
    return doc if len(steps) >= MIN_STEPS else None


def fingerprint(doc):
    return hashlib.sha1(json.dumps(doc["steps"], sort_keys=True).encode()).hexdigest()


def bundle(viewer_src, doc, title, out_path):
    data = json.dumps(doc, separators=(",", ":")).replace("</", "<\\/")
    page = viewer_src.replace("<title>Quantum Sim Viewer</title>",
                              "<title>%s - Quantum Sim Viewer</title>" % html.escape(title), 1)
    marker = '<script>\n"use strict";'
    if marker not in page:
        raise RuntimeError("viewer.html changed: cannot find the script marker to embed data")
    page = page.replace(marker, "<script>window.QSIM_DATA=%s;</script>\n%s" % (data, marker), 1)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(page)


INDEX_CSS = """
:root{--bg:#f6f7f9;--panel:#fff;--ink:#1c2230;--muted:#667085;--line:#d6dae1;--accent:#2f6fed}
@media (prefers-color-scheme:dark){:root{--bg:#12151c;--panel:#1b202b;--ink:#e6e9f0;--muted:#98a2b3;--line:#2f3646;--accent:#6c9bff}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.5 system-ui,sans-serif}
main{max-width:900px;margin:0 auto;padding:24px 16px}h1{font-size:20px;margin:0 0 4px}
h2{font-size:13px;text-transform:uppercase;letter-spacing:.05em;color:var(--muted);margin:28px 0 8px}
p.lead{color:var(--muted);margin:0}ul{list-style:none;margin:0;padding:0;background:var(--panel);border:1px solid var(--line);border-radius:10px}
li{display:flex;gap:12px;justify-content:space-between;padding:9px 14px;border-top:1px solid var(--line)}li:first-child{border-top:0}
a{color:var(--accent);text-decoration:none}a:hover{text-decoration:underline}.meta{color:var(--muted);white-space:nowrap;font-variant-numeric:tabular-nums}
"""


def write_index(out_dir, sections, skipped):
    parts = ['<!doctype html><html lang="en"><head><meta charset="utf-8">',
             '<meta name="viewport" content="width=device-width, initial-scale=1">',
             "<title>Quantum Sim Examples</title><style>%s</style></head><body><main>" % INDEX_CSS,
             "<h1>Quantum Sim Examples</h1>",
             '<p class="lead">Each page replays a circuit step by step: drag the slider or use the arrow keys.</p>']
    for heading, items in sections:
        if not items:
            continue
        parts.append("<h2>%s</h2><ul>" % html.escape(heading))
        for it in items:
            parts.append('<li><a href="%s">%s</a><span class="meta">%d qubits &middot; %d steps</span></li>'
                         % (html.escape(it["file"]), html.escape(it["title"]), it["qubits"], it["steps"]))
        parts.append("</ul>")
    if skipped:
        parts.append("<h2>No circuit to show</h2><ul>")
        for name, why in skipped:
            parts.append('<li><span>%s</span><span class="meta">%s</span></li>'
                         % (html.escape(name), html.escape(why)))
        parts.append("</ul>")
    parts.append("</main></body></html>")
    with open(os.path.join(out_dir, "index.html"), "w", encoding="utf-8") as f:
        f.write("".join(parts))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(ROOT, "viewer", "exports"))
    ap.add_argument("--binary", default=os.path.join(ROOT, "quantum_sim"))
    ap.add_argument("--only", default="", help="only run examples whose name contains this text")
    ap.add_argument("--seed", type=int, default=7, help="RNG seed so measurements are reproducible")
    ap.add_argument("--timeout", type=int, default=120, help="seconds allowed per example")
    args = ap.parse_args()

    if not os.access(args.binary, os.X_OK):
        sys.exit("error: %s not found; build first (see scripts/render_examples.sh)" % args.binary)
    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=True)
    tmp_dir = os.path.join(out_dir, ".tmp")
    os.makedirs(tmp_dir, exist_ok=True)
    viewer_src = read(os.path.join(ROOT, "viewer", "viewer.html"))

    jobs = []  # (section, name, title, commands, env, expects_circuit)
    for fn in sorted(os.listdir(os.path.join(ROOT, "examples"))):
        path = os.path.join(ROOT, "examples", fn)
        if fn.endswith(".qsim"):
            jobs.append(("Example circuits (examples/*.qsim)", fn[:-5], fn[:-5].replace("_", " "),
                         qsim_commands(path), {}, True))
    for name, title, cmds, env, expects in CATALOG:
        jobs.append(("Algorithms (final or representative run)", name, title,
                     [(c, None) for c in cmds], env, expects))
    for fn in sorted(os.listdir(os.path.join(ROOT, "examples"))):
        if fn.endswith(".md"):
            cmds = markdown_commands(os.path.join(ROOT, "examples", fn))
            if cmds:
                jobs.append(("Tutorial walkthroughs (commands replayed from examples/*.md)",
                             fn[:-3], fn[:-3].replace("_", " "), cmds, {}, True))

    sections, skipped = {}, []
    order = []
    used_labels = {}  # job name -> {label: times used}, to keep page titles unique
    for section, name, title, cmds, env, expects in jobs:
        if args.only and args.only not in name:
            continue
        if section not in sections:
            sections[section] = []
            order.append(section)
        prefix = os.path.join(tmp_dir, name)
        session, titles = with_exports(cmds, prefix)
        count = len(titles)
        ok, why = run_session(args.binary, session, env, args.seed, args.timeout)
        if not ok:
            print("[FAIL] %-28s %s" % (name, why))
            skipped.append((name, why))
            continue
        kept, seen = [], set()  # kept holds (doc, heading)
        for i in range(count):
            doc = load_export("%s_%03d.json" % (prefix, i))
            if doc is None:
                continue
            fp = fingerprint(doc)
            if fp not in seen:
                seen.add(fp)
                kept.append((doc, titles[i]))
        for i in range(count):
            try:
                os.remove("%s_%03d.json" % (prefix, i))
            except OSError:
                pass
        if not kept:
            reason = "classical in this simulator" if not expects else "no gate circuit recorded"
            print("[skip] %-28s %s" % (name, reason))
            skipped.append((name, reason))
            continue
        for k, (doc, heading) in enumerate(kept, 1):
            suffix = "" if len(kept) == 1 else "_%d" % k
            stem = name + suffix
            if len(kept) == 1:
                label = title
            elif heading:
                label = "%s: %s" % (title, heading)
            else:
                label = "%s (part %d of %d)" % (title, k, len(kept))
            used = used_labels.setdefault(name, {})
            used[label] = used.get(label, 0) + 1
            if used[label] > 1:
                label = "%s (%d)" % (label, used[label])
            with open(os.path.join(out_dir, stem + ".json"), "w", encoding="utf-8") as f:
                json.dump(doc, f, separators=(",", ":"))
            bundle(viewer_src, doc, label, os.path.join(out_dir, stem + ".html"))
            sections[section].append({
                "file": stem + ".html", "title": label,
                "qubits": max(s["state"]["n"] for s in doc["steps"]), "steps": len(doc["steps"]),
            })
            print("[ ok ] %-28s %d steps -> %s.html" % (stem, len(doc["steps"]), stem))

    try:
        os.rmdir(tmp_dir)
    except OSError:
        pass
    write_index(out_dir, [(h, sections[h]) for h in order], skipped)
    total = sum(len(v) for v in sections.values())
    print("\n%d visualizations written to %s" % (total, out_dir))
    print("Open %s" % os.path.join(out_dir, "index.html"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
