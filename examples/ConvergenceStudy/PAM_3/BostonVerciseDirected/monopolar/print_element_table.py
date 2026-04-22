"""Print a table of element counts per convergence strategy from existing logs."""

import os
import re

script_dir = os.path.dirname(os.path.abspath(__file__))

results = []
for entry in sorted(os.listdir(script_dir)):
    if not entry.startswith("Results_PAM_"):
        continue
    log_path = os.path.join(script_dir, entry, "ossdbs.log")
    if not os.path.isfile(log_path):
        continue

    name = entry.replace("Results_PAM_", "")
    before = after = None
    hp_applied = False

    with open(log_path) as f:
        for line in f:
            m = re.search(r"before material refinement:(\d+)", line)
            if m:
                before = int(m.group(1))
            m = re.search(r"after material refinement:(\d+)", line)
            if m:
                after = int(m.group(1))
            if "Applying HP" in line:
                hp_applied = True

    if before is not None and after is not None:
        results.append((name, before, after, hp_applied))

# Print table
hdr = (
    f"{'Strategy':<45s}  {'Initial':>10s}  {'After Mat.':>10s}  {'Delta':>12s}  {'HP'}"
)
print(hdr)
print("-" * len(hdr))
for name, before, after, hp in results:
    delta = f"+{after - before}" if after != before else ""
    hp_str = "yes" if hp else ""
    print(f"{name:<45s}  {before:>10d}  {after:>10d}  {delta:>12s}  {hp_str}")
