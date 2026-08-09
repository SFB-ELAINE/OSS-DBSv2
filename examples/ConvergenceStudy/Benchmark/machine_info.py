"""Collect the machine, OS and software context of a benchmark run.

Everything recorded here is something that plausibly moves the runtime:
CPU model and core count, memory, OS and kernel, the NGSolve/Netgen/NEURON
versions, and the threading environment. Two records are only comparable if
these agree, so they are stored alongside every timing.
"""

import importlib.metadata
import os
import platform
import re
import shutil
import subprocess
import sys

# Threading variables that change how much of the machine NGSolve and the
# underlying BLAS actually use. An unset value is recorded as None, which is
# not the same as "1" -- it means the library picked its own default.
# Import name -> distribution name(s), where they differ. netgen exposes no
# __version__ and ships as the "netgen-mesher" distribution.
DISTRIBUTION_NAMES = {
    "netgen": ("netgen-mesher",),
    "neuron": ("NEURON",),
}

THREAD_ENV_VARS = [
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NGSOLVE_NUM_THREADS",
]


def _run(cmd):
    """Return stripped stdout of a command, or None if it is unavailable."""
    if shutil.which(cmd[0]) is None:
        return None
    try:
        out = subprocess.run(
            cmd, capture_output=True, text=True, timeout=10, check=False
        )
    except (subprocess.SubprocessError, OSError):
        return None
    return out.stdout.strip() or None


def cpu_model():
    """Best-effort CPU model string for the current platform."""
    system = platform.system()
    if system == "Linux":
        try:
            with open("/proc/cpuinfo") as fp:
                for line in fp:
                    if line.startswith("model name"):
                        return line.split(":", 1)[1].strip()
        except OSError:
            pass
    elif system == "Darwin":
        return _run(["sysctl", "-n", "machdep.cpu.brand_string"])
    return platform.processor() or None


def physical_cores():
    """Physical core count, or None if it cannot be determined.

    Distinguished from the logical count because hyper-threading rarely
    helps FEM assembly and solving, so the physical count is usually the
    number that explains a timing difference.
    """
    system = platform.system()
    if system == "Linux":
        text = _run(["lscpu", "-p=Core,Socket"])
        if text:
            pairs = {
                line for line in text.splitlines() if line and not line.startswith("#")
            }
            if pairs:
                return len(pairs)
    elif system == "Darwin":
        value = _run(["sysctl", "-n", "hw.physicalcpu"])
        if value:
            return int(value)
    return None


def total_memory_gb():
    """Total system memory in GB, or None if it cannot be determined."""
    system = platform.system()
    if system == "Linux":
        try:
            with open("/proc/meminfo") as fp:
                for line in fp:
                    if line.startswith("MemTotal:"):
                        return round(int(line.split()[1]) / 1e6, 1)
        except OSError:
            pass
    elif system == "Darwin":
        value = _run(["sysctl", "-n", "hw.memsize"])
        if value:
            return round(int(value) / 1e9, 1)
    return None


def os_description():
    """Human-readable OS name, preferring the distribution on Linux."""
    system = platform.system()
    if system == "Linux":
        try:
            with open("/etc/os-release") as fp:
                match = re.search(
                    r'^PRETTY_NAME="?([^"\n]+)"?', fp.read(), re.MULTILINE
                )
                if match:
                    return match.group(1)
        except OSError:
            pass
        return "Linux"
    if system == "Darwin":
        return f"macOS {platform.mac_ver()[0]}"
    return f"{system} {platform.release()}"


def package_versions():
    """Versions of the packages that dominate the runtime.

    Falls back to the distribution metadata, because some packages (netgen)
    expose no ``__version__`` attribute even though they are installed.
    """
    versions = {}
    for name in ("ngsolve", "netgen", "neuron", "numpy", "scipy"):
        try:
            module = __import__(name)
        except ImportError:
            versions[name] = None
            continue
        version = getattr(module, "__version__", None) or getattr(
            module, "version", None
        )
        if not isinstance(version, str):
            for dist in (name, *DISTRIBUTION_NAMES.get(name, ())):
                try:
                    version = importlib.metadata.version(dist)
                    break
                except importlib.metadata.PackageNotFoundError:
                    version = None
        versions[name] = version
    return versions


def git_commit(repo_dir):
    """Short commit of the ossdbs checkout, with a dirty marker if modified."""
    head = _run(["git", "-C", repo_dir, "rev-parse", "--short", "HEAD"])
    if head is None:
        return None
    dirty = _run(
        ["git", "-C", repo_dir, "status", "--porcelain", "--untracked-files=no"]
    )
    return f"{head}-dirty" if dirty else head


def collect(repo_dir):
    """Return the full machine/software context as a plain dict."""
    return {
        "hostname": platform.node(),
        "os": os_description(),
        "platform": platform.system(),
        "kernel": platform.release(),
        "arch": platform.machine(),
        "cpu_model": cpu_model(),
        "cores_physical": physical_cores(),
        "cores_logical": os.cpu_count(),
        "memory_gb": total_memory_gb(),
        "python": sys.version.split()[0],
        "packages": package_versions(),
        "thread_env": {var: os.environ.get(var) for var in THREAD_ENV_VARS},
        "ossdbs_commit": git_commit(repo_dir),
    }
