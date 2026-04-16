#!/usr/bin/env python3
"""Build PyInstaller bundles for OSS-DBS entry points.

Builds one EXE per entry point from OSS-DBS, sharing a single dependency
folder via PyInstaller --onedir mode.

Each entry point is built as a separate --onedir target (no MERGE — it is
broken in PyInstaller >= 6.x).  The resulting per-entry folders are then
merged into one shared bundle so that DLLs and packages are not duplicated.

Usage:
    pip install pyinstaller
    python build_merged.py [--debug]

    --debug     Enable PyInstaller --log-level=DEBUG and print the
                generated .spec file contents

Output:
    dist/ossdbs_bundle/
        ossdbs.exe
        leaddbs2ossdbs.exe
        prepareaxonmodel.exe
        run_pathway_activation.exe
        (shared .dll/.so files, packages, etc.)
"""

import glob
import importlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

DEBUG = "--debug" in sys.argv

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

ENTRY_POINTS = {
    "ossdbs": "ossdbs/main.py",
    "leaddbs2ossdbs": "leaddbsinterface/convert_input_dictionary.py",
    "prepareaxonmodel": "leaddbsinterface/allocate_axons.py",
    "run_pathway_activation": "ossdbs/axon_processing/main.py",
}

HIDDEN_IMPORTS = [
    "ngsolve",
    "neuron",
    "scipy",
    "nibabel",
    "h5py",
    "dipy",
    "pandas",
    "matplotlib",
    "numpy",
]

COLLECT_ALL = [
    "ossdbs",
    "ngsolve",
    "neuron",
    "dipy",
    "matplotlib",
    "scipy",
    "nibabel",
    "pandas",
    "netgen",
]

ENV_PREFIX = Path(sys.prefix)

ADD_BINARIES_PATTERNS = [
    (str(ENV_PREFIX / "bin" / "*.dll"), "."),
    (str(ENV_PREFIX / "Library" / "bin" / "*.dll"), "."),
    (str(ENV_PREFIX / "Scripts" / "*.dll"), "."),
]

ADD_DATA = []

BUNDLE_NAME = "ossdbs_bundle"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _safe_path(p):
    return str(p).replace("\\", "/")


def _log_ok(msg):
    print(f"  [OK]      {msg}")


def _log_warn(msg):
    print(f"  [WARNING] {msg}")


def _log_err(msg):
    print(f"  [ERROR]   {msg}", file=sys.stderr)


def _log_info(msg):
    print(f"  [INFO]    {msg}")


# ---------------------------------------------------------------------------
# Pre-flight checks
# ---------------------------------------------------------------------------


def _check_pyinstaller():
    try:
        import PyInstaller

        _log_ok(f"PyInstaller {PyInstaller.__version__} found")
        return True
    except ImportError:
        _log_err("PyInstaller is not installed. Run: pip install pyinstaller")
        return False


def _check_entry_points():
    ok = True
    for name, path in ENTRY_POINTS.items():
        if Path(path).exists():
            _log_ok(f"{name:30s} -> {path}")
        else:
            _log_err(f"{name:30s} -> {path}  FILE NOT FOUND")
            ok = False
    return ok


def _check_main_functions():
    for name, path in ENTRY_POINTS.items():
        p = Path(path)
        if p.exists():
            source = p.read_text(encoding="utf-8", errors="replace")
            if "def main(" in source or "def main()" in source:
                _log_ok(f"{name:30s} -> main() found")
            elif "__name__" in source and "__main__" in source:
                _log_warn(f"{name:30s} -> no 'def main()' but has __main__ guard")
            else:
                _log_warn(f"{name:30s} -> no main() or __main__ guard found")


def _check_binary_patterns():
    any_binaries = False
    for pattern, _dest in ADD_BINARIES_PATTERNS:
        matches = glob.glob(pattern)
        if matches:
            _log_ok(f"'{pattern}' -> {len(matches)} file(s)")
            any_binaries = True
        else:
            _log_info(f"'{pattern}' -> no files (skipping)")
    if not any_binaries:
        _log_warn("No binary patterns matched.")


def _check_importable(label, modules):
    print(f"Checking {label}...")
    for mod in modules:
        try:
            importlib.import_module(mod)
            _log_ok(f"{mod}")
        except ImportError as e:
            _log_warn(f"{mod} -> not importable: {e}")


def _preflight():
    print("\n=== Pre-flight checks ===\n")
    print(f"Python env: {ENV_PREFIX}")

    print("Checking PyInstaller...")
    ok = _check_pyinstaller()

    print("Checking entry point source files...")
    ok = _check_entry_points() and ok

    print("Checking main() functions in entry points...")
    _check_main_functions()

    print("Checking --add-binary patterns...")
    _check_binary_patterns()

    _check_importable("hidden imports", HIDDEN_IMPORTS)
    _check_importable("collect_all packages", COLLECT_ALL)

    print()
    if ok:
        print("All critical checks passed.\n")
    else:
        print("FATAL errors detected.\n")
    return ok


# ---------------------------------------------------------------------------
# Resolve binary patterns into individual files
# Returns list of (abs_path, dest_dir) — 2-tuples with forward slashes
# ---------------------------------------------------------------------------


def _resolve_binaries():
    result = []
    for pattern, dest in ADD_BINARIES_PATTERNS:
        for m in glob.glob(pattern):
            result.append((_safe_path(os.path.abspath(m)), dest))
    return result


# ---------------------------------------------------------------------------
# Generate .spec file
#
# KEY INSIGHT about PyInstaller's tuple formats:
#   - Analysis(binaries=...) expects 2-tuples: (src_path, dest_dir)
#   - After Analysis runs, a.binaries contains 3-tuples: (dest, src, typecode)
#   - COLLECT expects 3-tuples: (dest, src, typecode)
#   - collect_all() returns 2-tuples for datas/binaries (src, dest)
#     which are the right format to APPEND to Analysis inputs but
#     NOT the right format for COLLECT
#
# So: we feed 2-tuples into Analysis, let PyInstaller convert them,
#     then only fix any leftover 2-tuples right before COLLECT.
# ---------------------------------------------------------------------------


def _generate_spec_for(name, script_path, resolved_bins, debug):
    L = []
    script_path_safe = _safe_path(script_path)
    debug_flag = "True" if debug else "False"

    L.append("# -*- mode: python ; coding: utf-8 -*-")
    L.append(f"# Auto-generated for entry point: {name}")
    L.append("")
    L.append("from PyInstaller.utils.hooks import collect_all")
    L.append("")

    # TOC sanitizer for COLLECT — ensures all entries are 3-tuples
    L.append("def _fix_toc(toc):")
    L.append("    fixed = []")
    L.append("    for entry in toc:")
    L.append("        n = len(entry)")
    L.append("        if n >= 3:")
    L.append("            fixed.append(entry[:3])")
    L.append("        elif n == 2:")
    L.append("            fixed.append((entry[0], entry[1], 'DATA'))")
    L.append("        elif n == 1:")
    L.append("            fixed.append((entry[0], entry[0], 'DATA'))")
    L.append("    return fixed")
    L.append("")

    # collect_all
    for pkg in COLLECT_ALL:
        var = pkg.replace("-", "_").replace(".", "_")
        L.append(f"{var}_d, {var}_b, {var}_h = collect_all('{pkg}')")
    if COLLECT_ALL:
        L.append("")

    # binaries as 2-tuples for Analysis (src_path, dest_dir)
    L.append("extra_binaries = [")
    for src, dest in resolved_bins:
        L.append(f"    ('{src}', '{dest}'),")
    L.append("]")
    L.append("")

    # datas
    if ADD_DATA:
        data_items = ", ".join(f"('{_safe_path(s)}', '{d}')" for s, d in ADD_DATA)
        datas = f"[{data_items}]"
    else:
        datas = "[]"

    # hidden imports
    hi_items = ", ".join(f"'{h}'" for h in HIDDEN_IMPORTS)

    # collect_all returns 2-tuples too — concat them for Analysis input
    ca_bins = " + ".join(
        f"{pkg.replace('-', '_').replace('.', '_')}_b" for pkg in COLLECT_ALL
    )
    ca_datas = " + ".join(
        f"{pkg.replace('-', '_').replace('.', '_')}_d" for pkg in COLLECT_ALL
    )
    ca_hi = " + ".join(
        f"{pkg.replace('-', '_').replace('.', '_')}_h" for pkg in COLLECT_ALL
    )

    # Analysis — all inputs are 2-tuples here
    L.append("a = Analysis(")
    L.append(f"    ['{script_path_safe}'],")
    L.append("    pathex=['.'],")
    L.append(f"    binaries=extra_binaries + {ca_bins},")
    L.append(f"    datas={datas} + {ca_datas},")
    L.append(f"    hiddenimports=[{hi_items}] + {ca_hi},")
    L.append("    hookspath=[],")
    L.append("    hooksconfig={},")
    L.append("    runtime_hooks=[],")
    L.append("    excludes=[],")
    L.append("    noarchive=False,")
    L.append(")")
    L.append("")

    # PYZ
    L.append("pyz = PYZ(a.pure)")
    L.append("")

    # EXE
    L.append("exe = EXE(")
    L.append("    pyz,")
    L.append("    a.scripts,")
    L.append("    [],")
    L.append("    exclude_binaries=True,")
    L.append(f"    name='{name}',")
    L.append(f"    debug={debug_flag},")
    L.append("    bootloader_ignore_signals=False,")
    L.append("    strip=False,")
    L.append("    upx=True,")
    L.append("    console=True,")
    L.append(")")
    L.append("")

    # COLLECT — sanitize to guarantee 3-tuples
    L.append("coll = COLLECT(")
    L.append("    exe,")
    L.append("    _fix_toc(a.binaries),")
    L.append("    _fix_toc(a.datas),")
    L.append("    strip=False,")
    L.append("    upx=True,")
    L.append("    upx_exclude=[],")
    L.append(f"    name='{name}',")
    L.append(")")
    L.append("")

    return "\n".join(L)


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------


def _build_entry(name, spec_content):
    spec_path = Path(f"{name}.spec")
    spec_path.write_text(spec_content, encoding="utf-8")

    try:
        compile(spec_content, str(spec_path), "exec")
    except SyntaxError as e:
        _log_err(f"Generated spec has a syntax error: {e}")
        _log_err("First 15 lines:")
        for i, line in enumerate(spec_content.splitlines()[:15], 1):
            _log_err(f"  {i:3d} | {line}")
        return False

    if DEBUG:
        print(f"\n--- {spec_path} ---")
        print(spec_content)
        print("--- end ---\n")

    cmd = [sys.executable, "-m", "PyInstaller", "--noconfirm"]
    if DEBUG:
        cmd.append("--log-level=DEBUG")
    cmd.append(str(spec_path))

    _log_info(f"Command: {' '.join(cmd)}")
    result = subprocess.run(cmd)

    if result.returncode != 0:
        _log_err(f"Build FAILED for '{name}' (exit code {result.returncode})")
        _log_err(f"Inspect spec: {spec_path.resolve()}")
        _log_err("First 15 lines of generated spec:")
        for i, line in enumerate(spec_content.splitlines()[:15], 1):
            _log_err(f"  {i:3d} | {line}")
        return False

    out_dir = Path("dist") / name
    if not out_dir.exists():
        _log_err(f"Expected output not found: {out_dir}")
        return False

    exe_files = [
        f for f in out_dir.glob(f"{name}*") if f.is_file() and f.suffix != ".manifest"
    ]
    if not exe_files:
        _log_err(f"No executable found in {out_dir}")
        return False

    for ef in exe_files:
        size_mb = ef.stat().st_size / (1024 * 1024)
        _log_ok(f"Built {ef.name} ({size_mb:.1f} MB)")
    return True


# ---------------------------------------------------------------------------
# Merge
# ---------------------------------------------------------------------------


def _copy_tree_items(src_dir, bundle, errors, entry_name=None):
    for item in src_dir.iterdir():
        dest = bundle / item.name
        try:
            if entry_name is not None:
                # Merging: only copy exe or missing files
                is_new = item.name.startswith(entry_name) or not dest.exists()
                if item.is_file() and is_new:
                    shutil.copy2(item, dest)
                elif item.is_dir() and not dest.exists():
                    shutil.copytree(item, dest)
            elif item.is_dir():
                shutil.copytree(item, dest)
            else:
                shutil.copy2(item, dest)
        except Exception as e:
            errors.append(f"Copy failed: {item} -> {dest}: {e}")


def _merge_bundles():
    dist = Path("dist")
    bundle = dist / BUNDLE_NAME
    names = list(ENTRY_POINTS.keys())

    print(f"\n=== Merging into {bundle} ===\n")

    if bundle.exists():
        shutil.rmtree(bundle)
    bundle.mkdir(parents=True)

    errors = []

    base = dist / names[0]
    if not base.exists():
        _log_err(f"Base directory not found: {base}")
        return False

    _log_info(f"Copying base from {names[0]}...")
    _copy_tree_items(base, bundle, errors)

    for name in names[1:]:
        src_dir = dist / name
        if not src_dir.exists():
            _log_err(f"Directory not found: {src_dir}")
            errors.append(f"Missing: {src_dir}")
            continue
        _log_info(f"Merging {name}...")
        _copy_tree_items(src_dir, bundle, errors, entry_name=name)

    if errors:
        for e in errors:
            _log_err(e)
        return False
    _log_ok("Merge complete")
    return True


# ---------------------------------------------------------------------------
# Verify
# ---------------------------------------------------------------------------


def _verify_bundle():
    bundle = Path("dist") / BUNDLE_NAME
    names = list(ENTRY_POINTS.keys())
    print("\n=== Verifying bundle ===\n")

    all_found = True
    for name in names:
        matches = [
            m
            for m in bundle.glob(f"{name}*")
            if m.is_file() and m.suffix != ".manifest"
        ]
        if matches:
            for ef in matches:
                size_mb = ef.stat().st_size / (1024 * 1024)
                _log_ok(f"{ef.name:40s} ({size_mb:.1f} MB)")
        else:
            _log_err(f"'{name}' executable NOT FOUND in bundle")
            all_found = False

    total = sum(f.stat().st_size for f in bundle.rglob("*") if f.is_file())
    total_mb = total / (1024 * 1024)
    file_count = sum(1 for _ in bundle.rglob("*") if _.is_file())
    print(f"\n  Bundle:     {bundle.resolve()}")
    print(f"  Total size: {total_mb:.1f} MB across {file_count} files")
    if not all_found:
        _log_err("Some executables are missing.")
    return all_found


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    """Build all entry points and merge into a single bundle."""
    print("=" * 60)
    print("  OSS-DBS PyInstaller Build (--onedir, no MERGE)")
    print("=" * 60)

    if not _preflight():
        sys.exit(1)

    resolved = _resolve_binaries()
    _log_info(f"Resolved {len(resolved)} individual binary files to bundle")

    names = list(ENTRY_POINTS.keys())
    failed = []

    for name in names:
        print(f"\n{'=' * 60}")
        print(f"  Building: {name}")
        print(f"{'=' * 60}\n")
        spec = _generate_spec_for(name, ENTRY_POINTS[name], resolved, DEBUG)
        if not _build_entry(name, spec):
            failed.append(name)
            _log_err(f"'{name}' failed — continuing with remaining entry points\n")

    if failed:
        print(f"\n{'=' * 60}")
        _log_err(f"FAILED: {', '.join(failed)}")
        _log_err("")
        _log_err("Re-run with --debug for full detail")
        print(f"{'=' * 60}")
        sys.exit(1)

    if not _merge_bundles():
        _log_err("Bundle merge failed.")
        sys.exit(1)

    if not _verify_bundle():
        sys.exit(1)

    print(f"\n{'=' * 60}")
    print("  BUILD SUCCESSFUL")
    print(f"{'=' * 60}\n")


if __name__ == "__main__":
    main()
