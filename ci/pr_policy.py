#!/usr/bin/env python3
"""Pull-request policy check for FEBio-Rust.

Fails (exit 1) unless:
  1. the PR title and every commit subject follow Conventional Commits
     (type(scope)!: description), merge commits excepted;
  2. there is no AI attribution in the PR title, PR description, commit messages or
     commit author/committer identities;
  3. VERSION is bumped exactly as the PR type requires
     (breaking -> major, feat -> minor, anything else -> patch);
  4. CHANGELOG.md has a "## [<new version>]" entry;
  5. the Rust crate version (rust/febio_rs_solver/Cargo.toml) equals VERSION.

In GitHub Actions the PR is read from $GITHUB_EVENT_PATH. Locally:
  ci/pr_policy.py --title "perf(x): ..." --body-file body.md --base origin/main --head HEAD
"""
import argparse, json, os, re, subprocess, sys

TYPES = "feat|fix|perf|refactor|docs|test|build|ci|chore|revert|style"
CONVENTIONAL = re.compile(rf"^(?P<type>{TYPES})(\([a-z0-9._/-]+\))?(?P<bang>!)?: \S.*$")
AI = re.compile(
    r"(?i)\b(claude|anthropic|chatgpt|openai|copilot|gemini)\b"
    r"|generated\s+(with|by)\s+\[?(ai|claude|chatgpt|copilot|an? ai)"
    r"|co-authored-by:[^\n]*(\bai\b|bot|claude|anthropic|copilot|openai)"
    r"|noreply@anthropic\.com|\U0001F916")
SEMVER = re.compile(r"^(\d+)\.(\d+)\.(\d+)$")

def git(*args):
    return subprocess.run(["git", *args], capture_output=True, text=True, check=True).stdout

def show(rev, path):
    r = subprocess.run(["git", "show", f"{rev}:{path}"], capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else None

def bump(v, level):
    ma, mi, pa = v
    return (ma + 1, 0, 0) if level == "major" else (ma, mi + 1, 0) if level == "minor" else (ma, mi, pa + 1)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--title"); ap.add_argument("--body-file"); ap.add_argument("--base"); ap.add_argument("--head")
    a = ap.parse_args()
    title, body, base, head = a.title, "", a.base, a.head
    if a.body_file:
        body = open(a.body_file).read()
    ev = os.environ.get("GITHUB_EVENT_PATH")
    if ev and not a.title:
        pr = json.load(open(ev))["pull_request"]
        title, body = pr["title"], pr.get("body") or ""
        base, head = pr["base"]["sha"], pr["head"]["sha"]
    errors, notes = [], []

    # 1. conventional title and commit subjects
    m = CONVENTIONAL.match(title or "")
    if not m:
        errors.append(f"PR title is not a conventional commit (`type(scope): description`, types: {TYPES}): {title!r}")
    log = git("log", "--format=%H%x1f%an <%ae>%x1f%cn <%ce>%x1f%B%x1e", f"{base}..{head}")
    commits = [c.strip("\n").split("\x1f") for c in log.split("\x1e") if c.strip()]
    for sha, author, committer, msg in commits:
        subject = msg.splitlines()[0] if msg else ""
        if not subject.startswith("Merge ") and not CONVENTIONAL.match(subject):
            errors.append(f"commit {sha[:8]} subject is not a conventional commit: {subject!r}")

    # 2. no AI attribution
    for where, text in [("PR title", title or ""), ("PR description", body)] + \
            [(f"commit {c[0][:8]} message", c[3]) for c in commits] + \
            [(f"commit {c[0][:8]} author/committer", c[1] + " " + c[2]) for c in commits]:
        hit = AI.search(text)
        if hit:
            errors.append(f"AI attribution in {where}: {hit.group(0)!r}")

    # 3. version bump according to the PR type
    breaking = bool(m and m.group("bang")) or "BREAKING CHANGE:" in body or any("BREAKING CHANGE:" in c[3] for c in commits)
    level = "major" if breaking else ("minor" if m and m.group("type") == "feat" else "patch")
    new = (show(head, "VERSION") or "").strip()
    old = (show(base, "VERSION") or "").strip()
    if not SEMVER.match(new):
        errors.append(f"VERSION missing or not MAJOR.MINOR.PATCH at the PR head: {new!r}")
    elif old:
        expect = ".".join(map(str, bump(tuple(map(int, SEMVER.match(old).groups())), level)))
        if new != expect:
            errors.append(f"VERSION must be bumped {old} -> {expect} ({level}, from the PR type), found {new}")
        else:
            notes.append(f"VERSION {old} -> {new} ({level})")
    else:
        notes.append(f"VERSION introduced: {new}")

    # 4. changelog entry
    changelog = show(head, "CHANGELOG.md") or ""
    if new and f"## [{new}]" not in changelog:
        errors.append(f"CHANGELOG.md has no '## [{new}]' entry")

    # 5. Rust crate version in sync
    cargo = show(head, "rust/febio_rs_solver/Cargo.toml") or ""
    cm = re.search(r'^version\s*=\s*"([^"]+)"', cargo, re.M)
    if cm and cm.group(1) != new:
        errors.append(f"rust/febio_rs_solver/Cargo.toml version {cm.group(1)} != VERSION {new}")

    report = ["### PR policy", ""] + [f"- ✅ {n}" for n in notes] + [f"- ❌ {e}" for e in errors]
    if not errors:
        report.append(f"- ✅ conventional title and {len(commits)} conventional commit(s), no AI attribution, changelog entry present")
    print("\n".join(report))
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        open(os.environ["GITHUB_STEP_SUMMARY"], "a").write("\n".join(report) + "\n")
    return 1 if errors else 0

if __name__ == "__main__":
    sys.exit(main())
