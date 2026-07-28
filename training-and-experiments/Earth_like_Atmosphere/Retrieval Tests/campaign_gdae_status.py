#!/usr/bin/env python3
"""Summarize G-DAE retrieval campaign progress."""

from __future__ import annotations

import argparse
import re
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from campaign_common import CAMPAIGN_DIR, TEST_IDS, iter_cases, normalize_test_id
from campaign_run_gdae_queue import result_exists


LOG_RE = re.compile(
    r"(?P<stamp>\d{8}_\d{6})_(?P<test_id>test_\d{2})_(?P<branch>phoenix|sphinx)_"
    r"gdae_(?P<f_spot>\d+\.\d+)spot-(?P<f_fac>\d+\.\d+)fac\.log$"
)
START_RE = re.compile(r"^# Started: (?P<value>.+)$", re.MULTILINE)
FINISH_RE = re.compile(r"^# Finished: (?P<value>.+)$", re.MULTILINE)
RETURN_RE = re.compile(r"^# Return code: (?P<value>-?\d+)$", re.MULTILINE)


@dataclass(frozen=True)
class Job:
    test_id: str
    branch: str
    f_spot: float
    f_fac: float
    model_name: str

    @property
    def key(self) -> tuple[str, str, str, str]:
        return (self.test_id, self.branch, f"{self.f_spot:.2f}", f"{self.f_fac:.2f}")


@dataclass
class LogInfo:
    path: Path
    started: datetime | None = None
    finished: datetime | None = None
    return_code: int | None = None

    @property
    def duration_seconds(self) -> float | None:
        if self.started is None or self.finished is None:
            return None
        return (self.finished - self.started).total_seconds()


def parse_dt(value: str) -> datetime | None:
    try:
        return datetime.fromisoformat(value.strip())
    except ValueError:
        return None


def parse_log(path: Path) -> LogInfo:
    text = path.read_text(encoding="utf-8", errors="replace")
    started_match = START_RE.search(text)
    finished_match = FINISH_RE.search(text)
    return_match = RETURN_RE.search(text)
    return LogInfo(
        path=path,
        started=parse_dt(started_match.group("value")) if started_match else None,
        finished=parse_dt(finished_match.group("value")) if finished_match else None,
        return_code=int(return_match.group("value")) if return_match else None,
    )


def build_jobs(test_ids: list[str]) -> list[Job]:
    jobs: list[Job] = []
    for test_id in test_ids:
        for branch in ("phoenix", "sphinx"):
            for case in iter_cases(branch):
                jobs.append(Job(test_id, branch, case.f_spot, case.f_fac, case.model_name("gdae")))
    return jobs


def latest_logs() -> dict[tuple[str, str, str, str], LogInfo]:
    logs: dict[tuple[str, str, str, str], LogInfo] = {}
    for path in sorted((CAMPAIGN_DIR / "logs").glob("*_gdae_*.log")):
        match = LOG_RE.search(path.name)
        if not match:
            continue
        key = (
            match.group("test_id"),
            match.group("branch"),
            f"{float(match.group('f_spot')):.2f}",
            f"{float(match.group('f_fac')):.2f}",
        )
        logs[key] = parse_log(path)
    return logs


def active_process_lines() -> list[str]:
    try:
        completed = subprocess.run(
            ["bash", "-lc", "ps -ef | grep -E 'campaign_run_gdae_queue|campaign_retrieval_mpi|mpirun' | grep -v grep || true"],
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError:
        return []
    return [line for line in completed.stdout.splitlines() if line.strip()]


def format_duration(seconds: float | None) -> str:
    if seconds is None:
        return "n/a"
    minutes, sec = divmod(int(round(seconds)), 60)
    hours, minutes = divmod(minutes, 60)
    if hours:
        return f"{hours}h {minutes}m"
    if minutes:
        return f"{minutes}m {sec}s"
    return f"{sec}s"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test-id", action="append", default=None)
    args = parser.parse_args()

    test_ids = [normalize_test_id(value) for value in args.test_id] if args.test_id else list(TEST_IDS)
    jobs = build_jobs(test_ids)
    logs = latest_logs()
    process_lines = active_process_lines()

    done = failed = running_or_open = pending = 0
    durations: list[float] = []
    current: list[str] = []
    failed_items: list[str] = []

    for job in jobs:
        log = logs.get(job.key)
        has_result = result_exists(job.test_id, job.branch, job.model_name)
        if log and log.return_code is not None:
            if log.return_code == 0 and has_result:
                done += 1
                duration = log.duration_seconds
                if duration is not None:
                    durations.append(duration)
            else:
                failed += 1
                failed_items.append(f"{job.test_id} {job.branch} {job.f_spot:.2f}/{job.f_fac:.2f} rc={log.return_code}")
        elif log:
            running_or_open += 1
            current.append(f"{job.test_id} {job.branch} {job.f_spot:.2f}/{job.f_fac:.2f} since {log.started or 'unknown'}")
        else:
            pending += 1

    avg = sum(durations) / len(durations) if durations else None
    remaining_after_current = max(len(jobs) - done - failed, 0)
    eta = avg * remaining_after_current if avg is not None else None

    print(f"Total GDAE jobs: {len(jobs)}")
    print(f"Done: {done}")
    print(f"Running/open: {running_or_open}")
    print(f"Failed: {failed}")
    print(f"Pending/no log yet: {pending}")
    print(f"Average completed duration: {format_duration(avg)}")
    print(f"ETA from completed average: {format_duration(eta)}")
    print(f"Active process lines: {len(process_lines)}")
    if current:
        print("Current/open jobs:")
        for item in current[:5]:
            print(f"  - {item}")
    if failed_items:
        print("Failed jobs:")
        for item in failed_items:
            print(f"  - {item}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
