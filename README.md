# nurbs-kernel

A **modern C++20 header-only NURBS library** that faithfully replicates all engineering algorithms from *The NURBS Book* (Piegl & Tiller, 2nd Ed., 1997).

## Goals

- Complete coverage of all algorithms in *The NURBS Book* Chapters 4–12
- Zero-dependency, header-only C++20
- Compile-time safety via `std::concepts` and `requires`
- Production-quality unit tests with Catch2

## Project Status

**Phase 0 — Architecting.** See [任务矩阵](./docs/TASK_MATRIX.md) for full algorithm breakdown.

## Quick Start

```bash
git clone https://github.com/mishuchu/nurbs-kernel.git
cd nurbs-kernel
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```

## Documentation

| Document | Contents |
|----------|----------|
| [TASK_MATRIX.md](./docs/TASK_MATRIX.md) | All algorithms, chapter by chapter |
| [MILESTONES.md](./docs/MILESTONES.md) | Three-phase delivery plan |

## Team

| Role | Agent | Status |
|------|-------|--------|
| Team Lead | @hermes/team-lead | Active |
| Architect | @hermes/architect | Pending |
| Core Engineer | @hermes/core-engineer | Pending |
| Test Engineer | @hermes/test-engineer | Pending |

## License

Apache 2.0
