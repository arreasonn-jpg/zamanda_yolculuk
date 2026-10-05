# Results policy

* `validated/` — outputs of stages that passed their acceptance tests
  (regression artifacts; regenerable via the matching runner).
* `exploratory/` — in-progress runs, parameter scans without acceptance
  criteria. Promote or delete; never cite as evidence.
* `archived/` — historical outputs, INCLUDING FAILURES. A FAIL result is
  valuable data and is never deleted.

File names keep their original `ztXXX...` identifiers so the archive stays
traceable to `stages/`.
