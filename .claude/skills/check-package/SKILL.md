---
name: check-package
description: Run devtools::check() and summarize ERRORs, WARNINGs, and NOTEs in a clean list
disable-model-invocation: false
---

Run `devtools::check()` on this R package using the Bash tool:

```bash
Rscript -e "devtools::check()"
```

Parse the output and report:
1. A one-line status: `✓ No issues` or `✗ N ERROR(s), N WARNING(s), N NOTE(s)`
2. A bulleted list of each ERROR, WARNING, and NOTE with its message (skip the surrounding boilerplate)
3. If there are ERRORs, suggest the most likely fix based on the message

Do not print the full raw output — just the structured summary.
