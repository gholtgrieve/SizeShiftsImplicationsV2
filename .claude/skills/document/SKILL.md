---
name: document
description: Run devtools::document() to regenerate NAMESPACE and man/*.Rd files, then report what changed
disable-model-invocation: false
---

Run `devtools::document()` on this R package:

```bash
Rscript -e "devtools::document()"
```

After it completes:
1. Run `git diff --stat NAMESPACE man/` to show what changed
2. Report a one-line summary: which `.Rd` files were added/modified/deleted and whether NAMESPACE changed
3. If NAMESPACE changed, show the diff so the user can verify exports look correct
4. Remind the user to stage `NAMESPACE` and `man/` before committing if they haven't already
