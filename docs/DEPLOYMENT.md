# PITAX online deployment

## Recommended path: Posit Connect Cloud

PITAX is an R/Shiny application, so GitHub Pages cannot run it by itself. The simplest GitHub-centered deployment is Posit Connect Cloud: the code stays in GitHub and Connect Cloud runs the Shiny server.

### Prepare the repository

1. Make sure PITAX runs locally.
2. From the project root run:

```r
source("scripts/prepare_connect_cloud.R")
```

This creates `manifest.json` with the R and package dependencies required by Connect Cloud.

3. Commit and push it:

```text
git add manifest.json scripts/prepare_connect_cloud.R docs/DEPLOYMENT.md
git commit -m "Prepare PITAX for Connect Cloud"
git push
```

### Publish

1. Sign in to Posit Connect Cloud.
2. Click **Publish** and choose **Shiny**.
3. Install/authorize the Posit Connect Cloud GitHub App when prompted.
4. Select repository `arielt-Gur-lab/PITAX`.
5. Select branch `main`.
6. Select `app.R` as the primary file.
7. Publish.
8. Leave automatic republish-on-push enabled if you want production to follow `main` automatically.

The deployed app receives a normal web URL and can be opened from other computers without installing R.

## Data behavior

PITAX currently treats each browser session as working state. AB1 files are uploaded by the user and results/projects are downloaded by the user. Do not rely on the hosting instance's local filesystem as permanent laboratory storage. The `.sangerproject` Save/Load workflow remains the portable persistence mechanism.

## Development workflow for v3

Keep the stable online deployment on `main`. Develop v3 on a separate branch and merge only after testing. This prevents unfinished multi-locus changes from replacing the stable online v2.x app.
