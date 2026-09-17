#!/usr/bin/env bash
# Mirror a private source repository into this public repository, then commit and push.
# Only files committed at the source repository's HEAD are copied.
set -euo pipefail

TARGET_REPO=$(git -C "$(dirname "${BASH_SOURCE[0]}")/.." rev-parse --show-toplevel)
SOURCE_REPO=${SOURCE_REPO:-"$(dirname "$TARGET_REPO")/mtdna_raw"}
COMMIT_MESSAGE=""
DRY_RUN=false
ALLOW_EXISTING_CHANGES=false

usage() {
    cat <<'EOF'
Usage: bash scripts/sync-private-source.sh [options]

Mirror the committed HEAD of a private source repository into this public
repository, create one commit, and push main to origin.

Options:
  --source PATH              Private source repository (default: ../mtdna_raw)
  --message MESSAGE          Commit message (default includes source HEAD)
  --dry-run                  Show the rsync changes without modifying or pushing
  --allow-existing-changes   Include existing public working-tree changes in the commit
  --help                     Show this help

Environment:
  SOURCE_REPO                Same as --source
EOF
}

while (($#)); do
    case "$1" in
        --source)
            SOURCE_REPO=${2:?"--source requires a path"}
            shift 2
            ;;
        --message)
            COMMIT_MESSAGE=${2:?"--message requires text"}
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --allow-existing-changes)
            ALLOW_EXISTING_CHANGES=true
            shift
            ;;
        --help)
            usage
            exit 0
            ;;
        *)
            printf 'Unknown option: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

SOURCE_REPO=$(git -C "$SOURCE_REPO" rev-parse --show-toplevel)
if [[ "$SOURCE_REPO" == "$TARGET_REPO" ]]; then
    printf 'Source and target repositories must be different.\n' >&2
    exit 2
fi

if [[ $(git -C "$TARGET_REPO" branch --show-current) != main ]]; then
    printf 'Refusing to sync: target must be on main.\n' >&2
    exit 2
fi

script_path="scripts/$(basename "${BASH_SOURCE[0]}")"
existing_changes=$(git -C "$TARGET_REPO" status --porcelain --untracked-files=all | grep -Fvx "?? $script_path" || true)
if [[ -n "$existing_changes" && "$ALLOW_EXISTING_CHANGES" != true ]]; then
    printf 'Refusing to commit existing target changes. Review or commit them first,\n' >&2
    printf 'or rerun with --allow-existing-changes to include them.\n' >&2
    exit 2
fi

source_head=$(git -C "$SOURCE_REPO" rev-parse HEAD)
source_short_head=$(git -C "$SOURCE_REPO" rev-parse --short HEAD)
temp_dir=$(mktemp -d)
trap 'rm -rf "$temp_dir"' EXIT
staging_dir="$temp_dir/source"

git clone --quiet --no-local --no-checkout "$SOURCE_REPO" "$staging_dir"
git -C "$staging_dir" checkout --quiet "$source_head"

rsync_args=(--archive --delete --exclude='.git' --exclude="/$script_path")
if [[ "$DRY_RUN" == true ]]; then
    rsync_args+=(--dry-run --itemize-changes)
fi

rsync "${rsync_args[@]}" "$staging_dir/" "$TARGET_REPO/"

if [[ "$DRY_RUN" == true ]]; then
    printf 'Dry run complete; no files were changed, committed, or pushed.\n'
    exit 0
fi

if [[ -z "$COMMIT_MESSAGE" ]]; then
    COMMIT_MESSAGE="sync: mirror mtdna_raw $source_short_head"
fi

git -C "$TARGET_REPO" add -A
if git -C "$TARGET_REPO" diff --cached --quiet; then
    printf 'Target already matches source HEAD %s; nothing to commit.\n' "$source_short_head"
    exit 0
fi

git -C "$TARGET_REPO" commit -m "$COMMIT_MESSAGE"
git -C "$TARGET_REPO" push origin main
printf 'Mirrored source HEAD %s, committed, and pushed main.\n' "$source_short_head"
