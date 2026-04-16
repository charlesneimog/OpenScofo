#!/usr/bin/env bash
set -e

# Re-exec with bash when launched as "sh installer.sh".
if [[ -z "${BASH_VERSION:-}" ]]; then
  exec bash "$0" "$@"
fi

SCRIPT_NAME="$(basename "$0")"
INSTALL_USER=""
INSTALL_HOME=""

# Try to relaunch in a terminal when started from GUI (double-click behavior varies by desktop).
if [[ ! -t 0 && -z "${OPEN_SCOFO_TERM_RELAUNCHED:-}" ]]; then
  export OPEN_SCOFO_TERM_RELAUNCHED=1
  if command -v x-terminal-emulator >/dev/null 2>&1; then
    exec x-terminal-emulator -e bash "$0"
  elif command -v gnome-terminal >/dev/null 2>&1; then
    exec gnome-terminal -- bash -lc '"$1"' _ "$0"
  elif command -v konsole >/dev/null 2>&1; then
    exec konsole -e bash "$0"
  elif command -v xfce4-terminal >/dev/null 2>&1; then
    exec xfce4-terminal --command="bash '$0'"
  elif command -v xterm >/dev/null 2>&1; then
    exec xterm -e bash "$0"
  fi
fi

log() {
  printf '[OpenScofo installer] %s\n' "$*"
}

die() {
  printf '[OpenScofo installer] ERROR: %s\n' "$*" >&2
  exit 1
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

expand_home_path() {
  local input="$1"
  if [[ "$input" == "\$HOME"* ]]; then
    printf '%s\n' "$INSTALL_HOME${input#\$HOME}"
  else
    printf '%s\n' "$input"
  fi
}

resolve_install_user() {
  if [[ "$(id -u)" -eq 0 ]]; then
    if [[ -n "${SUDO_USER:-}" && "${SUDO_USER}" != "root" ]]; then
      INSTALL_USER="$SUDO_USER"
    elif [[ -n "${PKEXEC_UID:-}" ]]; then
      INSTALL_USER="$(id -nu "$PKEXEC_UID" 2>/dev/null || true)"
    elif command -v logname >/dev/null 2>&1; then
      INSTALL_USER="$(logname 2>/dev/null || true)"
    fi

    if [[ -z "$INSTALL_USER" || "$INSTALL_USER" == "root" ]]; then
      die "Do not run as root directly. Run as your normal user, or use: sudo -E bash $SCRIPT_NAME"
    fi
  else
    INSTALL_USER="$(id -un)"
  fi

  INSTALL_HOME="$(getent passwd "$INSTALL_USER" | cut -d: -f6)"
  [[ -n "$INSTALL_HOME" ]] || die "Could not resolve home directory for user '$INSTALL_USER'"
}

run_as_install_user() {
  if [[ "$(id -u)" -eq 0 && "$INSTALL_USER" != "root" ]]; then
    if command -v sudo >/dev/null 2>&1; then
      sudo -u "$INSTALL_USER" "$@"
    elif command -v runuser >/dev/null 2>&1; then
      runuser -u "$INSTALL_USER" -- "$@"
    else
      die "Need sudo or runuser to execute install steps as '$INSTALL_USER'"
    fi
  else
    "$@"
  fi
}

needs_sudo_for_dest() {
  local dest="$1"
  case "$dest" in
    /usr/*|/opt/*|/etc/*|/lib/*|/lib64/*|/sbin/*|/bin/*|/var/*)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

safe_install_component() {
  local component="$1"
  local src_dir="$2"
  local dest_dir="$3"

  [[ -d "$src_dir" ]] || die "Missing payload directory for component '$component': $src_dir"

  if needs_sudo_for_dest "$dest_dir"; then
    log "Installing '$component' to '$dest_dir' (system path)"
    if [[ "$(id -u)" -eq 0 ]]; then
      mkdir -p "$dest_dir"
      cp -r "$src_dir/." "$dest_dir/"
    else
      sudo mkdir -p "$dest_dir"
      sudo cp -r "$src_dir/." "$dest_dir/"
    fi
  else
    log "Installing '$component' to '$dest_dir' as user '$INSTALL_USER'"
    if [[ "$(id -u)" -eq 0 && "$INSTALL_USER" != "root" ]]; then
      # When running under sudo, user directories may contain root-owned leftovers.
      # Copy as root, then fix ownership so runtime tools work for the target user.
      mkdir -p "$dest_dir"
      cp -r "$src_dir/." "$dest_dir/"
      chown -R "$INSTALL_USER":"$INSTALL_USER" "$dest_dir"
    else
      run_as_install_user mkdir -p "$dest_dir"
      run_as_install_user cp -r "$src_dir/." "$dest_dir/"
    fi
  fi
}

build_selection() {
  local raw="$1"
  local normalized
  normalized="$(printf '%s' "$raw" | tr '[:upper:]' '[:lower:]')"

  if [[ "$normalized" =~ ^[[:space:]]*a([[:space:]]*|ll[[:space:]]*)$ ]]; then
    printf 'csound\npuredata\nsupercollider\nvamp\n'
    return 0
  fi

  printf '%s\n' "$normalized" \
    | tr ',;' '\n' \
    | tr -s '[:space:]' '\n' \
    | sed '/^$/d' \
    | while IFS= read -r token; do
        case "$token" in
          1|csound) echo "csound" ;;
          2|pd|puredata) echo "puredata" ;;
          3|sc|supercollider) echo "supercollider" ;;
          4|vamp) echo "vamp" ;;
          *) die "Unknown selection token: '$token'" ;;
        esac
      done \
    | awk '!seen[$0]++'
}

main() {
  require_cmd tail
  require_cmd tar
  require_cmd cp
  require_cmd mktemp
  require_cmd getent

  resolve_install_user
  log "Using install user: $INSTALL_USER"

  local marker_line
  marker_line="$(grep -an '^__ARCHIVE_BELOW__$' "$0" | head -n1 | cut -d: -f1)"
  [[ -n "$marker_line" ]] || die "Archive marker not found in $SCRIPT_NAME"

  local archive_line=$((marker_line + 1))
  local temp_root
  temp_root="$(mktemp -d)"
  trap 'rm -rf "$temp_root"' EXIT

  log "Extracting embedded payload..."
  tail -n +"$archive_line" "$0" | tar -xz -C "$temp_root"

  local payload_root="$temp_root/payload"
  [[ -d "$payload_root" ]] || die "Extracted payload is missing 'payload/' root"

  cat <<'EOF'
Select components to install:
1) csound
2) puredata
3) supercollider
4) vamp
a) all
EOF

  printf 'Selection (example: 1; 3; 4): '
  local selection
  read -r selection
  [[ -n "$selection" ]] || die "No selection provided"

  mapfile -t chosen < <(build_selection "$selection")
  [[ ${#chosen[@]} -gt 0 ]] || die "No valid components selected"

  local csound_dest="$INSTALL_HOME/.local/lib/csound/6.0/plugins64"
  local puredata_dest="$INSTALL_HOME/Documents/Pd/externals/openscofo~"
  local supercollider_dest="$INSTALL_HOME/.local/share/SuperCollider/Extensions/OpenScofo"
  local vamp_dest="$INSTALL_HOME/.vamp"

  local component
  for component in "${chosen[@]}"; do
    case "$component" in
      csound)
        safe_install_component "$component" "$payload_root/csound" "$(expand_home_path "$csound_dest")"
        ;;
      puredata)
        safe_install_component "$component" "$payload_root/puredata" "$(expand_home_path "$puredata_dest")"
        ;;
      supercollider)
        safe_install_component "$component" "$payload_root/supercollider" "$(expand_home_path "$supercollider_dest")"
        ;;
      vamp)
        safe_install_component "$component" "$payload_root/vamp" "$(expand_home_path "$vamp_dest")"
        ;;
      *)
        die "Internal error: unsupported component '$component'"
        ;;
    esac
  done

  log "Installation complete."
}

main "$@"
exit 0

__ARCHIVE_BELOW__
