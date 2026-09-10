#!/bin/bash
# Khởi động web UI cho pipeline. Mặc định mở trong LAN, không cần đăng nhập.
#
#   bash web/start.sh              # cổng 8765
#   bash web/start.sh --port 9000
set -euo pipefail

cd "$(dirname "$0")/.."

PORT=8765
for ((i = 1; i <= $#; i++)); do
    if [[ "${!i}" == "--port" ]]; then
        j=$((i + 1))
        PORT="${!j}"
    fi
done

IP=$(hostname -I 2>/dev/null | awk '{print $1}')
echo "Trong máy này : http://localhost:${PORT}"
[ -n "$IP" ] && echo "Trong mạng LAN: http://${IP}:${PORT}"
echo

exec python3 web/server.py "$@"
