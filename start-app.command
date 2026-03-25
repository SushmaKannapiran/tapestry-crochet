#!/bin/bash
cd "$(dirname "$0")"

PORT=8080
echo ""
echo "  ╔══════════════════════════════════════════╗"
echo "  ║   Tapestry Crochet Pattern Converter     ║"
echo "  ║──────────────────────────────────────────║"
echo "  ║   App running at:                        ║"
echo "  ║   http://localhost:$PORT                  ║"
echo "  ║                                          ║"
echo "  ║   To install as app:                     ║"
echo "  ║   Chrome → ⋮ → Install App               ║"
echo "  ║   Safari → Share → Add to Home Screen    ║"
echo "  ║                                          ║"
echo "  ║   Press Ctrl+C to stop                   ║"
echo "  ╚══════════════════════════════════════════╝"
echo ""

open "http://localhost:$PORT/tapestry-converter.html" 2>/dev/null || xdg-open "http://localhost:$PORT/tapestry-converter.html" 2>/dev/null

python3 -m http.server $PORT --bind 0.0.0.0 2>/dev/null || python -m SimpleHTTPServer $PORT
