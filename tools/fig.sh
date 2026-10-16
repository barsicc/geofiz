#!/bin/bash
# Сохранить страницы PDF как рисунки: tools/fig.sh <код> <pdf> <стр> [<стр> ...]
# → рисунки/<код>/стрNNN.jpg  (110 dpi, качество 80; для палеток/планшетов задай DPI=150)
code=$1; pdf=$2; shift 2; dpi=${DPI:-110}
mkdir -p "рисунки/$code"
for p in "$@"; do
  pdftoppm -r "$dpi" -jpeg -jpegopt quality=80 -f "$p" -l "$p" -singlefile "$pdf" "рисунки/$code/стр$(printf %03d "$p")"
done
