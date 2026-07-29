#!/usr/bin/env bash
set -u

FILES=("cyclopentadiene" "benzene"  "pyrazine" #"pyridine"
  "chlorobenzene" "dichlorobenzene" "hexafluorobenzene" #"hexachlorobenzene"
  )
BUILDDIR="out"
DUMP="dump.txt"

# Größeren TeX-Eingabepuffer setzen
export buf_size=1000000

mkdir -p "$BUILDDIR"

: > "$DUMP"

OK=()
FAILED=()

for MAIN in "${FILES[@]}"; do
    echo "Compiling $MAIN.tex ..."

    TEXFILE="$MAIN.tex"
    TMPJOB="${MAIN}_build"

    FINALPDF="$BUILDDIR/$MAIN.pdf"
    TMPPDF="$BUILDDIR/$TMPJOB.pdf"
    BACKUPPDF="$BUILDDIR/$MAIN.pdf.bak.$$"

    {
        echo
        echo "========================================"
        echo "Compiling $MAIN.tex"
        echo "========================================"
    } >> "$DUMP"

    if [[ ! -f "$TEXFILE" ]]; then
        echo "ERROR: $TEXFILE does not exist. Skipping." >> "$DUMP"
        echo "  FAILED: $TEXFILE missing"
        FAILED+=("$MAIN")
        continue
    fi

    if [[ -f "$FINALPDF" ]]; then
        mv "$FINALPDF" "$BACKUPPDF"
    fi

    rm -f "$BUILDDIR/$TMPJOB".aux \
          "$BUILDDIR/$TMPJOB".log \
          "$BUILDDIR/$TMPJOB".out \
          "$BUILDDIR/$TMPJOB".toc \
          "$TMPPDF"

    if {
        pdflatex \
            -interaction=nonstopmode \
            -halt-on-error \
            -output-directory="$BUILDDIR" \
            -jobname="$TMPJOB" \
            "$TEXFILE" &&

        pdflatex \
            -interaction=nonstopmode \
            -halt-on-error \
            -output-directory="$BUILDDIR" \
            -jobname="$TMPJOB" \
            "$TEXFILE"
    } >> "$DUMP" 2>&1
    then
        mv "$TMPPDF" "$FINALPDF"
        rm -f "$BACKUPPDF"

        echo "  OK: $FINALPDF"
        OK+=("$MAIN")
    else
        echo "ERROR: Compilation failed for $MAIN.tex" >> "$DUMP"
        echo "  FAILED: $MAIN.tex"

        rm -f "$TMPPDF"

        if [[ -f "$BACKUPPDF" ]]; then
            mv "$BACKUPPDF" "$FINALPDF"
            echo "Restored previous $FINALPDF" >> "$DUMP"
        fi

        FAILED+=("$MAIN")
        continue
    fi
done

echo
echo "Summary:"
echo "--------"

if [[ ${#OK[@]} -gt 0 ]]; then
    echo "Successful:"
    printf '  %s\n' "${OK[@]}"
fi

if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo
    echo "Failed or missing:"
    printf '  %s\n' "${FAILED[@]}"
    echo
    echo "Details are in $DUMP"
    exit 1
else
    echo "All files compiled successfully."
    rm -f "$DUMP"
    rm $BUILDDIR/*aux
    rm $BUILDDIR/*log
fi
