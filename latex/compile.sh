#!/bin/bash
latex --shell-escape membrane.tex
latex --shell-escape membrane.tex
dvipdf membrane.dvi
rm membrane.aux membrane.auxlock  membrane.dvi membrane.log membraneNotes.bib 
