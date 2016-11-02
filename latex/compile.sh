#!/bin/bash
latex --shell-escape document.tex
latex --shell-escape document.tex
dvipdf document.dvi
rm document.aux document.auxlock  document.dvi document.log documentNotes.bib 
