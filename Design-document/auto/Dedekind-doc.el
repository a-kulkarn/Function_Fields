(TeX-add-style-hook
 "Dedekind-doc"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("hyperref" "bookmarks" "colorlinks" "breaklinks")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "url")
   (TeX-run-style-hooks
    "latex2e"
    "amsart"
    "amsart10"
    "color"
    "verbatim"
    "rotating"
    "amssymb"
    "xypic"
    "soul"
    "hyperref")
   (TeX-add-symbols
    "tr"
    "str"
    "p"
    "P"
    "su"
    "oo"
    "NM"))
 :latex)

