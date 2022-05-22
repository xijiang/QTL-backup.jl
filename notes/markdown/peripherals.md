---
title: MarkDown peripherals
author: Xijiang Yu
Institute: IHA, Norwegian University of Life Sciences
date: April 10, 2022
theme: Warsaw
output: beamer
---

# MarkDown

- plain text formatting syntax
  - easy to write and read
  - lightweight MarkUp (e.g., HTML)
  - portable, future proof. \pause

- Widely used
  - website, documents, notes, books, presentations, email messages, technical documentation.
  - for example
    - Microsoft teams
    - GitHub
    - Jupyter
    - Wordpress 
	- Reddit \pause

- Peripherals
  - MarkDown editors
  - Conversion tools
  
# Editors

- Any **text** editors
  - Emacs, nano, vim, notepad, vscode, jupyter-notebook, -lab

- A converter
  - `pandoc`

- WYSIWYG
  - Joplin, and many more
	- Many platforms: Windows, Mac, Linux, Mobile
	- synchronization, e.g., on (private) onedrive
	- Many functions, e.g., web clips.
  - Usually two panes: text-editor + display

# `pandoc` examples

```bash
# Note, pandoc can recognize format by .ext
pandoc source.md -t beamer -o slides.pdf  # to beamer presentation
pandoc source.md -o doc.pdf # to a pdf document
pandoc source.md -o doc.docx # to windows word document
pandoc source.md -o web.html # to html format
pandoc source.md # print html codes on screen
pandoc source.md -o book.epub # to an e-book
# and many more
# pandoc can also convert other files to MarkDown
# see pandoc manual
```
