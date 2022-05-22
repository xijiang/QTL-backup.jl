---
title: Advanced usage of MarkDown
author: Xijiang Yu
Institute: IHA, Norwegian University of Life Sciences
date: April 10, 2022
theme: Warsaw
subject: MarkDown tutorial
keywords: MarkDown, PDF, beamer
---

# YAML (YAML Ain't Markup Language)

- For example
```yaml
---
title: Advanced usage of MarkDown
author: Xijiang Yu
Institute: IHA, Norwegian University of Life Sciences
date: April 10, 2022
theme: Warsaw
subject: MarkDown tutorial
keywords: MarkDown, PDF, beamer
---
```

<!-- Use command pdfinfo to view pdf meta data -->

# Using rich HTML to markup

- Comments, for example
  - Can be used to trace modification
```html
<!-- This is a comment, which will not be displayed 
in target document.
There were some modifications of this MarkDown.
The reason was blah, and blah.
-->
```
\pause

- Ultimately, `git` can be used to keep a history of the document

# A model to keep your research diary
- Suppose you are using a Linux (alike) system, Linux or WSL
```bash
mkdir -p ~/Documents/diary/2022
cd ~/Documents/diary/2022
# edit 2022-04-10.md as today's diary
# put some one-line tags in the file, 
# e.g., #cattle2020, #pig2019, #salmon 2021
# you build up your diary day by day, and later
cd ~/Documents/diary
# to find all your files mentioned this project
grep -Ril '#cattle2020' .
```

# Writing a paper
- The lazy way
```bash
pandoc paper.md -o paper.docx
```

- The all MarkDown way
  - Build up your references `.bib`
  - Modify the YAML header
  - write and convert
  - See *`paper.md`* along with this presentation
  
# Writing a book
  
<!--

These are too advanced here.

\pause
- For example, centering objects

```html
<p align="center"> figure, table, or text
</p>
```

-  Floating objects
```html
<div "margin-bottom: 1rem; overflow-x autl;">
...
</div>
```
-->
