---
title: A breif introduction to MarkDown
subtitle: Master it in 5-10 minutes
author: Xijiang Yu
institute: IHA, NMBU
date: Sept. 9, 2022
theme: Warsaw
header-includes:
	- \usepackage{mhchem}
---
## What is MarkDown?
- plain text formatting syntax
  - plain text
  - easy to read and write
  - lightweight markup (e.g., HTML) language
- peripherals
  - editors
  - conversion tools
  
## Where is MarkDown
- GitHub, Slack, MediaWiki
- Teams, Emails
- Notes taking
- Books, Presentation, Documentation
- Report and test field
  - e.g., Jupyter

## What I am talking today
- the markdown syntax
  - a 5-10 min tour
- some peripherals
- some tricks \pause
- Not RMarkDown \pause

![Rmd vs md](https://d33wubrfki0l68.cloudfront.net/61d189fd9cdf955058415d3e1b28dd60e1bd7c9b/b739c/lesson-images/rmarkdownflow.png)

## Most used elements of an article
- structure, and paragraphs
- lists
- text formatting / emphasis
- tables
- figures
- quotes, paragraph
- codes, equations, hyper-link

## Headings
- Headings
```markdown
# Heading level 1
## Level 2
### Level 3
#### Level 4
##### Level 5
###### Level 6
```
- Paragraphs are separated by new lines

## Lists

::: {.columns}
:::: {.column width=.5}
```markdown
- unordered
  - an ordinary item
  - [ ] item not checked
  - [x] item checked
- ordered list
  1. item one
  2. item two
```
::::
:::: {.column width=.5}
- unordered
  - item
  - [ ] item not checked
  - [x] item checked
- ordered list
  1. item one
  2. item two
::::
:::

## Text emphasis
::: {.columns}
:::: {.column width=.5}
```markdown
- *Italic*
- **Bold**
- ***Bold and italic***
- ~~strikethrough~~
```
::::
:::: {.column width=.5}
- *Italic*
- **Bold**
- ***Bold and italic***
- ~~strikethrough~~
::::
:::

<!--
::: {.columns}
:::: {.column width=.5}
```markdown
```
::::
:::: {.column width=.5}
::::
:::
-->
## Quotes
::: {.columns}
:::: {.column width=.5}
```markdown
William Shakespeare:

> To be, or not to be, 
> that is the question:
> Whether 'tis nobler
> in the mind to suffer
```
::::
:::: {.column width=.5}
William Shakespeare:

> To be, or not to be, 
> that is the question:
> Whether 'tis nobler
> in the mind to suffer
::::
:::

## Display colorful codes
::: {.columns}
:::: {.column width=.5}
```markdown
opening three{`}julia
a = rand(3, 4) # 3x4 matrix
c = a'a        # 4x4 matrix
closing three{`}
```
::::
:::: {.column width=.5}
```julia
a = rand(3, 4) # 3×4 matrix
c = a'a        # 4×4 matrix
```
::::
:::

## Equations
- You type

```markdown
$e^{i\pi} + 1 = 0$, 
$\begin{bmatrix}1 & 2\\3 & 4\end{bmatrix}$
$\begin{pmatrix}1 & 2\\3 & 4\end{pmatrix}$
```
- You get

$e^{i\pi} + 1 = 0$, 
$\begin{bmatrix}1 & 2\\3 & 4\end{bmatrix}$,
$\begin{pmatrix}1 & 2\\3 & 4\end{pmatrix}$

## Tables
::: {.columns}
:::: {.column width=.5}
```markdown
| col 1   | col 2 | col 3   |
|--------:|:-----:|:--------|
| r- adj. | cntrd | l- adj. |
| more    | more  | more    |
|  ...    | ...   | ...     |
```
::::
:::: {.column width=.5}

\hspace{3cm}

| col 1 | col 2  | col 3   |
|--------:|:------:|:-------|
| r- adj. | cntrd | l- adj. |
|    more | more  | more    |
|     ... | ...   | ...     |
::::
:::

## Figures
::: {.columns}
:::: {.column width=.5}
```markdown
![A picture from NMBU](https://www.nmbu.no/
sites/default/files/styles/
banner_landscape/public/nmbu-054462.jpg)
```
::::
:::: {.column width=.5}

\hspace{3cm}

![A picture from NMBU](https://www.nmbu.no/sites/default/files/styles/banner_landscape/public/nmbu-054462.jpg)
::::
:::
 
## Some tricks
- revision
  - `<!-- -->`
  - git
- figure of desired size.
  - `![Joplin](joplin.png){width=30%}`

![NMBU](https://www.nmbu.no/sites/default/files/styles/banner_landscape/public/nmbu-054462.jpg?itok=ShjhhfkM){width=30%}

## Some tricks - c.
- YAML header
```yaml
---
title: "An introduction to MarkDown"
subtitle: "With many examples"
author: "Xijiang Yu"
institute: "NMBU"
date: "29/03/2022"
abstract: "Something you learn quick, may never forget"
keywords: 
  - key1
  - key2
tags:
  - tag1
  - tag2
---
```
## Using rich HTML to markup
```html
<p align="center"> fig, table, or text
</p>
```

- Floating object
```html
<div "margin-bottom: 1rem; overflow-x autl;">
...
</div>
```
<!--
## $\LaTeX$ Formula are rich
```latex
$\begin{bmatrix}...\end{bmatrix}$
$\begin{pmatrix}...\end{pmatrix}$
```
-->
<!-- https://ashki23.github.io/markdown-latex.html-->

## Peripherals
- editors
  - emacs, vscode, vim, nano, notepad, jupyter-lab/notebook
  - any text editor
- converters
  - e.g., pandoc
- notebooks
  - Joplin and many more

## Using MarkDown to keep diaries, notes and report

$\mathrm{Diary}\to\mathrm{Year}\to{Month}\to\mathrm{yyyy-mm-dd.md}$

```bash
grep -rnw '/path/to/' -e 'pattern'
# e.g.,
grep -rnw . -e '^##'
# will list all files with 2nd+ level titles, also the line number
grep -Ril '^##' . # will only show file name
```
