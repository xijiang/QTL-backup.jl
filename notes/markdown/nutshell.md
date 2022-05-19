---
geometry: margin=2cm
header-includes:
	- \usepackage{mhchem}
---

<!--
In this example, I used [a lua fiter](https://github.com/jdutant/columns) 
to create this two-column PDF output.
The command is 
- pandoc nutshell.md --lua-filter="path-to/columns.lua" -o nutshell.pdf
-->

::: {.columns}
:::: {.column width=.5}
```markdown
# MarkDown in a nutshell

- A five-minute crush course

## Introduction

MarkDown is a set of lightweight markup, e.g.,
HTML, language syntax expressed in plain text.
In this MarkDown cheatsheet, you will learn 
most often used elements of MarkDown.

## The syntex

### Paragraph

MarkDown paragraphs are separated by empty
line(s).  Hence, though this paragraph is 
of many lines, it is still one paragraph.
Headings are led by `#`. The number of `#`
indicates the level.

### Text formats

You can type *italic*, **bold**, ***bold and 
italic***, and ~~strikethrough~~ texts.

### Lists
- an unordered list
  - item one
  - [ ] an item not checked
  - [x] a checked item
- an ordered list
  1. item one
  1. item two

### Quotes

William Shakespeare:

> To be, or not to be, 
> that is the question:
> Whether 'tis nobler
> in the mind to suffer

### Including source codes

three-back-quote{`}julia
a = rand(3, 4) # -> 3×4 matrix
c = a'a        # -> 4×4 matrix
another-three-back-quote{`}

### Equations

$e^{i\pi} + 1 = 0$, 
$\begin{bmatrix}1 & 2 \\ 3 & 4\end{bmatrix}$,
$\begin{pmatrix}1 & 2 \\ 3 & 4\end{pmatrix}$
$\ce{Na2SO4 ->[H2O] 2Na+ + SO4^2-}$
```
::::
:::: {.column width=.5}
# MarkDown in a nutshell

- A five-minute crush course

## Introduction

MarkDown is a set of lightweight markup, e.g.,
HTML, language syntax expressed in plain text.
In this MarkDown cheatsheet, you will learn 
most often used elements of MarkDown.

## The syntex

### Paragraph

MarkDown paragraphs are separated by empty
line(s).  Hence, though this paragraph is 
of many lines, it is still one paragraph.
Headings are led by `#`. The number of `#`
indicates the level.

### Text formats

You can type *italic*, **bold**, ***bold and 
italic***, and ~~strikethrough~~ texts.

### Lists
- an unordered list
  - item one
  - [ ] an item not checked
  - [x] a checked item
- an ordered list
  1. item one
  1. item two

### Quotes

William Shakespeare:

> To be, or not to be, 
> that is the question:
> Whether 'tis nobler
> in the mind to suffer

### Including source codes
```julia
a = rand(3, 4) # -> 3×4 matrix
c = a'a        # -> 4×4 matrix
```

### Equations

$e^{i\pi} + 1 = 0$, 
$\begin{bmatrix}1 & 2 \\ 3 & 4\end{bmatrix}$, and
$\begin{pmatrix}1 & 2 \\ 3 & 4\end{pmatrix}$

$\ce{Na2SO4 ->[H2O] 2Na+ + SO4^2-}$
::::
:::

\newpage

```markdown
### Tables

| Col 1 | 2 | 3 |
| --: | :--: | :-- |
| right adj. | centered | left adj. |
| more | more | more |
```
### Tables

|      Col 1 | 2        | 3         |
|-----------:|:--------:|:----------|
| right adj. | centered | left adj. |
|       more | more     | more      |

```markdown
### Hyperlinks
- Jeg jobber p\aa\ [biovitenskap](https://www.nmbu.no/fakultet/biovit/om) om 
  [Norges milj\o- og biovitenskapelige universitet](https://nmbu.no).
```

### Hyperlinks
- Jeg jobber p\aa\ [biovitenskap](https://www.nmbu.no/fakultet/biovit/om) om [Norges milj\o- og biovitenskapelige universitet](https://nmbu.no).

```markdown
### Figures
![A picture of NMBU](https://www.nmbu.no/sites/default/files/styles/banner_landscape/
public/nmbu-054462.jpg){width=66%}
```

### Figures
![A picture of NMBU](https://www.nmbu.no/sites/default/files/styles/banner_landscape/public/nmbu-054462.jpg){width=80%}
